library(ggplot2) # 3.4.2
library(tidyr) # 1.3.1
library(dplyr) # 1.1.4

get_average_MSE <- function(dataset) {
  global_df <- data.frame()
  for (dir in dirs) {
    key <- basename(dir)
    files <- grep(paste('^',dataset,'_squared', sep=''), x = list.files(dir), value = TRUE)
    training_residual_errors_df <- data.frame()
    # gather the MSEs for all kfolds
    for (file in files) {
      tmp_df <- read.table(paste(dir, '/', file, sep=''), sep='\t', header=TRUE)
      training_residual_errors_df <- rbind(training_residual_errors_df, tmp_df)
    }
    # calculate the average MSE per variable
    tmp_df <- as.data.frame(colMeans(training_residual_errors_df))
    colnames(tmp_df) <- key
    tmp_df <- t(tmp_df)
    global_df <- rbind(global_df, tmp_df)
  }
  global_df$deep_split <- sapply(rownames(global_df), function(x) {as.integer(unlist(strsplit(x, '_'))[1])})
  global_df$min_size <- sapply(rownames(global_df), function(x) {as.integer(unlist(strsplit(x, '_'))[2])})
  global_df$key <- rownames(global_df)
  return(global_df)
}


plot_MSE <- function(MSE_df, outputname) {
  
  MSE_df <- MSE_df[MSE_df$min_size != 15,] # 15 gives too high errors
  
  plot_df <- pivot_longer(MSE_df, cols=variables, values_to="MSE", names_to="variable")
  plot_df <- as.data.frame(plot_df)
  my_palette <- c('#82c9b6', '#c98295', '#c9b682')
  
  p <- ggplot(data=plot_df, aes(x=min_size, y=MSE, group=deep_split, color=factor(deep_split))) +
    facet_wrap(~variable, ncol = 2, scales='free_y') +
    geom_point(shape=19, size=2) +
    geom_line(size=0.75) +
    scale_x_continuous(breaks=sort(unique(plot_df$min_size)), 
                       labels=sort(unique(plot_df$min_size))) +
    #scale_y_continuous(limits = c(0.4, NA)) +
    scale_color_manual(values = my_palette, name="Deep Split") +
    labs(x = "Module Min Size", y = "Mean Squared Error") +
    theme(strip.text.x = element_text(size=14),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8),
          axis.title = element_text(size=14),
          legend.text=element_text(size=12),
          legend.title = element_text(size=14)
    ) + 
    guides(color = guide_legend(override.aes = list(size=6, linetype = c(0,0,0))))
  ggsave(paste(main_dir, outputname, '.png', sep=''), p, device="png", dpi=600, 
         height=4.5, width=8, units=c("in"))
}



main_dir <- '../../results/wgcna_and_linear_modelling/grid_params/'
dirs <- list.dirs(main_dir, recursive = FALSE)
variables <- c('Fibrosis', 'Steatosis', 'Ballooning', 'NAS')


training_MSE_df <- get_average_MSE('training')
testing_MSE_df <- get_average_MSE('testing')
plot_MSE(training_MSE_df, 'training_errors')
plot_MSE(testing_MSE_df, 'testing_errors')


training_MSE_df$training_error <- rowMeans(training_MSE_df[,variables])
training_MSE_df[,variables] <- NULL
testing_MSE_df$testing_error <- rowMeans(testing_MSE_df[,variables])
testing_MSE_df[,variables] <- NULL


testing_MSE_df[,c('deep_split', 'min_size')] <- NULL
plot_df <- merge(training_MSE_df, testing_MSE_df, by = 'key')
plot_df <- pivot_longer(plot_df, cols=c('training_error', 'testing_error'), 
                        values_to="MSE", names_to="error")
plot_df <- as.data.frame(plot_df)
my_palette <- c('#82c9b6', '#c98295', '#c9b682')


labeller_function <- as_labeller(
  c('training_error' = 'Train Error', 
    'testing_error' = 'Test Error')
)


p <- ggplot(data=plot_df, aes(x=min_size, y=MSE, group=deep_split, color=factor(deep_split))) +
  facet_wrap(~factor(error, levels=c('training_error', 'testing_error')), nrow = 2, 
             scales='free_y', labeller = labeller_function) +
  geom_point(shape=19, size=2) +
  geom_line(size=0.75) +
  scale_x_continuous(breaks=sort(unique(plot_df$min_size)), 
                     labels=sort(unique(plot_df$min_size))) +
  #scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = my_palette, name="Deep Split") +
  labs(x = "Module Min Size", y = "Average MSE") +
  theme(strip.text.x = element_text(size=14),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title = element_text(size=14),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14)
  ) + 
  guides(color = guide_legend(override.aes = list(size=6, linetype = c(0,0,0))))

ggsave(paste(main_dir, 'average_errors_plot.png', sep=''), p, device="png", 
       dpi=600, height=5, width=8, units=c("in"))


idx <- which(testing_MSE_df$testing_error == min(testing_MSE_df$testing_error))
print(testing_MSE_df[idx,])
