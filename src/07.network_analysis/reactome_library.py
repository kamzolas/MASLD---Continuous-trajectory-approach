#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import igraph # 0.11.4


class OntoGraph(igraph.Graph):

    def __init__(self,  *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_vertex_parents(self, node):
        if type(node) == int:
            in_edges = self.incident(node, mode='in')
        else:
            node = self.vs.find(name=node)
            in_edges = self.incident(node.index, mode='in')
        parents = []
        for e in in_edges:
            parents.append(self.vs[self.es[e].source])
        return parents
    
    def get_vertex_ancestors(self, node):
        if type(node) == int:
            ancestors = self.bfs(node.index, mode='in')[0][1:]
        else:
            node = self.vs.find(name=node)
            ancestors = self.bfs(node.index, mode='in')[0][1:]
        return self.vs[ancestors]

    def get_vertex_children(self, node):
        if type(node) == int:
            out_edges = self.incident(node, mode='out')
        else:
            node = self.vs.find(name=node)
            out_edges = self.incident(node.index, mode='out')
        children = []
        for e in out_edges:
            children.append(self.vs[self.es[e].target])
        return children
    
    def get_vertex_descendants(self, node):
        if type(node) == int:
            descendants = self.bfs(node, mode='out')[0][1:]
        else:
            node = self.vs.find(name=node)
            descendants = self.bfs(node.index, mode='out')[0][1:]
        return self.vs[descendants]
    
    def get_distance_from_root(self, node):
        degrees = zip(self.vs['name'], self.degree(mode='in'))
        root_node = list(filter(lambda x: x[1]==0, degrees))[0][0]
        distance = self.distances(source=root_node, target=node)[0][0]
        return distance

    def get_root(self):
        degrees = zip(self.vs['name'], self.degree(mode='in'))
        root_node = list(filter(lambda x: x[1]==0, degrees))[0][0]
        return root_node