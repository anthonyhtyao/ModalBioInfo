# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:05:09 2016

@author: wanerchen
"""

import networkx as nx
import graphviz as gz
import matplotlib.pyplot as plt
import random
from random import randint 

k=10
G=nx.DiGraph()
D=gz.Digraph(comment='k-mer')

def k_mer(G,D, reads, k):
    lst={}
    for read in reads:
        for i in range(len(read)-k+1):
            if i<len(read)-k:
                lst[read[i:i+k]]=(read[i:i+k],read[i+1:i+k+1])
            G.add_node(read[i:i+k])
            read[i:i+k]
        for node in G.nodes():
            if lst.has_key(node):
                 (node1,node2)=lst[node]
                 G.add_edge(node1,node2)
                 D.edge(node1,node2)
    
    return G,D

def readfile(filename):
    read=[]
    file = open(filename,"r")
    for line in file:
      if line.startswith('@G') :
        read.append(next(file)[:-1])
    return read

reads=readfile("HB173_sn12.fastq")

G,D=k_mer(G,D,reads,k)

lab = { label:label for label in G }
nx.draw(G, labels = lab)
D.render('t.gv', view =True)
plt.draw()
plt.show()