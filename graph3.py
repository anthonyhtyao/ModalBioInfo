# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 23:52:56 2016

@author: wanerchen
"""

import time
import random
import itertools
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import graphviz as gv
import random


def file_generator(filename, table = None) :
  """ A generator for the  read from the given file. """
  if table == None :
    table = {}
  file = open(filename,"r")
  i = 0
  while True : 
    line = next(file)
    if i == 1 :
      yield line
    i += 1 
    i = i % 4

def graph(reads,k) : 
  G=nx.MultiDiGraph()
  for s in reads : 
    l = len(s)
    for i in range(0,l-k) :
      G.add_node(s[i:i+k])
      G.add_node(s[i+1:i+k+1])
      G.add_edge(s[i:i+k],s[i+1:i+k+1])
      
  rmRoots(G)
  rmLeaves(G)  
 
  
  #errorElim(G)
  #searchBetter(G)
  return G
  
def rmLeaves(G):
    leavesList = []
    for node in G.nodes():
        if G.out_degree(node) == 0:
            leavesList.append(node)
    while(leavesList):
        leaf = leavesList.pop()
        predecessors = G.predecessors(leaf)
        G.remove_node(leaf)
        for node in predecessors:
            if G.out_degree(node) == 0:
                leavesList.append(node)
    return G
    
def rmRoots(G):
    rootList = []
    for node in G.nodes():
        if G.in_degree(node) == 0:
            rootList.append(node)
    while(rootList):
        root = rootList.pop()
        successor = G.successors(root)
        G.remove_node(root)
        for node in successor:
            if G.in_degree(node) == 0:
                rootList.append(node)
    return G
    
def clean(G):                #only combine multi-edge to one edge
    G1=nx.MultiDiGraph()
    for node in G.nodes():
        if node not in G1:
            G1.add_node(node)
        lstOfsucc=G.successors(node)
        for node1 in lstOfsucc:
            if node1 not in G1:
                G1.add_node(node1)
            G1.add_edge(node,node1)
            G1[node][node1]['weight']=G.number_of_edges(node,node1)
           
    return G1

def cleanError(G):         # 
    G1=nx.MultiDiGraph()
    for node in G.nodes():
        if node not in G1:
            G1.add_node(node)
        lstOfsucc=G.successors(node)
        if len(lstOfsucc)>1:
            maxEdge=0
            node1=lstOfsucc[0]
            for suc in lstOfsucc:
                if G.number_of_edges(node,suc)>maxEdge: 
                      maxEdge=G.number_of_edges(node,suc)
                      node1=suc
        else:
            node1=lstOfsucc[0]
        if node1 not in G1:
            G1.add_node(node1)
        G1.add_edge(node,node1,amount=G.number_of_edges(node,node1))
       
    rmRoots(G1)      
    return G1
    
def cleanMutation(G):
    visited =[]
    P = {}	# dictionary of predecessors
    dist={}
    n=0
    Q =[]	
    edge=[]
    for node in G.nodes():
        if G.out_degree(node)>1 :
            break
    root=node
    visited.append(root)
    Q.append(root)
    P[root]=root
    dist[root]=0

    while len(Q)>0:
        node=Q.pop()
        if n == G.number_of_edges():
            break
        lstOfsucc=G.successors(node)
        for succ in lstOfsucc:
            if succ in visited:
                if (node,succ) not in edge and dist[succ]> dist[node]-G[node][succ]['weight']:
                     dist[succ]=dist[node]-G[node][succ]['weight']
                     P[succ]=node
                     edge.append((node,succ))
                     n+=1
            else :
                Q.append(succ)
                visited.append(succ)
                dist[succ]=-G[node][succ]['weight']
                P[succ]=node
                edge.append((node,succ))
                n+=1
            
    G1=nx.MultiDiGraph()
    for succ,pred in P.items():
        G1.add_nodes_from([succ,pred])
        G1.add_edge(pred,succ)
    rmLeaves(G1)
    return G1
        
    
                    
                
                
        
    
            
    
    
        
     
       
        
def graphgv(reads,k) : 
  G=gv.Digraph(format='svg')
  for s in reads : 
    l = len(s)
    for i in range(0,l-k) :
      G.node(s[i:i+k],s[i:i+k])
      G.node(s[i+1:i+k+1],s[i+1:i+k+1])
      G.edge(s[i:i+k],s[i+1:i+k+1])
  return G


    
def main() :
  reads = file_generator('HB175_sn12.fastq')
  print (reads)
  k = 30
  start = time.clock()
  graph1 = graph(reads,k)
  graph2 = clean(graph1)
  graph3 = cleanMutation(graph2)
  #graph2.render(filename='img/graph2')
  #print graph2['TCCTGGCGATGATGAGTGATGGGCGAACTG']['CCTGGCGATGATGAGTGATGGGCGAACTGA']['weight']
  nx.draw(graph3)
  A= nx.nx_agraph.to_agraph(graph3)
  A.layout(prog='dot')
  A.write("3.dot")
  A.draw("3.pdf")
  plt.draw()
  print ("TIME :", time.clock() - start)
  plt.show()

if __name__ == '__main__':
  main()