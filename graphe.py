#!/usr/bin/python -tt

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
  return G

def graphgv(reads,k) : 
  G=gv.Digraph(format='svg')
  for s in reads : 
    l = len(s)
    for i in range(0,l-k) :
      G.node(s[i:i+k],s[i:i+k])
      G.node(s[i+1:i+k+1],s[i+1:i+k+1])
      G.edge(s[i:i+k],s[i+1:i+k+1])
  return G

def error(reads) :
  random.seed(5) 
  read = []
  for s in reads : 
    l = len(s)
    n = random.randint(0, l/20)
    for i in range(0,n) : 
      j = random.randint(0,l)
      c = random.choice(['A','C','G','T'])
      s1 = s[:i]+c+s[i+1:]
      read = read+[s1]
  return read


def main() :
  reads = file_generator('snRNA12/HB173_sn12.fastq')
  print (reads)
  k = 10
  start = time.clock()
  #graph1 = graph(reads,k)
  graph2 = graphgv(reads,k)
  graph2.render(filename='img/graph2')
  #nx.draw(graph1)
  #plt.draw()
  print ("TIME :", time.clock() - start)
  #plt.show()

if __name__ == '__main__':
  main()
