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

def readfile(filename):
    read={}
    file = open(filename,"r")
    for line in file:
      if line.startswith('@G') :
        read[line]=next(file)
    return read

print readfile("HB173_sn12.fastq")