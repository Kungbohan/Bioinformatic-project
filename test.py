#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:40:40 2019

@author: kung
"""
import numpy as np
a=np.array([[0,0,0,0,0],
            [5,0,0,0,0],
            [9,10,0,0,0],
            [9,10,8,0,0],
            [8,9,7,3,0]])

name=['A','B','C','D','E']
import Neighbor_Joining
tree,final_distance=Neighbor_Joining.Neighbor_joining(a, name)

from ete3 import Tree
unrooted_tree = Tree(tree+';')
print (final_distance)

