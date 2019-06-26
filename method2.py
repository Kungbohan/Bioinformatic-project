#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 21:58:06 2019

@author: kung
"""

from alignment.sequence import Sequence
from alignment.vocabulary import Vocabulary
from alignment.sequencealigner import SimpleScoring, GlobalSequenceAligner
import numpy as np
import pandas as pd
# Create sequences to be aligned.

def split_sequence(s):
    ss_new=[]
    for i in range(len(s)):
        ss_new.append(s[i])
    return ss_new

sequence_family=pd.read_csv('amylase.txt')
sequence_family=np.array(sequence_family)
name=[]
spilt_pos=[]
for i in range( len(sequence_family)):
    if sequence_family[i][0][0]=='[':
        name.append(sequence_family[i][0][1:-1])
        spilt_pos.append(i)
sequence=[]
for i in spilt_pos:
    ss=sequence_family[i+1][0]
    for ii in range(i+2,i+9):
        ss=ss+sequence_family[ii][0]
    sequence.append(ss)
        
#%%
v = Vocabulary()
sequence_encoded=[]
for i in range(len(sequence)):
    sequence_encoded.append(v.encodeSequence(Sequence(split_sequence(sequence[i]))))

scoring = SimpleScoring(2, -1)
aligner = GlobalSequenceAligner(scoring, -2)

Matrix=np.zeros(9*9).reshape(9,9)
for i in range(len(sequence_encoded)):
    for j in range(i+1,len(sequence_encoded)):
        score, encodeds = aligner.align(sequence_encoded[i], 
                                        sequence_encoded[j], 
                                        backtrace=True)
        for encoded in encodeds:
            alignment = v.decodeSequenceAlignment(encoded)
            score=np.floor((100-alignment.percentIdentity())*
                        len(np.array(alignment))/100)
            print(i,j,score)           
        Matrix[i,j]=score

Matrix=Matrix.T    
#np.save('M.npy',Matrix)
#%%
#Matrix=np.load('M.npy',allow_pickle=True)
#%%
import Neighbor_Joining
tree,distance=Neighbor_Joining.Neighbor_joining(Matrix, name)

from ete3 import Tree
unrooted_tree = Tree(tree+';')
print (unrooted_tree)

#%%
To_plot=1
if To_plot==1:
    from ete3 import Tree, TreeStyle
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.branch_vertical_margin = 10 
    unrooted_tree.show(tree_style=ts)


