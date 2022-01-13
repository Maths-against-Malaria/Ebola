# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:13:59 2020

@author: helle
"""

"""
# code based on code of ??? ADE group.
This Code  Build a map which associate to a population name of the model  an
index, in order to create a list that will record the dynamics of the populations.

We first created a dictionary named "index" which will has as keys the names of
populations and as  values the associated index.

output desired ==> index : dictionary
"""
#from parameters_reinfection import *
#from parameters_original import NE, NP, NI, NL
import numpy as np
#from parameters_2 import *

# Name: compartments
# first index: _, t: tilde; s: star
# second index: _, p: home; h: hospital; i: isolation; 
#              j: after safe funeral; f: after normal funeral

def indexFunction(n):
    compartments = [('S',  0, '__'),
                ('E', n, '__'), ('E',n, 't_'),('E', n, 's_'),
                ('P', n, '__'), ('P',n, 't_'),('P', n, 's_'),
                ('I', n, '_p'), ('I',n, '_h'),('I', n, '_i'), ('I', n, 'sh'), ('I',n, 'sp'),
                ('F',  0, '__'), 
                ('B',  0, '_j'),('B',  0, '_f'),
                ('R',  0, '__')
                ]
    index = dict()
    ind = 0
    for i in compartments:
        notation = i[0] + i[2]
        if i[1] == 0:
            index[notation] = ind
            ind += 1
        else:
            for k in np.arange(1, i[1] + 1):
                notation_k = notation + str(k)
                index[notation_k] = ind
                ind += 1
    #print(index)
    return index

