# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:13:23 2020

@author: helle
"""

from parameters_1 import *
from index import *
import numpy as np


#index = indexFunction([NE, NP, NIp, NIh, NIi])

def popsum(pop, compartment, script1, script2, Nerls, index):
    if compartment in {'S', 'F', 'B', 'R'}:
        Nerl = 0
    if compartment == 'E':
        Nerl = Nerls[0]
    if compartment == 'P':
        Nerl = Nerls[1]
    if compartment == 'I':
        if script2 == 'p':
            Nerl = Nerls[2]
        if script2 == 'h':
            Nerl = Nerls[3]
        if script2 == 'i':
            Nerl = Nerls[4]
    x = 0
    #print(compartment + str(Nerl))
    if Nerl == 0:
         x = x + pop[index[compartment + script1 + script2]]
    else:      
        for k in range(1, Nerl+1):
            #print(compartment + script1 + script2 + str(k))
            x = x + pop[index[compartment + script1 + script2 + str(k)]]
    return x

def q (pop, t, t_iso, qmax, Nerls, index):
    if t < t_iso:
        q = 0
    else:
        #print(index)
        Q = (popsum(pop=pop, compartment ='P', script1='t', script2='_', Nerls=Nerls, index=index) + \
             popsum(pop=pop, compartment ='I', script1='_', script2='i', Nerls=Nerls, index=index))
        #print(Q)
        if Q <= qmax:
            q = 1
        else:
            q = qmax/Q
    return q

def la(pop, fiso, f_tb, betaP, betaIp, betaIh, betaF, ph, q, Nerls, index):
    ls0 =   betaP * (fiso * \
                popsum(pop=pop, compartment='P', script1='_', script2='_', Nerls=Nerls, index=index) + \
                popsum(pop=pop, compartment='P', script1='s', script2='_', Nerls=Nerls, index=index)) + \
            betaIp * \
                popsum(pop=pop, compartment='I', script1='s', script2='p', Nerls=Nerls, index=index)  + \
            betaIh * \
                popsum(pop=pop, compartment='I', script1='s', script2='h', Nerls=Nerls, index=index) + \
            betaIh * (1-ph) * (1-q) * (
                popsum(pop=pop, compartment='P', script1='t', script2='_', Nerls=Nerls, index=index) + \
                popsum(pop=pop, compartment='I', script1='_', script2='i', Nerls=Nerls, index=index)) # infections of individuals who could be traced back
    l =     betaP * (1 - fiso) * \
                popsum(pop=pop, compartment='P', script1='_', script2='_', Nerls=Nerls, index=index) + \
            betaIp * \
                popsum(pop=pop, compartment='I', script1='_', script2='p', Nerls=Nerls, index=index) + \
            betaIh * \
                popsum(pop=pop, compartment='I', script1='_', script2='h', Nerls=Nerls, index=index) + \
            betaF * pop[index['F__']] + \
            (1-f_tb) * ls0 # infections not traced back
    ls =    f_tb * ls0 # infections that are traced back
    return [l, ls]

def la_aliou2(pop, fiso, f_tb, betaP, betaIp, betaIh, betaF, ph, q, Nerls, index):
    l =     betaP *  \
                popsum(pop=pop, compartment='P', script1='_', script2='_', Nerls=Nerls, index=index) + \
            betaIp * \
                popsum(pop=pop, compartment='I', script1='_', script2='p', Nerls=Nerls, index=index) + \
            betaIh * \
                popsum(pop=pop, compartment='I', script1='_', script2='h', Nerls=Nerls, index=index) + \
            betaF * pop[index['F__']]
    ls =    0
    return [l, ls]

def la_Aliou(pop, fiso, f_tb, betaP, betaIp, betaIh, betaF, ph, q, Nerls, index):
    ls =   f_tb * betaP * (fiso * \
                popsum(pop=pop, compartment='P', script1='_', script2='_', Nerls=Nerls, index=index) + \
                popsum(pop=pop, compartment='P', script1='s', script2='_', Nerls=Nerls, index=index)) + \
            betaIp * \
                popsum(pop=pop, compartment='I', script1='s', script2='p', Nerls=Nerls, index=index)  + \
            betaIh * \
                popsum(pop=pop, compartment='I', script1='s', script2='h', Nerls=Nerls, index=index)
    l =     (1-f_tb) * betaP * (1 - fiso) * \
                popsum(pop=pop, compartment='P', script1='_', script2='_', Nerls=Nerls, index=index) + \
            betaIp * \
                popsum(pop=pop, compartment='I', script1='_', script2='p', Nerls=Nerls, index=index) + \
            betaIh * \
                popsum(pop=pop, compartment='I', script1='_', script2='h', Nerls=Nerls, index=index) + \
            betaF * pop[index['F__']]

    return [l, ls]


def lt(la, ls):
    lt = la + ls
    return lt
'''
def f_phi(t, k, t_iso, f_p1, f_h1, f_p2, f_h2):
    if t <= t_iso:
        f_p = f_p1
        f_h = f_h1
    else:
        f_p = f_p2[k]
        f_h = f_h2[k]
    f_i = 1-(f_p + f_h)
    x = [f_p, f_h, f_i]
    return x
'''

def f_phi(t, k, t_iso, f_p1, f_h1, f_p2, f_h2):
    if t <= 100:
        f_p = f_p2[0]
        f_h = f_h2[0]
    if t <= 200:
        f_p = f_p2[1]
        f_h = f_h2[1]
    if t <= 300:
        f_p = f_p2[2]
        f_h = f_h2[2]
    if t <= 400:
        f_p = f_p2[3]
        f_h = f_h2[3]
    if t > 400:
        f_p = f_p2[4]
        f_h = f_h2[4]
    f_i = 1-(f_p + f_h)
    x = [f_p, f_h, f_i]
    return x

def vac(pop, index, t, t_vac, Nvac):
    if t <= t_vac:
        x = 0
    if t > t_vac:
        x = min(pop[index['S__']], Nvac)
    return x
'''
# test: f_phi
    print(f_phi(300, 0))
    print(f_phi(300, 1))
    print(f_phi(300, 2))
    print(f_phi(300, 3))
    print(f_phi(30, 3))


# test: popsum
pop = [2,                               # S__
       3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, # E__
       4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, # Et_
       5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, # Es_
       6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, # P__
       7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, # Pt_
       8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, # Ps_
       9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, # I_p
       10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, # I_h
       11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, # I_i
       12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12, # Ish
       13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, # Isp
       14,15,16,17] #F__, B_j, B_f, R

pop2 = [0,                               # S__
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, # E__
       4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, # Et_
       5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, # Es_
       6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, # P__
       7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, # Pt_
       8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, # Ps_
       9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, # I_p
       10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, # I_h
       11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, # I_i
       12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12, # Ish
       13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, # Isp
       14,15,16,17] #F__, B_j, B_f, R
       
NIp = NIh = NIi

# test ls
print(ls(fiso = 0.2, pop=pop))
print(betaP * (0.2 * 6*16 + 8*16) + betaIp * 13*16 + betaIh * 12*16)


# test la
print(la(fiso = 0.2, pop=pop))
print(betaP * 0.8 * 6*16 + betaIp * 9*16 + betaIh * 10*16 + betaF * 14)




# test: popsum
print(pop)
print(popsum('S', '_', '_', 0, pop))  
print(popsum('E', '_', '_', NE, pop))
print(popsum('E', 't', '_', NE, pop))
print(popsum('E', 's', '_', NE, pop))
print(popsum('P', '_', '_', NP, pop))
print(popsum('P', 's', '_', NP, pop))
print(popsum('P', 't', '_', NP, pop))
print(popsum('I', '_', 'p', NIhome, pop))
print(popsum('I', '_', 'h', NIhosp, pop))
print(popsum('I', 's', 'p', NIhome, pop))
print(popsum('I', 's', 'h', NIhosp, pop))
print(popsum('I', '_', 'i', NIiso, pop))
print(popsum('F', '_', '_', 0, pop)) 
print(popsum('B', '_', 'j', 0, pop)) 
print(popsum('B', '_', 'f', 0, pop)) 
print(popsum('R', '_', '_', 0, pop)) 

print(f_phi(200, 0))
print(f_phi(300, 0))
print(f_phi(400, 1))
print(f_phi(500, 2))
print(f_phi(600, 3))
print(f_phi(250, 0)
print(f_phi(249, 0))
'''
