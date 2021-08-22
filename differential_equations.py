# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 12:35:02 2020

@author: helle
"""

#from parameters_2 import *
#from index import *
from functions import *
import numpy as np

def dS(pop, lt, N, index, vac):
    x = - (lt/N) * pop[index['S__']] - vac
    return x
    
def dE__1(pop, la, N, FE, index):
    x = (la/N) * pop[index['S__']] \
        - FE * pop[index['E__1']]
    return x
   
def dE__k(pop, j, FE, index): # k2 <= j <= NE
    x = FE * pop[index['E__' + str(j-1)]] \
        - FE * pop[index['E__' + str(j)]]
    return x


def dEs_1(pop, ls, N, FT, FE, index):
    x = (ls/N) * pop[index['S__']] \
        - (FT + FE) * pop[index['Es_1']]
    return x

def dEs_k(pop, j, FE, FT, index): # k2 <= j <= NE
    x = FE * pop[index['Es_' + str(j-1)]] \
        - (FT + FE) * pop[index['Es_' + str(j)]]
    return x


def dEt_1(pop, FT, FE, index):
    x =   FT * pop[index['Es_1']] \
        - FE * pop[index['Et_1']]
    return x

def dEt_k(pop, j, FT, FE, index): # k2 <= j <= NE
    x =   FT * pop[index['Es_' + str(j)]] \
        + FE * pop[index['Et_' + str(j-1)]] \
        - FE * pop[index['Et_' + str(j)]]
    return x


def dP__1(pop, FE, FP, NE, index):
    x =   FE * pop[index['E__' + str(NE)]] \
        - FP * pop[index['P__1']]
    return x
   
def dP__k(pop, j, FP, index): # k2 <= j <= NP
    x =   FP * pop[index['P__' + str(j-1)]] \
        - FP * pop[index['P__' + str(j)]]
    return x


def dPs_1(pop, FE, FT, FP, NE, index):
    x =   FE * pop[index['Es_' + str(NE)]] \
        - (FT + FP) * pop[index['Ps_1']]
    return x
   
def dPs_k(pop, j, FP, FT, index): # k2 <= j <= NP
    x =   FP * pop[index['Ps_' + str(j-1)]] \
        - (FT + FP) * pop[index['Ps_' + str(j)]]
    return x


def dPt_1(pop, FT, FE, FP, NE, index):
    x =   FT * pop[index['Ps_1']] \
        + FE * pop[index['Et_' + str(NE)]] \
        - FP * pop[index['Pt_1']]
    return x
   
def dPt_k(pop, j, FT, FP, index): # k2 <= j <= NP
    x =   FT * pop[index['Ps_' + str(j)]] \
        + FP * pop[index['Pt_' + str(j-1)]] \
        - FP * pop[index['Pt_' + str(j)]]
    return x


def dI_p1(pop, fp, FP, FI, NP, index):
    x =    fp * FP * pop[index['P__' + str(NP)]] \
              - FI * pop[index['I_p1']]
    return x

def dI_pk(pop, j, FI, index):
    x =   FI * pop[index['I_p' + str(j-1)]] \
        - FI * pop[index['I_p' + str(j)]]
    return x    

def dI_h1(pop, fh, FP, FI, NP, index):
    x =    fh * FP * pop[index['P__' + str(NP)]] \
              - FI * pop[index['I_h1']]
    return x

def dI_hk(pop, j, FI, index):
    x =   FI * pop[index['I_h' + str(j-1)]] \
        - FI * pop[index['I_h' + str(j)]]
    return x


def dIsp1(pop, fp, FP, FT, FI, NP, index):
    x =    fp * FP * pop[index['Ps_' + str(NP)]] \
              - (FT + FI) * pop[index['Isp1']]
    return x

def dIspk(pop, j, FI, FT, index): # 2 <=j < NIs
    x =   FI * pop[index['Isp' + str(j-1)]] \
        - (FT + FI) * pop[index['Isp' + str(j)]]
    return x 

def dIspNI(pop, FI, FT, NIp, index):
    x =   FI * pop[index['Isp' + str(NIp-1)]] \
        - FT * pop[index['Isp' + str(NIp)]]
    return x 


def dIsh1(pop, fh, FP, FT, FI, NP, index):
    x =    fh * FP * pop[index['Ps_' + str(NP)]] \
              - (FT + FI) * pop[index['Ish1']]
    return x

def dIshk(pop, j, FI, FT, index): # 2 <=j < NIh
    x =   FI * pop[index['Ish' + str(j-1)]] \
        - (FT + FI) * pop[index['Ish' + str(j)]]
    return x 

def dIshNI(pop, FI, FT, NIh, index):
    x =   FI * pop[index['Ish' + str(NIh-1)]] \
        - FT * pop[index['Ish' + str(NIh)]]
    return x 

def dI_i1(pop, fi, FP, FT, FI, NP, index):
    x =   fi * FP * pop[index['P__' + str(NP)]] \
        + fi * FP * pop[index['Ps_' + str(NP)]] \
        +      FP * pop[index['Pt_' + str(NP)]] \
        +      FT * pop[index['Isp1']] \
        +      FT * pop[index['Ish1']] \
        -      FI * pop[index['I_i1']]
    return x

def dI_ik(pop, j, FI, FT, index):
    x =   FI * pop[index['I_i' + str(j-1)]] \
        + FT * pop[index['Isp' + str(j)]] \
        + FT * pop[index['Ish' + str(j)]] \
        - FI * pop[index['I_i' + str(j)]]
    return x
        

def dR(pop,fdead_p, fdead_h, fdead_i, FI, NIp, NIh, NIi, index, vac):
    x = FI * ((1-fdead_p) * pop[index['I_p' + str(NIp)]] \
              + (1-fdead_h) * pop[index['I_h' + str(NIh)]] \
              + (1-fdead_i) * pop[index['I_i' + str(NIi)]])\
            + vac
    return x

def dF(pop, fdead_p, fdead_h, FI, FF, NIp, NIh, index,d_h, d_p):
    x = FI * ((1-d_p) * fdead_p * pop[index['I_p' + str(NIp)]] \
               +  (1-d_h) * fdead_h * pop[index['I_h' + str(NIh)]]) \
                - FF  * pop[index['F__']]
    return x

def dB_f(pop, FF, index):
    x = FF *   pop[index['F__']]
    return x

def dB_j(pop, fdead_i, FI, NIi, NIh, NIp, index, d_h, d_p, fdead_h, fdead_p):
    x = FI * (fdead_i * pop[index['I_i' + str(NIi)]] \
                + d_p * fdead_p * pop[index['I_p' + str(NIp)]] \
                + d_h * fdead_h * pop[index['I_h' + str(NIh)]])
                  # * (1-q(t)) * d_h
    return x
