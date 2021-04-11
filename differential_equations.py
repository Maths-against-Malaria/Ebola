# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 12:35:02 2020

@author: helle
"""

from parameters_1 import *
from index import *
from functions import *
import numpy as np

#lt_ = lt(0.2, pop)
#la_ = la(0.2, pop)
#ls_ = ls(0.2, pop)

#fp = 0.4            # fhome
#fh = 0.3            # fhosp
#fi = 1 - fp - fh    # fiso

def dS(pop, lt, N, index, vac):
    x = - (lt/N) * pop[index['S__']] - vac
    return x
    
def dE__1(pop, la, N, DE, index):
    x = (la/N) * pop[index['S__']] \
        - (1/DE) * pop[index['E__1']]
    return x
   
def dE__k(pop, j, DE, index): # k2 <= j <= NE
    x = (1/DE) * pop[index['E__' + str(j-1)]] \
        - (1/DE) * pop[index['E__' + str(j)]]
    return x


def dEs_1(pop, ls, N, DT, DE, index):
    x = (ls/N) * pop[index['S__']] \
        - (1/DT + 1/DE) * pop[index['Es_1']]
    return x

def dEs_k(pop, j, DE, DT, index): # k2 <= j <= NE
    x = (1/DE) * pop[index['Es_' + str(j-1)]] \
        - (1/DT + 1/DE) * pop[index['Es_' + str(j)]]
    return x


def dEt_1(pop, DT, DE, index):
    x =   (1/DT) * pop[index['Es_1']] \
        - (1/DE) * pop[index['Et_1']]
    return x

def dEt_k(pop, j, DT, DE, index): # k2 <= j <= NE
    x =   (1/DT) * pop[index['Es_' + str(j)]] \
        + (1/DE) * pop[index['Et_' + str(j-1)]] \
        - (1/DE) * pop[index['Et_' + str(j)]]
    return x



def dP__1(pop, DE, DP, NE, index):
    x =   (1/DE) * pop[index['E__' + str(NE)]] \
        - (1/DP) * pop[index['P__1']]
    return x
   
def dP__k(pop, j, DP, index): # k2 <= j <= NP
    x =   (1/DP) * pop[index['P__' + str(j-1)]] \
        - (1/DP) * pop[index['P__' + str(j)]]
    return x


def dPs_1(pop, DE, DT, DP, NE, index):
    x =   (1/DE) * pop[index['Es_' + str(NE)]] \
        - (1/DT + 1/DP) * pop[index['Ps_1']]
    return x
   
def dPs_k(pop, j, DP, DT, index): # k2 <= j <= NP
    x =   (1/DP) * pop[index['Ps_' + str(j-1)]] \
        - (1/DT + 1/DP) * pop[index['Ps_' + str(j)]]
    return x


def dPt_1(pop, DT, DE, DP, NE, index):
    x =   (1/DT) * pop[index['Ps_1']] \
        + (1/DE) * pop[index['Et_' + str(NE)]] \
        - (1/DP) * pop[index['Pt_1']]
    return x
   
def dPt_k(pop, j, DT, DP, index): # k2 <= j <= NP
    x =   (1/DT) * pop[index['Ps_' + str(j)]] \
        + (1/DP) * pop[index['Pt_' + str(j-1)]] \
        - (1/DP) * pop[index['Pt_' + str(j)]]
    return x



def dI_p1(pop, fp, DP, DI, NP, index):
    x =    fp * (1/DP) * pop[index['P__' + str(NP)]] \
              - (1/DI) * pop[index['I_p1']]
    return x

def dI_pk(pop, j, DI, index):
    x =   (1/DI) * pop[index['I_p' + str(j-1)]] \
        - (1/DI) * pop[index['I_p' + str(j)]]
    return x    

 
def dI_h1(pop, fh, DP, DI, NP, index):
    x =    fh * (1/DP) * pop[index['P__' + str(NP)]] \
              - (1/DI) * pop[index['I_h1']]
    return x

def dI_hk(pop, j, DI, index):
    x =   (1/DI) * pop[index['I_h' + str(j-1)]] \
        - (1/DI) * pop[index['I_h' + str(j)]]
    return x


def dIsp1(pop, fp, DP, DT, NP, index):
    x =    fp * (1/DP) * pop[index['Ps_' + str(NP)]] \
              - (1/DT + 1/DI) * pop[index['Isp1']]
    return x

def dIspk(pop, j, DI, DT, index): # 2 <=j < NIs
    x =   (1/DI) * pop[index['Isp' + str(j-1)]] \
        - (1/DT + 1/DI) * pop[index['Isp' + str(j)]]
    return x 

def dIspNI(pop, DI, DT, NIp, index):
    x =   (1/DI) * pop[index['Isp' + str(NIp-1)]] \
        - (1/DT) * pop[index['Isp' + str(NIp)]]
    return x 


def dIsh1(pop, fh, DP, DT, DI, NP, index):
    x =    fh * (1/DP) * pop[index['Ps_' + str(NP)]] \
              - (1/DT + 1/DI) * pop[index['Ish1']]
    return x

def dIshk(pop, j, DI, DT, index): # 2 <=j < NIh
    x =   (1/DI) * pop[index['Ish' + str(j-1)]] \
        - (1/DT + 1/DI) * pop[index['Ish' + str(j)]]
    return x 

def dIshNI(pop, DI, DT, NIh, index):
    x =   (1/DI) * pop[index['Ish' + str(NIh-1)]] \
        - (1/DT) * pop[index['Ish' + str(NIh)]]
    return x 

def dI_i1(pop, fi, DP, DT, DI, NP, index):
    x =   fi * (1/DP) * pop[index['P__' + str(NP)]] \
        + fi * (1/DP) * pop[index['Ps_' + str(NP)]] \
        +      (1/DP) * pop[index['Pt_' + str(NP)]] \
        +      (1/DT) * pop[index['Isp1']] \
        +      (1/DT) * pop[index['Ish1']] \
        -      (1/DI) * pop[index['I_i1']]      
    return x

def dI_ik(pop, j, DI, DT, index):
    x =   (1/DI) * pop[index['I_i' + str(j-1)]] \
        + (1/DT) * pop[index['Isp' + str(j)]] \
        + (1/DT) * pop[index['Ish' + str(j)]] \
        - (1/DI) * pop[index['I_i' + str(j)]]
    return x
        

def dR(pop,fdead_p, fdead_h, fdead_i, DI, NIp, NIh, NIi, index, vac):
    x = (1/DI)*((1-fdead_p) * pop[index['I_p' + str(NIp)]] \
              + (1-fdead_h) * pop[index['I_h' + str(NIh)]] \
              + (1-fdead_i) * pop[index['I_i' + str(NIi)]])\
            + vac
    return x

def dF(pop, fdead_p, fdead_h, DI, DF, NIp, NIh, index,f_di_h, f_di_p):
    x = (1/DI) * ((1-f_di_p) * fdead_p * pop[index['I_p' + str(NIp)]] \
               +  (1-f_di_h) * fdead_h * pop[index['I_h' + str(NIh)]]) \
                - (1/DF)  * pop[index['F__']]
    return x

def dB_f(pop, DF, index):
    x = (1/DF) *   pop[index['F__']] 
    return x

def dB_j(pop, fdead_i, DI, NIi, index, f_di_h, f_di_p):
    x = (1/DI) * (fdead_i * pop[index['I_i' + str(NIi)]] \
                + f_di_p * fdead_p * pop[index['I_p' + str(NIp)]] \
                + f_di_h * fdead_h * pop[index['I_h' + str(NIh)]])
                  # * (1-q(t)) * d_h

    return x
