# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:13:23 2020

@author: Kristina B. Helle, Aliou Bouba, Kristan A. Schneider
"""

# # General conditions

# paths: folders have to exist
# simulation results (to be used as plot input)
pathResult = 'result'

# duration of simulation, in days
days = 2 * 365

# initial population
N = 10_000

# # Diseease characteristics
# states: E (latent), P (podromal), I (fully infectous), F (dead, not yet buried)
# treatments (only applicable to state I): p (at home), h (in hospital), i (in isolation unit)

# number of Erlang stages in each disease state
NE = NP = NIp = NIh = NIi = 16
Nerls_ = [NE, NP, NIp, NIh, NIi]
n = 16

# duration of disease states in days
DE = 10
DP = 5
DIh = 5
DI = 5
DF = 2
DT = 21   # trace back time
D = [DE, DP, DIh, DI, DF, DT]


# basic reproduction number of the virus
R0 = 1.8

# contagiousness in different states (P, I, F) and conditions (p, h - for I), cIi = 0
cP = 0.3
cIh = 0.5
cI = 0.6
cF = 1
cc = [cP, cI, cIh, cF]

# initial number of sick (in state P) individuals
P0 = 10

# # countermeasures

# initialisation of countermeasures
# countermeasures in place: after day
t_iso_ = [30, 45, 60, 75, 90]   # scenarios to choose
t_iso = t_iso_[4]               # scenario chosen

# countermeasures in place: when number of infected at hospital has passed I_iso
I_iso_ = [1, 2, 5, 10, 20, 50, 100, 1000, N]   # scenarios to choose
I_iso = N                                      # scenario chosen

# form of countermeasures
# general reduction of contacts with sick (when countermeasures are in place)
fc = 1

# fraction of treatment when no countermeasures are in place
f_ph0 = [0.5, 0.5]

# fraction of treatment when countermeasures are in place (p: at home, h: in hospital; remaining in isolation)
f_p1 = [0.5, 0.45, 0.35, 0.25, 0.1]    # scenarios to choose
f_h1 = [0.5, 0.35, 0.25, 0.15, 0.1]    # scenarios to choose
k = 0
f_ph1 = [f_p1[k], f_h1[k]]             # scenario chosen


# probability that a person who had contact with somebody who gets into isolation
# (f_iso = 1- (f_p2 + f_h2)) is traced back
f_tb_ = [0, 0.2, 0.4, 0.6, 0.8, 1]
f_tb = f_tb_[0]

# probability of safe funeral for those dying in different treatments (p: at home, h: in hospital)
# when countermeasures are not in place ([0] at home, [1] in hospital)
d_ph0 = [0,0]
# when countermeasures are in place
d_p1 = [0, 0.04, 0.08, 0.12, 0.16]    # scenarios to choose
d_h1 = [0, 0.2, 0.4, 0.6, 0.8]        # scenarios to choose
l = 0
d_ph1 = [d_p1[l], d_h1[l]]            # scenario chosen

# fraction of infected who die under the respective treatment
fdead_p = [0.6, 0.90]
fdead_h = [0.4, 0.6]
fdead_i = [0.1, 0.3]
m = 1 # m=0: moderate mortality; m=1: severe mortality
fdead = [fdead_p[m], fdead_h[m], fdead_i[m]]    # chosen mortality

# --------- The following parameters were not used in the simulations of the article ------------------
# these parameters are discussed in the article
# capacity for trace back
cmax_ = [0, N/100_000*10, N/100_000*20,N/100_000*50,N/100_000*100,N/100_000*200, N]     # scenarios to choose
cmax = cmax_[6]                                                                         # scenario chosen

# capacity of isolation wards
qmax_ = [0, N/100_000*10, N/100_000*20,N/100_000*50,N/100_000*100,N/100_000*200,N/100_000*500, N]     # scenarios to choose
qmax = qmax_[7]                                                                                       # scenario chosen
# goodness of isolation beyond isolation wards (0 - not better than else at home/hospital)
ph = 0.6

# these parameters are not discussed in the article model, the related part of the model is very simplified
# vaccination
# start date
t_vac=days
# daily number
N_vac=0
