# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 18:20:23 2020

@author: helle
"""

import numpy as np
from scipy.integrate import solve_ivp

from parameters_1 import *  # basic scenario (all parameters)
from index import *
from functions import *
from differential_equations import *
from solve_function import *

# fiso = 0, 0.2...1; tb=0;
# f_iso__f_tb0_qmaxN
#modelEbola(days=730, k=0, f_tb=0, qmax=N)
#modelEbola(days=730, k=1, f_tb=0, qmax=N)
#modelEbola(days=730, k=2, f_tb=0, qmax=N)
#modelEbola(days=730, k=3, f_tb=0, qmax=N)
#modelEbola(days=730, k=4, f_tb=0, qmax=N)

modelEbola(days=1500, k=0, f_tb=0.8, qmax=N)
#modelEbola(days=730, k=1, f_tb=0.8, qmax=N)
#modelEbola(days=730, k=2, f_tb=0.8, qmax=N)
#modelEbola(days=730, k=3, f_tb=0.8, qmax=N)
#modelEbola(days=730, k=4, f_tb=0.8, qmax=N)

# fiso = 0, 0.2...1; tb=0.2;
# f_iso__f_tb0.2_qmaxN_ph0.6
#modelEbola(days=1500, k=0, f_tb=0.2, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=1, f_tb=0.2, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=2, f_tb=0.2, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=3, f_tb=0.2, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=4, f_tb=0.2, ph = 0.6, qmax=N)

# fiso = 0, 0.2...1; tb=0.8
# f_iso__f_tb0.8_qmaxN_ph0.6
#modelEbola(days=1500, k=0, f_tb=0.8, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=1, f_tb=0.8, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=2, f_tb=0.8, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=3, f_tb=0.8, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=4, f_tb=0.8, ph = 0.6, qmax=N)

# fiso =0.2; tb = 0, 0.2, ...1; qmax=N
# f_tb__f_iso0.2_qmaxN_ph0.6
#modelEbola(days=1500, k=1, f_tb=0.0, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=1, f_tb=0.2, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=1, f_tb=0.4, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=1, f_tb=0.6, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=1, f_tb=0.8, ph = 0.6, qmax=N)

# fiso =0.8; tb = 0, 0.2, ...1; qmax=N
# f_tb__f_iso0.8_qmaxN_ph0.6
#modelEbola(days=1500, k=4, f_tb=0.0, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=4, f_tb=0.2, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=4, f_tb=0.4, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.6, qmax=N)
#modelEbola(days=1500, k=4, f_tb=0.8, ph = 0.6, qmax=N)

# qmax = 1,10,100,1000,10000; ph=0.7
# Qmax__f_iso0.8_f_tb0.6_ph0.6
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.6, qmax=1)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.6, qmax=10)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.6, qmax=100)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.6, qmax=1000)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.6, qmax=10000)

# qmax = 100; ph=0, 0.2, 0.4, 0.6, 0.8, 1
# ph__f_iso0.8_f_tb0.6_Qmax100
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.0, qmax=100)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.2, qmax=100)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.4, qmax=100)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.6, qmax=100)
#modelEbola(days=1500, k=4, f_tb=0.6, ph = 0.8, qmax=100)

#-------
# 'Qmax__f_iso0.8_f_tb1_ph0.6'
#modelEbola(days=1500, k=4, f_tb=1, ph = 0.6, qmax=1)
#modelEbola(days=1500, k=4, f_tb=1, ph = 0.6, qmax=10)
#modelEbola(days=1500, k=4, f_tb=1, ph = 0.6, qmax=100)
#modelEbola(days=1500, k=4, f_tb=1, ph = 0.6, qmax=1000)
#modelEbola(days=1500, k=4, f_tb=1, ph = 0.6, qmax=10000)

# 'Qmax__f_iso0.8_f_tb1_ph0.2'



# 'f_iso__f_QmaxN_f_tb0'
#modelEbola(days=1500, f_tb=0, qmax=N, k=0)
#modelEbola(days=1500, f_tb=0, qmax=N, k=1)
#modelEbola(days=1500, f_tb=0, qmax=N, k=2)
#modelEbola(days=1500, f_tb=0, qmax=N, k=3)
#modelEbola(days=1500, f_tb=0, qmax=N, k=4)

# 'f_dih__f_iso0.8_QmaxN100_f_tb0.6'
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.2)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.4)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.6)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.8)

# 'f_dip__f_iso0.8_QmaxN100_f_tb0.6_f_dih0.6'
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.6, f_di_p =0)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.6, f_di_p =0.2)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.6, f_di_p =0.4)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.6, f_di_p =0.6)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, f_di_h = 0.6, f_di_p =0.8)

# 'Nvac__f_iso0.8_QmaxN100_f_tb0.6'
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/5000)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/2500)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/1000)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/500)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/250)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/100)

# 'Nvac__f_iso0.2_QmaxN100_f_tb0.2'
#modelEbola(days=1500, f_tb=0.2, k=1, qmax=100, Nvac=N/5000)
#modelEbola(days=1500, f_tb=0.2, k=1, qmax=100, Nvac=N/2500)
#modelEbola(days=1500, f_tb=0.2, k=1, qmax=100, Nvac=N/1000)
#modelEbola(days=1500, f_tb=0.2, k=1, qmax=100, Nvac=N/500)
#modelEbola(days=1500, f_tb=0.2, k=1, qmax=100, Nvac=N/250)
#modelEbola(days=1500, f_tb=0.2, k=1, qmax=100, Nvac=N/100)

# 't_vac_f_iso0.8_QmaxN100_f_tb0.6_Nvac10'
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/1000, t_vac=300)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/1000, t_vac=400)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/1000, t_vac=500)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/1000, t_vac=600)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/1000, t_vac=700)
#modelEbola(days=1500, f_tb=0.6, k=4, qmax=100, Nvac=N/1000, t_vac=800)