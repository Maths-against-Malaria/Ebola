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
modelEbola(days=1500, k=4, f_tb=1, ph = 0.2, qmax=1)
modelEbola(days=1500, k=4, f_tb=1, ph = 0.2, qmax=10)
modelEbola(days=1500, k=4, f_tb=1, ph = 0.2, qmax=100)
modelEbola(days=1500, k=4, f_tb=1, ph = 0.2, qmax=1000)
modelEbola(days=1500, k=4, f_tb=1, ph = 0.2, qmax=10000)




