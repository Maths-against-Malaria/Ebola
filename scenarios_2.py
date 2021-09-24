# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 18:20:23 2020

@author: helle
"""

import numpy as np
from scipy.integrate import solve_ivp

#from parameters_2 import *  # basic scenario (all parameters)
#from index import *
#from functions import *
#from differential_equations import *
from solve_function import *

# in parameters set  NE = 16, NP = 16, NIp = 16, NIh = 16, NIi = 16
# in parameters set  NE = 1, NP = 1, NIp = 1, NIh = 1, NIi = 1

'''
# no countermeasures
modelEbola(k = 0, f_tb = 0, l = 0)

# fiso...f_tb=0; d_ = 0
modelEbola(k = 1, f_tb = 0, l = 0)
modelEbola(k = 2, f_tb = 0, l = 0)
modelEbola(k = 3, f_tb = 0, l = 0)
modelEbola(k = 4, f_tb = 0, l = 0)

# fiso...; f_tb=0.8; d_ = 0
modelEbola(k = 1, f_tb = 0.8, l = 0)
modelEbola(k = 2, f_tb = 0.8, l = 0)
modelEbola(k = 3, f_tb = 0.8, l = 0)
modelEbola(k = 4, f_tb = 0.8, l = 0)

# qmax
modelEbola(k = 3, f_tb = 0.8, l = 0, qmax = 0)
modelEbola(k = 3, f_tb = 0.8, l = 0, qmax = 20)
modelEbola(k = 3, f_tb = 0.8, l = 0, qmax = 40)
modelEbola(k = 3, f_tb = 0.8, l = 0, qmax = 60)
modelEbola(k = 3, f_tb = 0.8, l = 0, qmax = 80)
modelEbola(k = 3, f_tb = 0.8, l = 0, qmax = 100)

# ftb
modelEbola(k = 3, f_tb = 0.0, l = 0)
modelEbola(k = 3, f_tb = 0.2, l = 0)
modelEbola(k = 3, f_tb = 0.4, l = 0)
modelEbola(k = 3, f_tb = 0.6, l = 0)
modelEbola(k = 3, f_tb = 0.6, l = 0)


modelEbola(k = 3, f_tb = 0.6, l = 0, cmax = 0)
modelEbola(k = 3, f_tb = 0.6, l = 0, cmax = 5)
modelEbola(k = 3, f_tb = 0.6, l = 0, cmax = 10)
modelEbola(k = 3, f_tb = 0.6, l = 0, cmax = 15)
modelEbola(k = 3, f_tb = 0.6, l = 0, cmax = 20)


modelEbola(k = 3, f_tb = 0.6, l = 0, cmax = 10)
modelEbola(k = 3, f_tb = 0.6, l = 1, cmax = 10)
modelEbola(k = 3, f_tb = 0.6, l = 2, cmax = 10)
modelEbola(k = 3, f_tb = 0.6, l = 3, cmax = 10)
modelEbola(k = 3, f_tb = 0.6, l = 4, cmax = 10)
modelEbola(k = 3, f_tb = 0.6, l = 5, cmax = 10)

modelEbola(k = 1, f_tb = 0, l = 0)
modelEbola(k = 1, f_tb = 0, l = 1)
modelEbola(k = 1, f_tb = 0, l = 2)
modelEbola(k = 1, f_tb = 0, l = 3)
modelEbola(k = 1, f_tb = 0, l = 4)
modelEbola(k = 1, f_tb = 0, l = 5)

modelEbola(k = 3, f_tb = 0.6, l = 3, cmax = 10, t_iso = 150)
modelEbola(k = 3, f_tb = 0.6, l = 3, cmax = 10, t_iso = 120)
modelEbola(k = 3, f_tb = 0.6, l = 3, cmax = 10, t_iso = 90)
modelEbola(k = 3, f_tb = 0.6, l = 3, cmax = 10, t_iso = 60)
modelEbola(k = 3, f_tb = 0.6, l = 3, cmax = 10, t_iso = 30)
modelEbola(k = 3, f_tb = 0.6, l = 3, cmax = 10, t_iso = 0)


# parameters3
modelEbola(m=0)
modelEbola(m=1)

modelEbola(m=0, t_iso = 7)
modelEbola(m=0, t_iso = 14)
modelEbola(m=0, t_iso = 21)
modelEbola(m=0, t_iso = 28)

modelEbola(m=1, t_iso = 7)
modelEbola(m=1, t_iso = 14)
modelEbola(m=1, t_iso = 21)
modelEbola(m=1, t_iso = 28)
'''
modelEbola(m=1, t_iso = 14)
modelEbola(m=1, t_iso = 28)
modelEbola(m=1, t_iso = 42)
modelEbola(m=1, t_iso = 56)
modelEbola(m=1, t_iso = 70)
