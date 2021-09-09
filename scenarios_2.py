# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 18:20:23 2020

@author: helle
"""

import numpy as np
from scipy.integrate import solve_ivp

from parameters_2 import *  # basic scenario (all parameters)
from index import *
from functions import *
from differential_equations import *
from solve_function import *

# in parameters set  NE = 16, NP = 16, NIp = 16, NIh = 16, NIi = 16
# in parameters set  NE = 1, NP = 1, NIp = 1, NIh = 1, NIi = 1

'''
# fiso = 0, 0.2...1; f_tb=0; d_ = 0
modelEbola(k=0, f_tb=0, d_p = 0, d_h = 0)
modelEbola(k=1, f_tb=0, d_p = 0, d_h = 0)
modelEbola(k=2, f_tb=0, d_p = 0, d_h = 0)
modelEbola(k=3, f_tb=0, d_p = 0, d_h = 0)
modelEbola(k=4, f_tb=0, d_p = 0, d_h = 0)
'''
'''
#tiso = 0, 60, 120, 180, 240, 300; fiso = 0.8, f_tb = 0, d_ = 0
modelEbola(t_iso = 0, k=4, f_tb=0, d_p = 0, d_h = 0)
modelEbola(t_iso = 60, k=4, f_tb=0, d_p = 0, d_h = 0)
modelEbola(t_iso = 120, k=4, f_tb=0, d_p = 0, d_h = 0)
modelEbola(t_iso = 180, k=4, f_tb=0, d_p = 0, d_h = 0)
modelEbola(t_iso = 240, k=4, f_tb=0, d_p = 0, d_h = 0)
modelEbola(t_iso = 300, k=4, f_tb=0, d_p = 0, d_h = 0)
'''
'''
# f_tb=0, 0.2, ...; fiso = 0.8; d_ = 0
modelEbola(f_tb=0, d_p = 0, d_h = 0)
modelEbola(f_tb=0.2, d_p = 0, d_h = 0)
modelEbola(f_tb=0.4, d_p = 0, d_h = 0)
modelEbola(f_tb=0.6, d_p = 0, d_h = 0)
modelEbola(f_tb=0.8, d_p = 0, d_h = 0)
''''''
# cmax=N/10000, N/1000,...; fiso = 0.8; d_ = 0
modelEbola(f_tb=0.8, d_p = 0, d_h = 0, cmax = N/10000)
modelEbola(f_tb=0.8, d_p = 0, d_h = 0, cmax = N/1000)
modelEbola(f_tb=0.8, d_p = 0, d_h = 0, cmax = N/100)
modelEbola(f_tb=0.8, d_p = 0, d_h = 0, cmax = N/10)
modelEbola(f_tb=0.8, d_p = 0, d_h = 0, cmax = N)
'''

# d_p, d_h;
modelEbola(f_tb=0.8, l=0)
modelEbola(f_tb=0.8, l=1)
modelEbola(f_tb=0.8, l=2)
modelEbola(f_tb=0.8, l=3)
modelEbola(f_tb=0.8, l=4)
