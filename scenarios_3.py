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
# only isolation
modelEbola(k=0)
modelEbola(k=1)
modelEbola(k=2)
modelEbola(k=3)
modelEbola(k=4)

# combined with traceback
modelEbola(k=0, f_tb =0.6)
modelEbola(k=1, f_tb =0.6)
modelEbola(k=2, f_tb =0.6)
modelEbola(k=3, f_tb =0.6)
modelEbola(k=4, f_tb =0.6)

# combined with safe funeral
modelEbola(k=0, f_tb =0.6, l=1)
modelEbola(k=1, f_tb =0.6, l=1)
modelEbola(k=2, f_tb =0.6, l=1)
modelEbola(k=3, f_tb =0.6, l=1)
modelEbola(k=4, f_tb =0.6, l=1)

# combined with earlier intervention
modelEbola(k=0, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=1, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=2, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=3, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=4, f_tb =0.6, l=1, t_iso = 15)
'''
#---------------------------------
'''
# only isolation
modelEbola(k=0)
modelEbola(k=1)
modelEbola(k=2)
modelEbola(k=3)
modelEbola(k=4)

# combined with traceback
modelEbola(k=0)
modelEbola(k=1, f_tb =0.6)
modelEbola(k=2, f_tb =0.6)
modelEbola(k=3, f_tb =0.6)
modelEbola(k=4, f_tb =0.6)

# combined with safe funeral
modelEbola(k=0, f_tb =0.6, l=1)
modelEbola(k=1, f_tb =0.6, l=1)
modelEbola(k=2, f_tb =0.6, l=1)
modelEbola(k=3, f_tb =0.6, l=1)
modelEbola(k=4, f_tb =0.6, l=1)

# combined with earlier intervention
modelEbola(k=0, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=1, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=2, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=3, f_tb =0.6, l=1, t_iso = 15)
modelEbola(k=4, f_tb =0.6, l=1, t_iso = 15)
'''

modelEbola(k=3, f_tb =0.6, l=1, DIp = 5)
modelEbola(k=3, f_tb =0.6, l=1, DIp = 3)
