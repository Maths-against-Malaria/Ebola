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



# basic
#modelEbola(t_iso = 1000)
#modelEbola()
#modelEbola(k=1)
#modelEbola(k=2)
#modelEbola(k=3)


# longer
#modelEbola(days=2000, t_iso=1000)
#modelEbola(days=2000)
#modelEbola(days=2000, k=1)
#modelEbola(days=2000, k=2)
#modelEbola(days=2000, k=3)

# with control over traceback
#modelEbola(days=2000,f_tb=0, k=0)
#modelEbola(days=2000,f_tb=0.2, k=0)
#modelEbola(days=2000,f_tb=0.4, k=0)
#modelEbola(days=2000,f_tb=0.6, k=0)
#modelEbola(days=2000,f_tb=0.8, k=0)
#modelEbola(days=2000,f_tb=1, k=0)

#modelEbola(days=2000,f_tb=0, k=3)
#modelEbola(days=2000,f_tb=0.2, k=3)
#modelEbola(days=2000,f_tb=0.4, k=3)
#modelEbola(days=2000,f_tb=0.6, k=3)
#modelEbola(days=2000,f_tb=0.8, k=3)
#modelEbola(days=2000,f_tb=1, k=3)

# limited quarantine wards
#modelEbola(days=2000,f_tb=1, k=3, ph=0.5, qmax=0)
#modelEbola(days=2000,f_tb=1, k=3, ph=0.5, qmax=500)
#modelEbola(days=2000,f_tb=1, k=3, ph=0.5, qmax=1000)
#modelEbola(days=2000,f_tb=1, k=3, ph=0.5, qmax=2000)
#modelEbola(days=2000,f_tb=1, k=3, ph=0.5, qmax=4000)

# 1 Erlang stage: compare to poster
#modelEbola(days=1000,f_tb=1, qmax=10000, NE = 1, NP = 1, NIp = 1, NIh = 1, NIi = 1, t_iso = 1000)
#modelEbola(days=1000,f_tb=1, qmax=10000, NE = 1, NP = 1, NIp = 1, NIh = 1, NIi = 1, k=0)
#modelEbola(days=1000,f_tb=1, qmax=10000, NE = 1, NP = 1, NIp = 1, NIh = 1, NIi = 1, k=1)
#modelEbola(days=1000,f_tb=1, qmax=10000, NE = 1, NP = 1, NIp = 1, NIh = 1, NIi = 1, k=2)
#modelEbola(days=1000,f_tb=1, qmax=10000, NE = 1, NP = 1, NIp = 1, NIh = 1, NIi = 1, k=3)

#modelEbola(days=2000,f_tb=1, qmax=10000, NE = 12, NP = 12, NIp = 12, NIh = 12, NIi = 12, t_iso = 1000)
#modelEbola(days=2000,f_tb=1, qmax=10000, NE = 12, NP = 12, NIp = 12, NIh = 12, NIi = 12, k=0)
#modelEbola(days=2000,f_tb=1, qmax=10000, NE = 12, NP = 12, NIp = 12, NIh = 12, NIi = 12, k=1)
#modelEbola(days=2000,f_tb=1, qmax=10000, NE = 12, NP = 12, NIp = 12, NIh = 12, NIi = 12, k=2)
#modelEbola(days=2000,f_tb=1, qmax=10000, NE = 12, NP = 12, NIp = 12, NIh = 12, NIi = 12, k=3)

#modelEbola(days=2000, t_iso = 2000, f_tb=0)  # fiso = 0 all the time
#modelEbola(days=2000, f_tb=0, k=0) # fiso = 0.2 (after day 250)
#modelEbola(days=2000, f_tb=0, k=1) # fiso = 0.4 (after day 250)
#modelEbola(days=2000, f_tb=0, k=2) # fiso = 0.6 (after day 250)
#modelEbola(days=2000, f_tb=0, k=3) # fiso = 0.8 (after day 250)

modelEbola(days=2000,t_iso=2000)
#modelEbola(days=2000, k=0)
#modelEbola(days=2000, k=1)
#modelEbola(days=2000, k=2)
#modelEbola(days=2000, k=3)