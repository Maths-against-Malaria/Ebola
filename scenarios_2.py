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
#from solve_function import *

from functions2 import *
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

modelEbola(m=1, I_iso = 14)
modelEbola(m=1, I_iso = 28)
modelEbola(m=1, I_iso = 42)
modelEbola(m=1, I_iso = 56)
modelEbola(m=1, I_iso = 70)

'''
'''
# ASTMH poster Nov 2021
#days= 365, t_iso = 30
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

'''

# only isolation
modelEbola(k=0)
modelEbola(k=1)
modelEbola(k=2)
modelEbola(k=3)
modelEbola(k=4)


# combined with traceback
modelEbola(k=3, f_tb =0.0)
modelEbola(k=3, f_tb =0.2)
modelEbola(k=3, f_tb =0.4)
modelEbola(k=3, f_tb =0.6)
modelEbola(k=3, f_tb =0.8)

# combined with safe funeral
modelEbola(k=3, f_tb =0.6, l=0)
modelEbola(k=3, f_tb =0.6, l=1)
modelEbola(k=3, f_tb =0.6, l=2)
modelEbola(k=3, f_tb =0.6, l=3)
modelEbola(k=3, f_tb =0.6, l=4)

# combined with earlier intervention
modelEbola(k=3, f_tb =0.6, l=3, t_iso = 15)
modelEbola(k=3, f_tb =0.6, l=3, t_iso = 30)
modelEbola(k=3, f_tb =0.6, l=3, t_iso = 45)
modelEbola(k=3, f_tb =0.6, l=3, t_iso = 60)
modelEbola(k=3, f_tb =0.6, l=3, t_iso = 75)

# only safe funeral
modelEbola(k=0, f_tb =0, l=0)
modelEbola(k=0, f_tb =0, l=1)
modelEbola(k=0, f_tb =0, l=2)
modelEbola(k=0, f_tb =0, l=3)
modelEbola(k=0, f_tb =0, l=4)

# DIp
modelEbola(k=3, f_tb =0.6, l=1, DIp = 5)
modelEbola(k=3, f_tb =0.6, l=1, DIp = 3)

modelEbola(k=3, f_tb =0.6, l=0, DIp = 5)
modelEbola(k=3, f_tb =0.6, l=0, DIp = 3)

# the following scenarios should all be the same (they are)
modelEbola(k=0, f_tb =0)
modelEbola(k=0, f_tb =0.8)
modelEbola(k=0, f_tb =0.8,t_iso=0)
modelEbola(k=0, f_tb =0.8,t_iso=730)
modelEbola(k=4, f_tb =0.8,t_iso=730)


modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, qmax = qmax_[0])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, qmax = qmax_[1])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, qmax = qmax_[2])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, qmax = qmax_[3])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, qmax = qmax_[4])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, qmax = qmax_[5])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, qmax = qmax_[6])
modelEbola(k=0, l = 1)

modelEbola(k=3, f_tb =0, t_iso = 60, l = 1)
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, cmax = cmax_[0])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, cmax = cmax_[1])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, cmax = cmax_[2])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, cmax = cmax_[3])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, cmax = cmax_[4])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, cmax = cmax_[5])
modelEbola(k=3, f_tb =0.6, t_iso = 60, l = 1, cmax = cmax_[6])

modelEbola(k=3, f_tb =0.6, t_iso = 45, l = 1, cmax = cmax_[6], qmax = qmax_[6])


k = 3
l = 1
f_ph1 = [f_p1[k], f_h1[k]]
d_ph1 = [d_p1[l], d_h1[l]]
modelEbola(f_ph1=f_ph1, f_tb =0.6, t_iso = 45, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7])
modelEbola(f_ph1=f_ph1, f_tb =0.6, t_iso = 45, d_ph1=d_ph1, cmax = cmax_[5], qmax = qmax_[6])
modelEbola(f_ph1=f_ph1, f_tb =0.6, t_iso = 45, d_ph1=d_ph1, cmax = cmax_[4], qmax = qmax_[5])
modelEbola(f_ph1=f_ph1, f_tb =0.6, t_iso = 45, d_ph1=d_ph1, cmax = cmax_[3], qmax = qmax_[4])
modelEbola(f_ph1=f_ph1, f_tb =0.6, t_iso = 45, d_ph1=d_ph1, cmax = cmax_[2], qmax = qmax_[3])
modelEbola(f_ph1=f_ph1, f_tb =0.6, t_iso = 45, d_ph1=d_ph1, cmax = cmax_[1], qmax = qmax_[2])
modelEbola(f_ph1=f_ph1, f_tb =0.6, t_iso = 45, d_ph1=d_ph1, cmax = cmax_[0], qmax = qmax_[1])


modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 730, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 1)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 730, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 2)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 730, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 5)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 730, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 10)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 730, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 20)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 730, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 50)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 730, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 100)

modelEbola(f_ph1=[f_p1[0], f_h1[0]], f_tb =0, t_iso = 30, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= N)
modelEbola(f_ph1=[f_p1[1], f_h1[1]], f_tb =0, t_iso = 30, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= N)
modelEbola(f_ph1=[f_p1[2], f_h1[2]], f_tb =0, t_iso = 30, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= N)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0, t_iso = 30, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= N)
modelEbola(f_ph1=[f_p1[4], f_h1[4]], f_tb =0, t_iso = 30, d_ph1=d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= N)


#----------------
from plot import getQ, pathIn

x1=modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 730, d_ph1 = d_ph1, cmax = cmax_[0], qmax = qmax_[0], I_iso= I_iso_[0])
x2=modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 730, d_ph1 = d_ph1, cmax = cmax_[1], qmax = qmax_[0], I_iso= I_iso_[1])
x3=modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 730, d_ph1 = d_ph1, cmax = cmax_[2], qmax = qmax_[0], I_iso= I_iso_[2])
x4=modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 730, d_ph1 = d_ph1, cmax = cmax_[3], qmax = qmax_[0], I_iso= I_iso_[3])
x5=modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 730, d_ph1 = d_ph1, cmax = cmax_[4], qmax = qmax_[0], I_iso= I_iso_[4])
x6=modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 730, d_ph1 = d_ph1, cmax = cmax_[5], qmax = qmax_[0], I_iso= I_iso_[5])


modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], t_iso= t1, I_iso = N)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], t_iso= t2, I_iso = N)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], t_iso= t5, I_iso = N)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], t_iso= t10, I_iso = N)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], t_iso= t20, I_iso = N)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], t_iso= t50, I_iso = N)


modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 365, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 2)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 365, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 5)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 365, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 10)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 365, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 20)
modelEbola(f_ph1=[f_p1[3], f_h1[3]], f_tb =0.6, t_iso = 365, d_ph1 = d_ph1, cmax = cmax_[6], qmax = qmax_[7], I_iso= 50)
'''

'''
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "DOP853",f_tb = f_tb_[0])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "DOP853",f_tb = f_tb_[1])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "DOP853",f_tb = f_tb_[2])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "DOP853",f_tb = f_tb_[3])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "DOP853",f_tb = f_tb_[4])

# does not work
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "Radau",f_tb = f_tb_[0])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "Radau",f_tb = f_tb_[1])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "Radau",f_tb = f_tb_[2])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "Radau",f_tb = f_tb_[3])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "Radau",f_tb = f_tb_[4])

modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "BDF",f_tb = f_tb_[0])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "BDF",f_tb = f_tb_[1])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "BDF",f_tb = f_tb_[2])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "BDF",f_tb = f_tb_[3])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "BDF",f_tb = f_tb_[4])

modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "LSODA",f_tb = f_tb_[0])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "LSODA",f_tb = f_tb_[1])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "LSODA",f_tb = f_tb_[2])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "LSODA",f_tb = f_tb_[3])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], method = "LSODA",f_tb = f_tb_[4])


modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], method = "Radau", f_tb = f_tb_[4])
modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], method = "BDF", f_tb = f_tb_[4])
modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], method = "LSODA", f_tb = f_tb_[4])

modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[0])
modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[1])
modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[2])
modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[3])
modelEbola(days = 120, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[4])
'''
modelEbola(days = 365, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[3],d_ph1=[d_p1[4],d_h1[4]],fdead = [fdead_p[1],fdead_h[1],fdead_i[1]],t_iso=50, DT=10)
modelEbola(days = 365, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[3],d_ph1=[d_p1[4],d_h1[4]],fdead = [fdead_p[1],fdead_h[1],fdead_i[1]],t_iso=50, DT=15)
modelEbola(days = 365, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[3],d_ph1=[d_p1[4],d_h1[4]],fdead = [fdead_p[1],fdead_h[1],fdead_i[1]],t_iso=50, DT=20)
modelEbola(days = 365, f_ph1 = [f_p1[3], f_h1[3]], f_tb = f_tb_[3],d_ph1=[d_p1[4],d_h1[4]],fdead = [fdead_p[1],fdead_h[1],fdead_i[1]],t_iso=50, DT=25)

