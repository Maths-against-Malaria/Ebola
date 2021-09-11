
from plot import *
from plotnames2 import *

'''
# -------------- fiso -----------------------
plotEbolaAll(names=names0 + names_fiso, savename='fiso',Nerls = [16,16,16,16,16],
             lab=['none','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
             col=colsA, q_max=False, tb=False)

plotEbolaAll(names=names0 + names_fiso2, savename='fiso_ftb0.8',Nerls = [16,16,16,16,16],
             lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
             col=colsA, q_max=False, tb=True)

# ----------- qmax ----------------
plotEbolaAll(names=names_qmax, savename='fiso_qmax',Nerls = [16,16,16,16,16],
             lab=['qmax = 0','qmax = 20','qmax = 40','qmax = 60','qmax = 80','qmax = 100'],
             col=colsA, q_max=True, tb=True, days=365)

# -------------- ftb --------------------
plotEbolaAll(names=names_ftb, savename='ftb',Nerls = [16,16,16,16,16],
             lab=['f_tb=0','f_tb=0.2','f_tb=0.4','f_tb=0.6','f_tb=0.8'],
             col=colsA, q_max=False, tb=True)

# ------------- cmax ------------------
plotEbolaAll(names=names_cmax, savename='cmax',Nerls = [16,16,16,16,16],
             lab=['cmax = 0','cmax = 5','cmax = 10','cmax = 15','cmax = 20',],
             col=colsA, q_max=False, tb=True)
             
# ------------- d ---------------------------
plotEbolaAll(names=names_d, savename='d',Nerls = [16,16,16,16,16],
             lab=['d_p = 0','d_p = 0.1','d_p = 0.2','d_p = 0.4', 'd_p = 0.6', 'd_p = 0.8',],
             col=colsA, q_max=False, tb=True)

plotEbolaAll(names=names_d2, savename='d_tb0',Nerls = [16,16,16,16,16],
             lab=['d_p = 0','d_p = 0.1','d_p = 0.2','d_p = 0.4', 'd_p = 0.6', 'd_p = 0.8',],
             col=colsA, q_max=False, tb=True)

# -------------- tiso ------------------------
plotEbolaAll(names=names_tiso, savename='tiso',Nerls = [16,16,16,16,16],
             lab=['t_iso=150','t_iso=120','t_iso=90','t_iso=60','t_iso=30', 't_iso=0'],
             col=colsA, q_max=False, tb=True)
'''
# --------- all ----------------------
plotEbolaAll(names=names0 + [names_fiso[2]] + [names_ftb[3]] + [names_cmax[2]] + [names_d[3]] + [names_tiso[3]], savename='all',Nerls = [16,16,16,16,16],
             lab=['none','f_iso = 0.4','f_tb = 0.6','cmax = 10','d = 0.4', 't_iso = 60'],
             col=colsA, q_max=False, tb=True)
