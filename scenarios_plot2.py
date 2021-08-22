
from plot import *
from plotnames2 import *

# -------------- fiso -----------------------

plotEbolaAll(names=names_fiso_1, savename='fiso_1',Nerls = [1,1,1,1,1],
             lab=['f_iso=0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
             col=colsA, q_max=False, tb=False)

plotEbolaAll(names=names_fiso_16, savename='fiso_16', Nerls = [16,16,16,16,16],
             lab=['f_iso=0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
             col=colsA, q_max=False, tb=False)

# -------------- tiso --------------------
plotEbolaAll(names=names_tiso_1, savename='tiso_1',Nerls = [1,1,1,1,1],
             lab=['t_iso=0','t_iso=60','t_iso=120','t_iso=180','t_iso=240', 't_iso=300'],
             col=colsA, q_max=False, tb=False)

plotEbolaAll(names=names_tiso_16, savename='tiso_16',Nerls = [16,16,16,16,16],
             lab=['t_iso=0','t_iso=60','t_iso=120','t_iso=180','t_iso=240', 't_iso=300'],
             col=colsA, q_max=False, tb=False)

# -------------- f_tb --------------------
plotEbolaAll(names=names_tb_1, savename='tb_1',Nerls = [1,1,1,1,1],
             lab=['f_tb=0','f_tb=0.2','f_tb=0.4','f_tb=0.6','f_tb=0.8'],
             col=colsA, q_max=False, tb=True)

plotEbolaAll(names=names_tb_16, savename='tb_16',Nerls = [16,16,16,16,16],
             lab=['f_tb=0','f_tb=0.2','f_tb=0.4','f_tb=0.6','f_tb=0.8'],
             col=colsA, q_max=False, tb=True)

# ---------- cmax -------------------
plotEbolaAll(names=names_cmax_1, savename='cmax_1',Nerls = [1,1,1,1,1],
             lab=['cmax=N/10000','cmax=N/1000','cmax=N/100','cmax=N/10','cmax=N/1',],
             col=colsA, q_max=False, tb=True)
plotEbolaAll(names=names_cmax_16, savename='cmax_16',Nerls = [16,16,16,16,16],
             lab=['cmax=N/10000','cmax=N/1000','cmax=N/100','cmax=N/10','cmax=N/1',],
             col=colsA, q_max=False, tb=True)

# ------------- d_p, d_h ---------------------
plotEbolaAll(names=names_d_1, savename='d_1',Nerls = [1,1,1,1,1],
             lab=['d_p = 0','d_p = 0.1','d_p = 0.2','d_p = 0.3','d_p = 0.4','d_p = 0.5'],
             col=colsA, q_max=False, tb=True)
plotEbolaAll(names=names_d_16, savename='d_16',Nerls = [16,16,16,16,16],
             lab=['d_p = 0','d_p = 0.1','d_p = 0.2','d_p = 0.3','d_p = 0.4','d_p = 0.5'],
             col=colsA, q_max=False, tb=True)