
from functions2 import *
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

# --------- all ----------------------
plotEbolaAll(names=names0 + [names_fiso[2]] + [names_ftb[3]] + [names_cmax[2]] + [names_d[3]] + [names_tiso[3]], savename='all',Nerls = [16,16,16,16,16],
             lab=['none','f_iso = 0.4','f_tb = 0.6','cmax = 10','d = 0.4', 't_iso = 60'],
             col=colsA, q_max=False, tb=True)

plotEbolaAll(names=names3_base, savename='base3',Nerls = [16,16,16,16,16],
             lab=['m=0','m=1'],
             col=colsA, q_max=False, tb=True)

plotEbolaAll(names=names3_tiso_m1, savename='base3',Nerls = [16,16,16,16,16],
             lab=['tiso = 7','tiso = 14','tiso = 21','tiso = 28',],
             col=colsA, q_max=False, tb=True)

plotEbolaAll(names=names3_Iiso, savename='base3_Iiso',Nerls = [16,16,16,16,16],
             lab=['I_iso = 10','I_iso = 20','I_iso = 30','I_iso = 40','I_iso = 50',],
             col=colsA, q_max=False, tb=True)
'''
# --------ASTMH poster Nov 2021
'''
plotEbolaScenarios(
    namesAllScenarios= [names_ASTMH21_fiso, names_ASTMH21_fiso_tb, names_ASTMH21_fiso_tb_sf, names_ASTMH21_fiso_tb_sf_tiso],
                        savename = 'ASTMH21_1',
                        lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
                        Nerls = [16,16,16,16,16,],
                        pathIn=pathIn,
                        pathOut=pathOut,
                        col=colsA,
)

plotEbolaScenarios(
    namesAllScenarios= [names_ASTMH21_2_fiso, names_ASTMH21_2_fiso_tb, names_ASTMH21_2_fiso_tb_sf, names_ASTMH21_2_fiso_tb_sf_tiso],
                        savename = 'ASTMH21_3',
                        lab=['0.0','0.2','0.4','0.6','0.8'],
                        #lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
                        Nerls = [16,16,16,16,16,],
                        pathIn=pathIn,
                        pathOut=pathOut,
                        col=colsA,
)

names_ASTMH21_3_fiso = names_ASTMH21_2_fiso
names_ASTMH21_3_fiso_tb = names_ASTMH21_2_fiso_tb
names_ASTMH21_3_fiso_tb[0] = names_ASTMH21_3_fiso[0]
names_ASTMH21_3_fiso_tb_sf=names_ASTMH21_2_fiso_tb_sf
names_ASTMH21_3_fiso_tb_sf[0] = names_ASTMH21_3_fiso[0]
names_ASTMH21_3_fiso_tb_sf_tiso=names_ASTMH21_2_fiso_tb_sf_tiso
names_ASTMH21_3_fiso_tb_sf_tiso[0] = names_ASTMH21_3_fiso[0]
plotEbolaScenarios(
    namesAllScenarios= [names_ASTMH21_3_fiso, names_ASTMH21_3_fiso_tb, names_ASTMH21_3_fiso_tb_sf, names_ASTMH21_3_fiso_tb_sf_tiso],
                        savename = 'ASTMH21_3',
                        lab=['0.0','0.2','0.4','0.6','0.8'],
                        #lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
                        Nerls = [16,16,16,16,16,],
                        pathIn=pathIn,
                        pathOut=pathOut,
                        col=colsA,
)

plotEbolaScenarios(
    namesAllScenarios= [names_ASTMH21_4_fiso, names_ASTMH21_4_fiso_tb, names_ASTMH21_4_fiso_tb_sf, names_ASTMH21_4_fiso_tb_sf_tiso],
                        savename = 'ASTMH21_4',
                        lab=['a','b','c','d','e'],
                        #lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
                        Nerls = [16,16,16,16,16,],
                        pathIn=pathIn,
                        pathOut=pathOut,
                        col=colsA,
)

plotEbolaScenarios(
    namesAllScenarios= [names_ASTMH21_4_fiso, names_ASTMH21_4_fiso_tb, names_ASTMH21_4_fiso_tb_sf, names_ASTMH21_4_sf],
                        savename = 'ASTMH21_5',
                        lab=['a','b','c','d','e'],
                        #lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
                        Nerls = [16,16,16,16,16,],
                        pathIn=pathIn,
                        pathOut=pathOut,
                        col=colsA,
)

plotEbolaAll(names=names_DIp, savename='DIp',Nerls = [16,16,16,16,16],
             lab=['DIp =5 sf','DIp = 3 sf','DIp =5','DIp = 3', ],
             col=colsA, q_max=False, tb=True)


plotEbolaScenarios(
    namesAllScenarios= [names_ASTMH21_4_fiso, names_ASTMH21_4_fiso_bugfix],
                        savename = 'ASTMH21_bugfix',
                        lab=['a','b','c','d','e'],
                        #lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
                        Nerls = [16,16,16,16,16,],
                        pathIn=pathIn,
                        pathOut=pathOut,
                        col=colsA,
)
plotEbolaAll(names=names_test_fiso_ftb_tiso, savename='DIp',Nerls = [16,16,16,16,16],
             lab=['a', 'b', 'c', 'd', 'e'],
             col=colsA, q_max=False, tb=True)
             
plotEbolaAll(names=names_qmax2, savename='Qmax',Nerls = [16,16,16,16,16],
             lab=['qmax = 0', 'qmax = 1', 'qmax = 2','qmax = 5','qmax = 10','qmax = 20'],
             col=colsA, q_max=True, tb=True)

plotEbolaAll(names=names_cmax1, savename='Cmax',Nerls = [16,16,16,16,16],
             lab=['no tb', 'cmax = 0', 'cmax = 1', 'cmax = 2','cmax = 5','cmax = 10','cmax = 20', 'unlimited'],
             col=colsA, q_max=True, tb=True)
             
plotEbolaAll(names=names_Aliou_all, savename='Aliou_all',Nerls = [16,16,16,16,16],
             lab=['mix'],
             col=colsA, q_max=True, tb=True)

plotEbolaAll(names=names_qcmax, savename='qcmax',
             lab=['a','b', 'c', 'd', 'e', 'f', 'g'],
              q_max=True, tb=True)
plotEbolaAll(names=names_qcmax, savename='qcmax2',
             lab=['a','b', 'c', 'd', 'e', 'f', 'g'],
             q_max=True, tb=True)

plotEbolaAll(names=names_fiso_Iiso10, savename='fiso_Iiso10',
             lab=['f_iso=0.0','f_iso=0.2','f_iso=0.4','f_iso=0.6','f_iso=0.8'],
              q_max=True, tb=False)

plotEbolaAll(names=names_Iiso, savename='Iiso',
             lab=['Iiso = 1', 'Iiso = 2','Iiso = 5','Iiso = 10','Iiso = 20','Iiso = 50','Iiso = 100',],
              q_max=True, tb=False)

plotEbolaAll(names=names_Iiso , savename='Iiso',n = 16,
             lab=['Iiso = 1', 'Iiso = 2','Iiso = 5','Iiso = 10','Iiso = 20','Iiso = 50'],
              q_max=True, tb=False)

print(np.where(getQ(pathIn=pathIn, name=names_Iiso[1], i=1, days = 730) < 730)[0][0]) # Iiso=1 => t_iso = 4
print(np.where(getQ(pathIn=pathIn, name=names_Iiso[2], i=1, days = 730) < 730)[0][0]) #7
print(np.where(getQ(pathIn=pathIn, name=names_Iiso[3], i=1, days = 730) < 730)[0][0]) #13
print(np.where(getQ(pathIn=pathIn, name=names_Iiso[4], i=1, days = 730) < 730)[0][0]) #25
print(np.where(getQ(pathIn=pathIn, name=names_Iiso[5], i=1, days = 730) < 730)[0][0]) # Iiso = 20 => t_iso = 47

plotEbolaAll(names= names_DT  , savename='DT',n=16,
             lab = ['Tr_time = 30','Tr_time = 25', 'Tr_time = 20', 'Tr_time = 15', 'Tr_time = 10', 'Tr_time = 5'],
             #lab=['a','b', 'c', 'd', 'e', 'f', 'g', 'h'],
             q_max=False, tb=True, sf=True,
             legendout=True, days = 180)
'''
#---------------fiso (+ftb=0.8+l4)------------------
plotEbolaAll(names=names4, savename='fisoftb0.8l4m1tiso90', n = 16,
             lab=['none', 'f_iso=0.2', 'f_iso=0.4', 'f_iso=0.6', 'f_iso=0.8'],
             col=colsA, q_max=False, tb=True, sf=True)
