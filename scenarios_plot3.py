from functions2 import *
from plotnames3 import *

def plotEbolaAll_(names,
              savename,
              lab,
              tb,
              sf):
    plotEbolaAllCumulative(names = names,
                   savename = savename,
                   lab = lab,
                   n=16,
                   col=colsA,
                   q_max=False,
                   tb=tb,
                   sf=sf,
                   legendlong=False)

lab_DT = ['$D_T=25$', '$D_T=20$', '$D_T=15$', '$D_T=10$', '$D_T=5$']
lab_tiso = ['$t_{Iso}=90$', '$t_{Iso}=75$', '$t_{Iso}=60$', '$t_{Iso}=45$', '$t_{Iso}=30$']
lab_tb = ['none', '$f_{Tr}=0.2$', '$f_{Tr}=0.4$', '$f_{Tr}=0.6$', '$f_{Tr}=0.8$']
lab_fiso = ['none', '$f_{Iso}=0.2$', '$f_{Iso}=0.4$', '$f_{Iso}=0.6$', '$f_{Iso}=0.8$']


#plotEbolaAll_(names = [names1__2_fisoftb0l0m1[0],names1__2_fisoftb0l0m1[4], names2__3_fisoftb0l0m0[0],names2__3_fisoftb0l0m0[4]],
#              savename = 'm1_m0',
#              lab = ['severe_iso0', 'severe_iso08', 'mild_iso0', 'mild_iso08'], tb = False, sf = False)
'''
#--------------- 1) fiso (only) -------------------------------------
plotEbolaAll_(names = names1__2_fisoftb0l0m1,
              savename = 'fisoftb0l0m1_',
              lab = lab_fiso, tb = False, sf = True)

plotEbolaAll_(names = names2__3_fisoftb0l0m0,
              savename = 'fisoftb0l0m0_',
              lab = lab_fiso, tb = False, sf = True)

# ------------- 2*) fiso (only+safe funeral) ---------------------
plotEbolaAll_(names = names1_4_fisoftb0l4m1,
              savename = 'fisoftb0l4m1_',
              lab = lab_fiso, tb = False, sf = True)

plotEbolaAll_(names = names2_5_fisoftb0l4m0,
              savename = 'fisoftb0l4m0_',
              lab = lab_fiso, tb = False, sf = True)

# --------------5) fiso (ftb=0.8) ---------------------

plotEbolaAll_(names = names5_8_fisoftb08l0m1,
              savename = 'fisoftb08l0m1_',
              lab = lab_fiso, tb = True, sf = True)

plotEbolaAll_(names = names6_9_fisoftb08l0m0,
              savename = 'fisoftb08l0m0_',
              lab = lab_fiso, tb = True, sf = True)

#-----------7) fiso (ftb=0.8+l4)

plotEbolaAll_(names = names7_10_fisoftb08l4m1,
              savename = 'fisoftb08l4m1_',
              lab = lab_fiso, tb = True, sf = True)

plotEbolaAll_(names = names8_11_fisoftb08l4m0,
              savename = 'fisoftb08l4m0_',
              lab = lab_fiso, tb = True, sf = True)

#--------------11)ftb (fiso=0.8)

plotEbolaAll_(names = names11_14_ftbfiso08l0m1,
              savename = 'ftbfiso08l0m1_',
              lab = lab_tb, tb = True, sf = True)
plotEbolaAll_(names = names12_15_ftbfiso08l0m0,
              savename = 'ftbfiso08l0m0_',
              lab = lab_tb, tb = True, sf = True)

#-----------------13)ftb (fiso=0.8+l4)
plotEbolaAll_(names = names13_16_ftbfiso08l4m1,
              savename = 'ftbfiso08l4m1_',
              lab = lab_tb, tb = True, sf = True)
plotEbolaAll_(names = names14_17_ftbfiso08l4m0,
              savename = 'ftbfiso08l4m0_',
              lab = lab_tb, tb = True, sf = True)

#----------17)tiso (fiso=0.6+ftb=0.6+l4)
plotEbolaAll_(names = names17_18_tisofiso06ftb06l4m1,
              savename = 'tisofiso06ftb06l4m1_',
              lab = lab_tiso, tb = True, sf = True)
plotEbolaAll_(names = names18_19_tisofiso06ftb06l4m0,
              savename = 'tisofiso06ftb06l4m0_',
              lab = lab_tiso, tb = True, sf = True)

#-----------19)tbtime (fiso=0.6+ftb=0.6+l4+tiso50+m=1)
plotEbolaAll_(names = names19_20_DTfiso06ftb06l4m1,
              savename = 'DTfiso06ftb06l4m1_',
              lab = lab_DT, tb = True, sf = True)
plotEbolaAll_(names = names20_21_DTfiso06ftb06l4m0,
              savename = 'DTfiso06ftb06l4m0_',
              lab = lab_DT, tb = True, sf = True)

'''

