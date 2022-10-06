# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:13:23 2020

@author: Kristina B. Helle, Aliou Bouba, Kristan A. Schneider
"""

from functions import *
from parameters import *


# labels for the scenarios to be used for plotting
lab_DT = ['$D_T=25$', '$D_T=20$', '$D_T=15$', '$D_T=10$', '$D_T=5$']
lab_tiso = ['$t_{Iso}=90$', '$t_{Iso}=75$', '$t_{Iso}=60$', '$t_{Iso}=45$', '$t_{Iso}=30$']
lab_tb = ['$f_{Tr}=0.0$', '$f_{Tr}=0.2$', '$f_{Tr}=0.4$', '$f_{Tr}=0.6$', '$f_{Tr}=0.8$']
lab_fiso = ['none', '$f_{Iso}=0.0$', '$f_{Iso}=0.2$', '$f_{Iso}=0.4$', '$f_{Iso}=0.6$', '$f_{Iso}=0.8$']
lab_fiso0 = ['none', '$f_{Iso}=0.2$', '$f_{Iso}=0.4$', '$f_{Iso}=0.6$', '$f_{Iso}=0.8$']

# ---------------------- Fig 2 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, no safe funeral: l=0, severe mortality: m=1)
print('Fig2: fisoftb0l0m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]])

names2_fisoftb0l0m1 =[
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = names2_fisoftb0l0m1,
              savename = 'Fig2',
              lab = lab_fiso0, # labels for varied fiso
              tb = False, sf = True)

# ---------------------- Fig 3 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, safe funeral scenario 4: l=4, severe mortality: m=1)

print('Fig3: fisoftb0l0m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[0], d_h1[0]])
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], d_ph1 = [d_p1[4], d_h1[4]])

names3_fisoftb0l4m1 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_6(names = names3_fisoftb0l4m1,
                   savename = 'Fig3',
                   lab = lab_fiso,
                    tb = False, sf = True)

# ---------------------- Fig 4 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, no safe funeral: l=0, severe mortality: m=1)

print('Fig4: fisoftb08l0m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4])

names4_fisoftb08l0m1 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
]


plotEbolaAll_5(names = names4_fisoftb08l0m1,
              savename = 'Fig4',
              lab = lab_fiso0,
               tb = True, sf = True)


# ---------------------- Fig 5 -----------------------------------------------------------------------------------------
# trace back fraction (ftb) varied (fraction of treatment scenario 4 (fiso=0.8), safe funeral scenario 4: l=4, severe mortality: m=1)

print('Fig5: ftbfiso08l4m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]])

names5_ftbfiso08l4m1 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.2_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.4_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.6_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]


plotEbolaAll_5(names = names5_ftbfiso08l4m1,
              savename = 'Fig5',
              lab = lab_tb,
               tb = True, sf = True,
               col=colsA1) # color scheme without baseline scenario

# ---------------------- Fig 6 -----------------------------------------------------------------------------------------
# trace back time (DT) varied (fraction of treatment scenario 4 (fiso=0.8), safe funeral scenario 4: l=4, severe mortality: m=1)

print('Fig6: DTfiso08l4m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=25)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=20)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=15)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=10)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=5)

names6_DTfiso08l4m1=[
'730_16_10000_[10, 5, 5, 5, 2, 25]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 20]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 15]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 10]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 5]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = names6_DTfiso08l4m1,
              savename = 'Fig6',
              lab = lab_DT,
               tb = True, sf = True,
               col=colsA1)

# ---------------------- Fig 7 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, safe funeral: l=4, severe mortality: m=1)

print('Fig7: fisoftb08l4m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]])
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]])

names7_fisoftb08l4m1 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_6(names = names7_fisoftb08l4m1,
              savename = 'Fig7',
              lab = lab_fiso,
               tb = True, sf = True)


# ---------------------- Fig 8 -----------------------------------------------------------------------------------------
# onset time of intervention (tiso) varied (isolation scenario 4 (fiso=0.8), trace back fraction: ftb=0.8, safe funeral: l=4, severe mortality: m=1)

print('Fig8: tisofiso08tb08l4m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=75)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=60)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=45)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=30)

names8_tisofiso08tb08l4m1 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_75_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_60_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_45_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_30_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = names8_tisofiso08tb08l4m1,
              savename = 'Fig8',
              lab = lab_tiso,
               tb = True, sf = True,
               col=colsA1)


# ---------------------- Fig S1 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, no safe funeral: l=0, moderate mortality: m=0)
print('FigS1: fisoftb0l0m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])

namesS1_fisoftb0l0m0 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = namesS1_fisoftb0l0m0,
              savename = 'FigS1',
              lab = lab_fiso0,
              tb = False, sf = True)

# ---------------------- Fig S2 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, safe funeral scenario 4: l=4, moderate mortality: m=0)
print('FigS2: fisoftb0l4m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])

namesS2_fisoftb0l4m0 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
]
plotEbolaAll_6(names = namesS2_fisoftb0l4m0,
                   savename = 'FigS2',
                   lab = lab_fiso,
               tb = False, sf = True)

# ---------------------- Fig S3 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, no safe funeral: l=0, moderate mortality: m=0)

print('FigS3: fisoftb08l0m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])

namesS3_fisoftb08l0m0 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
]


plotEbolaAll_5(names = namesS3_fisoftb08l0m0,
              savename = 'FigS3',
              lab = lab_fiso0,
               tb = True, sf = True)

# ---------------------- Fig S4 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, safe funeral: l=4, moderate mortality: m=0)

print('FigS4: fisoftb08l4m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])

namesS4_fisoftb08l4m0 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.5, 0.5]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.45, 0.35]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.35, 0.25]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.25, 0.15]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_6(names = namesS4_fisoftb08l4m0,
              savename = 'FigS4',
              lab = lab_fiso,
               tb = True, sf = True)

# ---------------------- Fig S5 -----------------------------------------------------------------------------------------
# trace back fraction (ftb) varied (fraction of treatment scenario 4 (fiso=0.8), safe funeral: l=4, moderate mortality: m=0)

print('FigS5: ftbfiso08l4m0')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]])

namesS5_ftbfiso08l4m0 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.2_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.4_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.6_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = namesS5_ftbfiso08l4m0,
              savename = 'FigS5',
              lab = lab_tb,
               tb = True, sf = True,
               col=colsA1)

# ---------------------- Fig S6 -----------------------------------------------------------------------------------------
# trace back time (DT) varied (fraction of treatment scenario 4 (fiso=0.8), safe funeral: l=4, moderate mortality: m=0)

print('FigS6: DTfiso08l4m0')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], DT=25)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], DT=20)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], DT=15)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], DT=10)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], DT=5)

namesS6_DTfiso08l4m0 =[
'730_16_10000_[10, 5, 5, 5, 2, 25]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 20]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 15]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 10]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 5]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = namesS6_DTfiso08l4m0,
              savename = 'FigS6',
              lab = lab_DT,
               tb = True, sf = True,
               col=colsA1)

# ---------------------- Fig S7 -----------------------------------------------------------------------------------------
# onset time of intervention (tiso) varied (isolation scenario 4 (fiso=0.8), trace back fraction: ftb=0.8, safe funeral: l=4, moderate mortality: m=0)

print('FigS7: tisofiso08tb08l4m0')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=75)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=60)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=45)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=30)

names_tisofiso08tb08l4m0 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_75_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_60_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_45_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.6, 0.4, 0.1]_1.8_[0.3, 0.6, 0.5, 1]_10_30_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0.16, 0.8]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = names_tisofiso08tb08l4m0,
              savename = 'FigS7',
              lab = lab_tiso,
               tb = True, sf = True,
               col=colsA1)

# ---------------------- Fig S8 -----------------------------------------------------------------------------------------
# trace back varied (fraction of treatment scenario 4 (fiso=0.8), no safe funeral: l=0, severe mortality: m=1)

print('FigS8: ftbfiso08l0m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[0], d_h1[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[0], d_h1[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]])
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]])

namesS8_ftbfiso08l0m1 = [
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0.2_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0.4_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0.6_0.6_10000_10000_1_730_0_',
'730_16_10000_[10, 5, 5, 5, 2, 21]_[0.9, 0.6, 0.3]_1.8_[0.3, 0.6, 0.5, 1]_10_90_10000_[0.5, 0.5]_[0.1, 0.1]_[0, 0]_[0, 0]_0.8_0.6_10000_10000_1_730_0_',
]

plotEbolaAll_5(names = namesS8_ftbfiso08l0m1,
              savename = 'FigS8',
              lab = lab_tb, tb = True, sf = True,
              col=colsA1)



