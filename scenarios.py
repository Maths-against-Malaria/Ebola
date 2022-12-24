# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:13:23 2020

@author: Kristina B. Helle, Aliou Bouba, Kristan A. Schneider
"""

from functions import *
from parameters import *


#--------- set mortality to mild -------------------------
fdead0 = [fdead_p[0], fdead_h[0], fdead_i[0]] # mild
fdead1 = [fdead_p[1], fdead_h[1], fdead_i[1]] # severe

# ---------------------- Fig 2 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, no safe funeral: l=0, mild mortality: m=0)
print('Fig2: fisoftb0l0m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], fdead = fdead0)

# ---------------------- Fig 3 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, safe funeral scenario 4: l=4, mild mortality: m=0)
print('Fig3: fisoftb0l1m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[0], d_h1[0]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)

# ---------------------- Fig 4 -----------------------------------------------------------------------------------------
# fraction of trace back (ftb) varied (isolation fraction: fiso=0.8, safe funeral: l=1, mild mortality: m=0)
print('Fig4: ftbfiso08l1m0')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)

# ---------------------- Fig 5 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, safe funeral scenario 4: l=4, mild mortality: m=0)
print('Fig5: fisoftb08l1m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead0)

# ---------------------- Fig S1 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, no safe funeral: l=0, severe mortality: m=1)
print('FigS1: fisoftb0l0m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], fdead = fdead1)

# ---------------------- Fig S2 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0, safe funeral: l=4, severe mortality: m=1)
print('FigS2: fisoftb0l1m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[0], d_h1[0]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)

# ---------------------- Fig S3 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, no safe funeral: l=0, mild mortality: m=0)
print('FigS3: fisoftb08l0m0')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], fdead = fdead0)

# ---------------------- Fig S4 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, no safe funeral: l=0, severe mortality: m=1)
print('FigS4: fisoftb08l0m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], fdead = fdead1)

# ---------------------- Fig S5 -----------------------------------------------------------------------------------------
# fraction of trace back (ftb) varied (isolation fraction: fiso=0.8, safe funeral: l=1, severe mortality: m=1)
print('FigS5: ftbfiso08l1m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)

# ---------------------- Fig S6 -----------------------------------------------------------------------------------------
# fraction of trace back (ftb) varied (isolation fraction: fiso=0.8, no safe funeral: l=0, severe mortality: m=1)
print('FigS6: ftbfiso08l0m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)

# ---------------------- Fig S7 -----------------------------------------------------------------------------------------
# duration of trace back (DT) varied (isolation fraction: fiso=0.8, trace back fraction: ftb=0.8; safe funeral: l=1, mild mortality: m=0)
print('FigS7: DTfiso08ftb08l1m0')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=25, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=20, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=15, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=10, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=5, fdead = fdead0)

# ---------------------- Fig S8 -----------------------------------------------------------------------------------------
# duration of trace back (DT) varied (isolation fraction: fiso=0.8, trace back fraction: ftb=0.8; safe funeral: l=1, severe mortality: m=1)
print('FigS8: DTfiso08ftb08l1m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=25, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=20, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=15, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=10, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], DT=5, fdead = fdead1)

# ---------------------- Fig S9 -----------------------------------------------------------------------------------------
# fraction in treatments (fiso) varied (trace back fraction: ftb=0.8, safe funeral: l=1, severe mortality: m=1)
print('FigS9: fisoftb08l1m1')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = fdead1)

# ---------------------- Fig S10 -----------------------------------------------------------------------------------------
# onset of measures (tiso) varied (isolation fraction: fiso=0.8; trace back fraction: ftb=0.8, safe funeral: l=1, mild mortality: m=0)
print('FigS10: tisofiso08ftb08l1m0')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=90, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=75, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=60, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=45, fdead = fdead0)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=30, fdead = fdead0)

# ---------------------- Fig S11 -----------------------------------------------------------------------------------------
# onset of measures (tiso) varied (isolation fraction: fiso=0.8; trace back fraction: ftb=0.8, safe funeral: l=1, severe mortality: m=1)
print('FigS11: tisofiso08ftb08l1m1')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=90, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=75, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=60, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=45, fdead = fdead1)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], t_iso=30, fdead = fdead1)
