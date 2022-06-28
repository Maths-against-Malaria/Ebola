import numpy as np
from scipy.integrate import solve_ivp

########################################################
# from Aliou Bouba VI 2022


#20)tbtime (fiso=0.6+ftb=0.6+l4+tiso90+m=0)
print('20_21)tbtime (fiso=0.6+ftb=0.6+l4+tiso90+m=0)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90, DT=25)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90, DT=20)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90, DT=15)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90, DT=10)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90, DT=5)

#18)tiso (fiso=0.6+ftb=0.6+l4)
print('18_19)tiso (fiso=0.6+ftb=0.6+l4)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=75)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=60)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=45)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=30)

#16)tiso (fiso=0.6+ftb=0.6)
print('16_-)tiso (fiso=0.6+ftb=0.6)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=75)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=60)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=45)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=30)

#14)ftb (fiso=0.8+l4)
print('14_17)ftb (fiso=0.8+l4)')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#12)ftb (fiso=0.8)
print('12_15)ftb (fiso=0.8)')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#10)ftb (fiso=0.6)
print('10_13)ftb (fiso=0.6)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[1], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[2], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#8)fiso (ftb=0.8+l4)
print('8_11)fiso (ftb=0.8+l4)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#6)fiso (ftb=0.8)
print('6_9)fiso (ftb=0.8)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#4)fiso (ftb=0.6)
print('4_/)fiso (ftb=0.6)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#2')fiso (only+safe funeral)
print('2*_5)fiso (only+safe funeral)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#2)fiso (only)
print('2_3)fiso (only)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[0],fdead_h[0], fdead_i[0]], t_iso=90)

#MILD MORTALITY


#19)tbtime (fiso=0.6+ftb=0.6+l4+tiso50+m=1)
print('19_20)tbtime (fiso=0.6+ftb=0.6+l4+tiso50+m=1)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90, DT=25)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90, DT=20)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90, DT=15)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90, DT=10)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90, DT=5)

#17)tiso (fiso=0.6+ftb=0.6+l4)
print('17_18)tiso (fiso=0.6+ftb=0.6+l4)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=75)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=60)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=45)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=30)

#15)tiso (fiso=0.6+ftb=0.6)
print('15_)tiso (fiso=0.6+ftb=0.6)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=75)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=60)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=45)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=30)

#13)ftb (fiso=0.8+l4)
print('13_16)ftb (fiso=0.8+l4)')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)

#11)ftb (fiso=0.8)
print('11_14)ftb (fiso=0.8)')
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[1], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[2], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)

#9)ftb (fiso=0.6)
print('9_12)ftb (fiso=0.6)')
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[1], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[2], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)

#7)fiso (ftb=0.8+l4)
print('7_10)fiso (ftb=0.8+l4)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)

#5)fiso (ftb=0.8)
print('5_8)fiso (ftb=0.8)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[4], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)

#3)fiso (ftb=0.6)
print('3_6)fiso (ftb=0.6)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[3], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)

#1')fiso (only+safe funeral)
print('1*_4)fiso (only+safe funeral)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[4], d_h1[4]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)

#1)fiso (only)
print('1_2)fiso (only)')
modelEbola(f_ph1 = [f_p1[0], f_h1[0]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[1], f_h1[1]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[2], f_h1[2]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[3], f_h1[3]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
modelEbola(f_ph1 = [f_p1[4], f_h1[4]], f_tb =f_tb_[0], d_ph1 = [d_p1[0], d_h1[0]], fdead = [fdead_p[1],fdead_h[1], fdead_i[1]], t_iso=90)
#'''