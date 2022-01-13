import matplotlib.pyplot as plt
import numpy as np
#from scipy.integrate import solve_ivp
#from plot import *
#from scenarios_plot import names
'''
popSum = np.empty(shape=[len(names), s[0], s[1]])
popSum[0] = popSum0
for i in range(1, len(names)):
    # pops_i = np.loadtxt(path +'/ebola_' + names[i] + '.txt')
    pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
    # plt.plot(np.sum(pops_i, axis=0), color=col[i])
    popSum[i] = popsum2d(pops_i, Nerls=Nerls)

days=1500
pathIn = 'results'
name= '1500_16_16_16_16_16_16_10000_7_5.5_5_1.5_9_0.7_0.5_0.3_2_0.6_0.9_0.8_1_1_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1]_4_0.8_0.0_1000_'
#name = '1500_16_16_16_16_16_16_10000_7_5.5_5_1.5_9_0.7_0.5_0.3_2_0.6_0.9_0.8_1_1_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1]_4_0.8_0.0_1000_'
#name = '1500_16_16_16_16_16_16_10000_7_5.5_5_1.5_9_0.7_0.5_0.3_2_0.6_0.9_0.8_1_1_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1]_4_0.8_0.8_1000_'
pathV = pathIn + '/ebolaVar_' + name + '.txt'
q_ = np.loadtxt(pathV)
q__=np.delete(q_, np.where(q_ == [-1.00000000e+04, -1.00000000e+04]), axis=0)
q___=np.interp(x=np.arange(days), xp=q__[:,0], fp=q__[:,1])
plt.plot(q___)
plt.show()

dummy= 1
NEPIL=2
days=3
b=[[dummy for i in np.arange(15 + 1 + NEPIL)] for j in np.arange(days + 1)]
a=[[0,0]]
#print(a.class)
#print(b.class)
pop=[10,0]

def f(t,out):
   out[0] = 0.1*out[0] + 0.9*out[1]
   out[1] = 0.9*out[0] + 0.1*out[1]
   return out

def fevents():
   y[1]*10

sol=solve_ivp(f, [0,100], pop, t_eval=np.arange(0,100), dense_output=True, event=fevents())

print(sol)


'''

'''
pathIn= 'results'
nameV = 'ebolaVar_1500_16_16_16_16_16_16_10000_7_5.5_5_1.5_9_0.7_0.5_0.3_2_0.6_0.9_0.8_1_1_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1]_4_0.8_0.6_1000_.txt'
nameI = 'ebolaInt_1500_16_16_16_16_16_16_10000_7_5.5_5_1.5_9_0.7_0.5_0.3_2_0.6_0.9_0.8_1_1_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1]_4_0.8_0.6_1000_.txt'
name =  'ebola_1500_16_16_16_16_16_16_10000_7_5.5_5_1.5_9_0.7_0.5_0.3_2_0.6_0.9_0.8_1_1_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1]_4_0.8_0.6_1000_.txt'

ints = np.loadtxt(pathIn + '/' + nameI)
ints= ints[:-2]
var = np.loadtxt(pathIn + '/' + nameV)
popul = np.loadtxt(pathIn + '/' + name)
pop= popsum2d(popul, Nerls=Nerls)
Ii = [pop[9,int(i)] for i in ints]
Pt = [pop[5,int(i)] for i in ints]
vars = [var[int(i),3] for i in ints]
varInterpol = np.interp(x=np.arange(1500), xp=ints, fp=vars)
varsN = np.multiply(vars, -1) + 1
#plt.plot(ints, np.multiply(Ii, varsN), linestyle= '--', color = colsA[1])
#plt.plot(ints, np.multiply(Ii, vars), color = colsA[1])
#plt.plot(ints, np.multiply(Pt, varsN), linestyle= '--', color = colsA[3])
#plt.plot(ints, np.multiply(Pt, vars), color = colsA[3])

plt.plot(np.multiply(pop[9,:],varInterpol), linestyle= '--', color = colsA[1])
plt.plot(np.multiply(pop[5,:], varInterpol), linestyle= '--', color = colsA[3])
#plt.plot(ints, np.multiply(Pt, np.sum([-1], vars)))
plt.show()

'''
'''
val = np.loadtxt(pathIn + '/' + name)
vals = popsum2d(pops=val, Nerls=Nerls)
vars = [vals[int(i),4] for i in ints]
plt.plot(vals)
plt.show()

ebola_1000_1_1_1_1_1_10000_9_5_5_2_9_0.7_0.5_0.3_2_0.3_0.6_0.5_1_25_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1, 0]_4_0_0.8_10000_100.0_0_0_1_1000_0_
ebola_1000_1_1_1_1_1_10000_9_5_5_2_9_0.7_0.5_0.3_1.7_0.5_0.7_0.6_1_25_250_0.5_0.5_[0.5, 0.5, 0.5, 0.3, 0.1, 0.1]_[0.5, 0.3, 0.1, 0.1, 0.1, 0]_0_0_0.8_1000.0_10.0_0_0_1_500_0_.txt not found.
'''
for i in np.arange(10):
    print(np.floor(i/3))

