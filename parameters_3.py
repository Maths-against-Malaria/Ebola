# duration of simulation, in days
days = 2 * 365

# number of Erlang stages for the compartments (compartment names below)
NE = NP = NIp = NIh = NIi = 16
Nerls_ = [NE, NP, NIp, NIh, NIi]
n = 16

# initial population
N = 10_000 # 1_000_000

# duration of diseas states in days
## general meaning of E (latent), P (podromal), I (fully infectous), F (dead, not yet buried)
DE = 10
DP = 5
DIh = 5
DI= 5 # 3
DF = 2
DT = 15   # trace back time 9

D = [DE, DP, DIh, DI, DF, DT]

# fraction of infected who dies under the conditions: p (at home), h (in hospital), i (in isolation unit)
fdead_p = [0.6, 0.95] # 0.7 # f dead home # 0.5 - 0.95
fdead_h = [0.35, 0.6] # 0.5 # f dead hosp # 0.35 - 0.6
fdead_i = [0.1, 0.3]  # 0.3 # f dead iso  # 0.1- 0.3
m = 0
fdead = [fdead_p[m],fdead_h[m], fdead_i[m]]

# basic reproduction number of the virus
R0 = 2      # basic reproduction number

# contagiousness in different states (P, I, F) and conditions (p, h - for I), cIi = 0
cP = 0.3
cIh = 0.5   # c_I hosp
cI = 0.6    # c_I home
cF = 1      # relative contagiousness after death before being buried

cc = [cP,cIh,cI,cF]

# initial number of sick (in state P) individuals
P0 = 10

# countermeasures
# day when isolation starts
t_iso_ = [7, 14, 21, 30, 45, 60]
t_iso = t_iso_[5]
# number of infected at hospital where isolation starts (first criterion relevant)
I_iso_=[1,2,5,10,20,50,100,1000,N]
I_iso = N

# general reduction of contact with sick
fc = 1 # 0.2

# fraction at home (p) or in hospital (h) before this day
#f_p1 = 0.5
#f_h1 = 0.5
f_ph0 = [0.5, 0.5]
# fraction at home, in hospital, in isolation after this day (4 scenarios provided)
f_p1 = [0.5, 0.5, 0.5, 0.3, 0.1]
f_h1 = [0.5, 0.3, 0.1, 0.1, 0.1]
# index of the scenario from the ones above
k = 3
f_ph1 = [f_p1[k], f_h1[k]]
# probability that a person who had contact with somebody who gets into isolation (f_iso = 1- (f_p2 + f_h2)) is traced back
f_tb_ = [0,0.2,0.4,0.6,0.8,1]
f_tb = f_tb_[2]

# capacity for trace back
#cmax = N    #/100_000 * 200
cmax_ = [0, N/100_000*10, N/100_000*20,N/100_000*50,N/100_000*100,N/100_000*200, N]
cmax = cmax_[6]
# goodness of isolation beyond isolation wards (0 - not better than else at home/hospital)
ph = 0.6

# capacity of isolation wards
#qmax = N    #'/100_000 * 500
#qmax_ = [0, N/100_000*250, N/100_000*500, N/100_000*750, N/100_000*1000]
qmax_ = [0, N/100_000*10, N/100_000*20,N/100_000*50,N/100_000*100,N/100_000*200,N/100_000*500, N]
qmax = qmax_[7]

# probability of safe funeral for those dying
#d_p1 = 0.0  # at home
#d_h1 = 0.0  # in hospital (no isolation)
d_ph0 = [0,0]
d_p1 = [0, 0.04, 0.08, 0.12, 0.16]# at home
d_h1 = [0, 0.2, 0.4, 0.6, 0.8] # in hospital (no isolation)
l = 0
d_ph1 = [d_p1[l], d_h1[l]]

# vaccination
# start date
t_vac=days  # 500

# daily number
N_vac=0


