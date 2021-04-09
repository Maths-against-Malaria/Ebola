# duration of simulation, in days
days = 1000
# number of Erlang stages for the compartments (compartment names below)
NE = NP = NIp = NIh = NIi = 16
# initial population 
N = 10_000
# duration of diseas states in days
## general meaning of E (latent), P (podromal), I (fully infectous), F (dead, not yet buried)
DE = 7
DP = 5.5
DI = 5
DF = 1.5 # time death to be buried
DT = 9   # trace back time

# fraction of infected who dies under the conditions: p (at home), h (in hospital), i (in isolation unit)
fdead_p = 0.7 # f dead home
fdead_h = 0.5 # f dead hosp
fdead_i = 0.3 # f dead iso
# basic reproduction number of the virus
R0 = 2      # basic reproduction number
# contagiousness in different states (P, I, F) and conditions (p, h - for I), cIi = 0
cP = 0.6
cIp = 0.9   # c_I home
cIh = 0.8   # c_I hosp
cF = 1      # relative contagiousness after death before being buried       
# initial number of sick (in state P) individuals
P0 = 1     # P(0)

# countermeasures
# day when isolation starts
t_iso = 250
# fraction at home (p) or in hospital (h) before this day
f_p1 = 0.5
f_h1 = 0.5
# fraction at home, in hospital, in isolation after this day (4 scenarios provided)
f_p2 = [0.5, 0.5, 0.3, 0.1]
f_h2 = [0.3, 0.1, 0.1, 0.1]
# index of the scenario from the ones above
k = 0 
# probability that a person who had contact with somebody who gets into isolation (f_iso = 1- (f_p2 + f_h2)) is traced back
f_tb = 1
# goodness of isolation beyond isolation wards (0 - not better than else at home/hospital)
ph = 0 
# capacity of isolation wards
qmax = N 


