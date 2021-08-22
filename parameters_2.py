# duration of simulation, in days
days = 1000
# number of Erlang stages for the compartments (compartment names below)
NE = NP = NIp = NIh = NIi = 1

# initial population
N = 10_000
# duration of diseas states in days
## general meaning of E (latent), P (podromal), I (fully infectous), F (dead, not yet buried)
DE = 9
DP = 5
DI = 5
DF = 2
DT = 9   # trace back time

# fraction of infected who dies under the conditions: p (at home), h (in hospital), i (in isolation unit)
fdead_p = 0.7 # f dead home
fdead_h = 0.5 # f dead hosp
fdead_i = 0.3 # f dead iso

# basic reproduction number of the virus
R0 = 2      # basic reproduction number

# contagiousness in different states (P, I, F) and conditions (p, h - for I), cIi = 0
cP = 0.3
cIh = 0.5   # c_I hosp
cI = 0.6    # c_I home
cF = 1      # relative contagiousness after death before being buried

# initial number of sick (in state P) individuals
P0 = 25

# countermeasures
# day when isolation starts
t_iso = 120

# general reduction of contact with sick
fc = 1 # 0.2

# fraction at home (p) or in hospital (h) before this day
f_p1 = 0.5
f_h1 = 0.5

# fraction at home, in hospital, in isolation after this day (4 scenarios provided)
f_p2 = [0.5, 0.5, 0.5, 0.3, 0.1, 0.1]
f_h2 = [0.5, 0.3, 0.1, 0.1, 0.1, 0]

# index of the scenario from the ones above
k = 4

# probability that a person who had contact with somebody who gets into isolation (f_iso = 1- (f_p2 + f_h2)) is traced back
f_tb = 0.8

# goodness of isolation beyond isolation wards (0 - not better than else at home/hospital)
ph = 0.8

# capacity of isolation wards
qmax = N

# capacity for trace back
cmax = N/100

# probability of safe funeral for those dying
d_p = 0.5  # at home
d_h = 0.75  # in hospital (no isolation)

# vaccination
# start date
t_vac=days  # 500

# daily number
Nvac=0


