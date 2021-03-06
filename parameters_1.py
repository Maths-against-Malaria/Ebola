days = 1000
NE = NP = NIp = NIh = NIi = 16

N = 10_000

DE = 7
DP = 5.5
DI = 5
DF = 1.5 # time death to be buried
DT = 9   # trace back time

fdead_p = 0.7 # f dead home
fdead_h = 0.5 # f dead hosp
fdead_i = 0.3 # f dead iso

R0 = 2      # basic reproduction number

cP = 0.6
cIp = 0.9   # c_I home
cIh = 0.8   # c_I hosp
cF = 1      # relative contagiousness after death before being buried       

P0 = 1     # P(0)

t_iso = 250

f_p1 = 0.5
f_h1 = 0.5

#after day 250
f_p2 = [0.5, 0.5, 0.5, 0.3, 0.1]
f_h2 = [0.5, 0.3, 0.1, 0.1, 0.1]
k = 0 # scenario

f_tb = 1# trace back
ph = 0 # goodness of isolation beyond isolation wards
qmax = N # capacity of isolation wards


