# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:13:23 2020

@author: helle
"""

from parameters_3 import *
import numpy as np
from scipy.integrate import solve_ivp

# index of compartments (depending on number of Erlang stages per compartment - assuming one general number of Erlang stages everywhere
def indexFunction(n):
    compartments = [('S',  0, '__'),
                ('E', n, '__'), ('E',n, 't_'),('E', n, 's_'),
                ('P', n, '__'), ('P',n, 't_'),('P', n, 's_'),
                ('I', n, '_p'), ('I',n, '_h'),('I', n, '_i'), ('I', n, 'sh'), ('I',n, 'sp'),
                ('F',  0, '__'),
                ('B',  0, '_j'),('B',  0, '_f'),
                ('R',  0, '__')
                ]
    index = dict()
    ind = 0
    for i in compartments:
        notation = i[0] + i[2]
        if i[1] == 0:
            index[notation] = ind
            ind += 1
        else:
            for k in np.arange(1, i[1] + 1):
                notation_k = notation + str(k)
                index[notation_k] = ind
                ind += 1
    #print(index)
    return index



# sum of population in compartment_script1_script1
def popsum(pop, compartment, script1, script2, n, index):
    x=0
    if compartment in {'S', 'F', 'B', 'R'}:
        x = x + pop[index[compartment + script1 + script2]]
    else:
        for k in range(1, n + 1):
            x = x + pop[index[compartment + script1 + script2 + str(k)]]
    return x

# general contact reduction, dependent on time (countermeasures being in place after t_iso)
def fct(t, t_iso, fc):
    if t < t_iso:
        fct = 1
    else:
        fct = fc
    return fct

# safe funeral, dependent on time  (countermeasures being in place after t_iso)
def d_ph(t, t_iso, d_ph0, d_ph1):
    if t <= t_iso:
        d_p = d_ph0[0]
        d_h = d_ph0[1]
    else:
        d_p = d_ph1[0]
        d_h = d_ph1[1]
    x = [d_p, d_h]
    return x

# fraction in isolation - and elsewhere -  (countermeasures being in place after t_iso)
def f_phi(t, t_iso, f_ph0, f_ph1):
    if t <= t_iso:
        f_p = f_ph0[0]
        f_h = f_ph0[1]
    else:
        f_p = f_ph1[0]
        f_h = f_ph1[1]
    f_i = 1 - (f_p + f_h)
    x = [f_p, f_h, f_i]
    return x

# persons to be in quarantine: potential (q) and actual (Q) (all in compartments Pt_ and I_i; after t_iso; limited by capacity qmax)
def q(pop, t, t_iso, qmax, n, index):
    if t < t_iso:
        q = 0
        Q = 0
    else:
        Q = (popsum(pop=pop, compartment='P', script1='t', script2='_', n=n, index=index) + \
             popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index))
        if Q <= qmax:
            q = 1
        else:
            q = qmax / Q
    return [q,Q]

# persons to be traced back: potential (c) and actual (C)
# only after t_iso
# depends on number of people entering status Pt or Ii
def c(pop, t, t_iso, cmax, n, index, FT, FP, f_iso):
    if t < t_iso:
        c = 0
        C = 0
    else:
        C = 1 / FT * (
                FP * (
                f_iso * (
                pop[index['P__' + str(n)]] + pop[index['Ps_' + str(n)]]) +
                pop[index['Pt_' + str(n)]])) + \
            popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index) + \
            popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index)
        if C <= cmax:
            c = 1
        else:
            c = cmax / C
    return [c, C]

# force of infection
def la(pop, fiso, f_tb, betaP, betaIp, betaIh, betaF, ph, q, c, fc, n, index):
    # infections that may be found later by tracing back
    ls0 = fc * betaP * (fiso * \
                        popsum(pop=pop, compartment='P', script1='_', script2='_', n=n, index=index) + \
                        popsum(pop=pop, compartment='P', script1='s', script2='_', n=n, index=index)) + \
          betaIp * \
          popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index) + \
          betaIh * \
          popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index) + \
          betaIh * (1 - ph) * (1 - q) * (
                  popsum(pop=pop, compartment='P', script1='t', script2='_', n=n, index=index) + \
                  popsum(pop=pop, compartment='I', script1='_', script2='i', n=n,
                         index=index))
    # infections that will not be traced back
    l = fc * betaP * (1 - fiso) * \
        popsum(pop=pop, compartment='P', script1='_', script2='_', n=n, index=index) + \
        betaIp * \
        popsum(pop=pop, compartment='I', script1='_', script2='p', n=n, index=index) + \
        betaIh * \
        popsum(pop=pop, compartment='I', script1='_', script2='h', n=n, index=index) + \
        betaF * pop[index['F__']] + \
        (1 - f_tb * c) * ls0
    # infections that will be traced back
    ls = f_tb * c * ls0
    return [l, ls]


def lt(la, ls):
    lt = la + ls
    return lt


def vac(pop, index, t, t_vac, N_vac):
    if t <= t_vac:
        x = 0
    if t > t_vac:
        x = min(pop[index['S__']], N_vac)
    return x

# differential equations ---------------------------------------------------------------

def dS(pop, lt, N, index, vac):
    x = - (lt / N) * pop[index['S__']] - vac
    return x


def dE__1(pop, la, N, FE, index):
    x = (la / N) * pop[index['S__']] \
        - FE * pop[index['E__1']]
    return x


def dE__k(pop, j, FE, index):  # k2 <= j <= NE
    x = FE * pop[index['E__' + str(j - 1)]] \
        - FE * pop[index['E__' + str(j)]]
    return x


def dEs_1(pop, ls, N, FT, FE, index):
    x = (ls / N) * pop[index['S__']] \
        - (FT + FE) * pop[index['Es_1']]
    return x


def dEs_k(pop, j, FE, FT, index):  # k2 <= j <= NE
    x = FE * pop[index['Es_' + str(j - 1)]] \
        - (FT + FE) * pop[index['Es_' + str(j)]]
    return x


def dEt_1(pop, FT, FE, index):
    x = FT * pop[index['Es_1']] \
        - FE * pop[index['Et_1']]
    return x


def dEt_k(pop, j, FT, FE, index):  # k2 <= j <= NE
    x = FT * pop[index['Es_' + str(j)]] \
        + FE * pop[index['Et_' + str(j - 1)]] \
        - FE * pop[index['Et_' + str(j)]]
    return x


def dP__1(pop, FE, FP, n, index):
    x = FE * pop[index['E__' + str(n)]] \
        - FP * pop[index['P__1']]
    return x


def dP__k(pop, j, FP, index):  # k2 <= j <= NP
    x = FP * pop[index['P__' + str(j - 1)]] \
        - FP * pop[index['P__' + str(j)]]
    return x


def dPs_1(pop, FE, FT, FP, n, index):
    x = FE * pop[index['Es_' + str(n)]] \
        - (FT + FP) * pop[index['Ps_1']]
    return x


def dPs_k(pop, j, FP, FT, index):  # k2 <= j <= NP
    x = FP * pop[index['Ps_' + str(j - 1)]] \
        - (FT + FP) * pop[index['Ps_' + str(j)]]
    return x


def dPt_1(pop, FT, FE, FP, n, index):
    x = FT * pop[index['Ps_1']] \
        + FE * pop[index['Et_' + str(n)]] \
        - FP * pop[index['Pt_1']]
    return x


def dPt_k(pop, j, FT, FP, index):  # k2 <= j <= NP
    x = FT * pop[index['Ps_' + str(j)]] \
        + FP * pop[index['Pt_' + str(j - 1)]] \
        - FP * pop[index['Pt_' + str(j)]]
    return x


def dI_p1(pop, fp, FP, FI, n, index):
    x = fp * FP * pop[index['P__' + str(n)]] \
        - FI * pop[index['I_p1']]
    return x


def dI_pk(pop, j, FI, index):
    x = FI * pop[index['I_p' + str(j - 1)]] \
        - FI * pop[index['I_p' + str(j)]]
    return x


def dI_h1(pop, fh, FP, FI, n, index):
    x = fh * FP * pop[index['P__' + str(n)]] \
        - FI * pop[index['I_h1']]
    return x


def dI_hk(pop, j, FI, index):
    x = FI * pop[index['I_h' + str(j - 1)]] \
        - FI * pop[index['I_h' + str(j)]]
    return x


def dIsp1(pop, fp, FP, FT, FI, n, index):
    x = fp * FP * pop[index['Ps_' + str(n)]] \
        - (FT + FI) * pop[index['Isp1']]
    return x


def dIspk(pop, j, FI, FT, index):  # 2 <=j < NIs
    x = FI * pop[index['Isp' + str(j - 1)]] \
        - (FT + FI) * pop[index['Isp' + str(j)]]
    return x


def dIspNI(pop, FI, FT, n, index):
    x = FI * pop[index['Isp' + str(n - 1)]] \
        - FT * pop[index['Isp' + str(n)]]
    return x


def dIsh1(pop, fh, FP, FT, FI, n, index):
    x = fh * FP * pop[index['Ps_' + str(n)]] \
        - (FT + FI) * pop[index['Ish1']]
    return x


def dIshk(pop, j, FI, FT, index):  # 2 <=j < NIh
    x = FI * pop[index['Ish' + str(j - 1)]] \
        - (FT + FI) * pop[index['Ish' + str(j)]]
    return x


def dIshNI(pop, FI, FT, n, index):
    x = FI * pop[index['Ish' + str(n - 1)]] \
        - FT * pop[index['Ish' + str(n)]]
    return x


def dI_i1(pop, fi, FP, FT, FI, n, index):
    x = fi * FP * pop[index['P__' + str(n)]] \
        + fi * FP * pop[index['Ps_' + str(n)]] \
        + FP * pop[index['Pt_' + str(n)]] \
        + FT * pop[index['Isp1']] \
        + FT * pop[index['Ish1']] \
        - FI * pop[index['I_i1']]
    return x


def dI_ik(pop, j, FI, FT, index):
    x = FI * pop[index['I_i' + str(j - 1)]] \
        + FT * pop[index['Isp' + str(j)]] \
        + FT * pop[index['Ish' + str(j)]] \
        - FI * pop[index['I_i' + str(j)]]
    return x


def dR(pop, fdead_p, fdead_h, fdead_i, FI, FIp, n, index, vac):
    x = FIp * (1 - fdead_p) * pop[index['I_p' + str(n)]] \
        + FI * ((1 - fdead_h) * pop[index['I_h' + str(n)]] \
                + (1 - fdead_i) * pop[index['I_i' + str(n)]]) \
        + vac
    return x


def dF(pop, fdead_p, fdead_h, FI, FIp, FF, n, index, d_h, d_p):
    x = FIp * (1 - d_p) * fdead_p * pop[index['I_p' + str(n)]] \
        + FI * (1 - d_h) * fdead_h * pop[index['I_h' + str(n)]] \
        - FF * pop[index['F__']]
    return x


def dB_f(pop, FF, index):
    x = FF * pop[index['F__']]
    return x


def dB_j(pop, fdead_i, FI, FIp, n, index, d_h, d_p, fdead_h, fdead_p):
    x = FIp * d_p * fdead_p * pop[index['I_p' + str(n)]] \
        + FI * (fdead_i * pop[index['I_i' + str(n)]] \
                + d_h * fdead_h * pop[index['I_h' + str(n)]])
    # * (1-q(t)) * d_h
    return x


# model --------------------------------------------------------------------------
def modelEbola(
        days=days,
        n=n,
        N=N,
        D = D,
        fdead = fdead,
        R0=R0,
        cc = cc,
        P0=P0,  # P(0)
        t_iso=t_iso,
        I_iso = I_iso,
        f_ph0 = f_ph0,
        f_ph1 = f_ph1,
        d_ph0 = d_ph0,
        d_ph1 = d_ph1,
        f_tb=f_tb,
        ph=ph,
        qmax=qmax,
        cmax=cmax,
        fc = fc,
        t_vac=t_vac,
        N_vac=N_vac,
        nameIn='',
        pathOut='results'):
    name = str(days) + '_' \
            + str(n) + '_' \
            + str(N) + '_' \
            + str(D) + '_' \
            + str(fdead) + '_' \
            + str(R0) + '_' \
            + str(cc) + '_' \
            + str(P0) + '_' \
            + str(t_iso) + '_' \
            + str(I_iso) + '_' \
            + str(f_ph0) + '_' \
            + str(f_ph1) + '_' \
            + str(d_ph0) + '_' \
            + str(d_ph1) + '_' \
            + str(f_tb) + '_' \
            + str(ph) + '_' \
            + str(qmax) + '_' \
            + str(cmax) + '_' \
            + str(fc) + '_' \
            + str(t_vac) + '_' \
            + str(N_vac) + '_' \
            + nameIn

    print("'" + name + "',")

    # compute values that do not change by time (or population)
    index = indexFunction(n)

    #cD = R0 / (cP * DP + cI * DI + cF * DF)
    cD = R0 /(cc[0] * D[1] + cc[2] * D[3] + cc[3] * D[4])
    cD = float(cD)

    betaP =     cc[0] * cD
    betaIh =    cc[1] * cD
    betaIp =    cc[2] * cD
    betaF =     cc[3] * cD

    FE = n/D[0] # epsilon
    FP = n/D[1] # gamma
    FI = n/D[2] # delta
    FIp = n/D[3] # delta at home
    FF = 1/D[4] # phi
    FT = 1/D[5] # alpha

    ####################################################
    #rec = [[-10000 for i in np.arange(2)] for j in np.arange(days + 1)]
    rec = [[-10000 for i in np.arange(10 + len(index))] for j in np.arange(days + 1)]

    def f(t, pop):
        # Initialize
        out = [0 for i in np.arange(len(index))]

        ## Values that change by time or population but are not differential equations
        # compute once per time and then use constants

        # are countermeasures in place? after time t_iso or after number of cases in I_h + I_i >= I_iso
        t_iso_ = t_iso
        if (t < t_iso):
            if popsum(pop=pop, compartment ='I', script1='_', script2='h', n=n, index=index) + \
                    popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index) + \
                    popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index) + \
                    pop[index['F__']] + pop[index['B_j']] + pop[index['B_f']] + pop[index['R__']]\
                    >= I_iso:
                t_iso_ = t-1

        fdead_i_ = fdead[2]
        fdead_h_ = fdead[1]
        fdead_p_ = fdead[0]

        fc_ = fct(t=t, t_iso=t_iso_, fc=fc)
        f_phi_ = f_phi(t=t, t_iso=t_iso_, f_ph0=f_ph0, f_ph1=f_ph1)
        d_ph_ = d_ph(t=t, t_iso=t_iso_, d_ph0 = d_ph0, d_ph1 = d_ph1)
        qq = q(pop=pop, t=t, t_iso=t_iso_, qmax=qmax, n=n, index=index)
        q_ = qq[0]
        cc = c(pop=pop, t=t, t_iso=t_iso_, cmax=cmax, n=n, index=index, FT=FT, FP=FP, f_iso=f_phi_[2])
        c_ = cc[0]
        la__ = la(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph,
                  q=q_, c = c_, fc = fc_, n=n, index=index)
        la_ = la__[0] # lambda
        ls_ = la__[1] # lambda^*
        lt_ = la_ + ls_  #lt(la=la_, ls=ls_)

        vac_ = vac(pop=pop, index=index, t=t,t_vac=t_vac, N_vac=N_vac)
        rec[int(t)] = [t] + [t_iso_] + [q_] + [fc_] + f_phi_ + [c_] + la__ + [vac_] + np.ndarray.tolist(pop)

        ## Differential equations for compartments
        # Susceptible
        out[index['S__']] = dS(pop=pop, lt=lt_, N=N, index=index, vac=vac_)

        # Latent
        out[index['E__1']] = dE__1(pop=pop, la=la_, N=N, FE=FE, index=index)

        for i in range(2, n + 1):
            out[index['E__' + str(i)]] = dE__k(pop=pop, j=i, FE=FE, index=index)

        out[index['Es_1']] = dEs_1(pop=pop, ls=ls_, N=N, FT=FT, FE=FE, index=index)

        for i in range(2, n + 1):
            out[index['Es_' + str(i)]] = dEs_k(pop=pop, j=i, FE=FE, FT=FT, index=index)

        out[index['Et_1']] = dEt_1(pop=pop, FT=FT, FE=FE, index=index)

        for i in range(2, n + 1):
            out[index['Et_' + str(i)]] = dEt_k(pop=pop, j=i, FE=FE, FT=FT, index=index)

        # Podromal
        out[index['P__1']] = dP__1(pop=pop, FE=FE, FP=FP, n=n, index=index)

        for i in range(2, n + 1):
            out[index['P__' + str(i)]] = dP__k(pop=pop, j=i, FP=FP, index=index)

        out[index['Ps_1']] = dPs_1(pop=pop, FE=FE, FT=FT, FP=FP, n=n, index=index)

        for i in range(2, n + 1):
            out[index['Ps_' + str(i)]] = dPs_k(pop=pop, j=i, FP=FP, FT=FT, index=index)

        out[index['Pt_1']] = dPt_1(pop=pop, FT=FT, FE=FE, FP=FP, n=n, index=index)

        for i in range(2, n + 1):
            out[index['Pt_' + str(i)]] = dPt_k(pop=pop, j=i, FP=FP, FT=FT, index=index)

        # Fully Infectious
        # I_p home
        out[index['I_p1']] = dI_p1(pop=pop, fp=f_phi_[0], FP=FP, FI=FIp, n=n, index=index)

        for i in range(2, n + 1):
            out[index['I_p' + str(i)]] = dI_pk(pop=pop, j=i, FI=FIp, index=index)

        # I_h hosp
        out[index['I_h1']] = dI_h1(pop=pop, fh=f_phi_[1], FP=FP, FI=FI, n=n, index=index)

        for i in range(2, n + 1):
            out[index['I_h' + str(i)]] = dI_hk(pop=pop, j=i, FI=FI, index=index)

        # Isp * home
        out[index['Isp1']] = dIsp1(pop=pop, fp=f_phi_[0], FP=FP, FT=FT, FI = FI, n=n, index=index)
        for i in range(2, n):
            out[index['Isp' + str(i)]] = dIspk(pop=pop, j=i, FI=FI, FT=FT, index=index)
        if n > 1:
            out[index['Isp' + str(n)]] = dIspNI(pop=pop, FI=FI, FT=FT, n=n, index=index)

        # Ish * hosp
        out[index['Ish1']] = dIsh1(pop=pop, fh=f_phi_[1], FP=FP, FT=FT, FI=FI, n=n, index=index)

        for i in range(2, n):
            out[index['Ish' + str(i)]] = dIshk(pop=pop, j=i, FI=FI, FT=FT, index=index)
        if n > 1:
            out[index['Ish' + str(n)]] = dIshNI(pop=pop, FI=FI, FT=FT, n=n, index=index)

        # I_i iso
        out[index['I_i1']] = dI_i1(pop=pop, fi=f_phi_[2], FP=FP, FT=FT, FI=FI, n=n, index=index)

        for i in range(2, n + 1):
            out[index['I_i' + str(i)]] = dI_ik(pop=pop, j=i, FI=FI, FT=FT, index=index)

        # Recovered, Dead
        out[index['R__']] = dR(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, fdead_i=fdead_i_, FI=FI, FIp = FIp, n=n, index=index, vac=vac_)

        out[index['F__']] = dF(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, FI=FI, FIp=FIp, FF=FF, n=n, index=index, d_h=d_ph_[1], d_p=d_ph_[0])

        out[index['B_f']] = dB_f(pop=pop, FF=FF, index=index)

        out[index['B_j']] = dB_j(pop=pop, fdead_i=fdead_i_, FI=FI, FIp=FIp, n = n, index=index, d_h=d_ph_[1], d_p=d_ph_[0], fdead_h=fdead_h_, fdead_p=fdead_p_)

        # return
        return out

    ###################################################
    #                 Solve ODEs                      #
    ###################################################
    # Initial values
    pop0 = [0 for i in np.arange(len(index))]

    # All individuals are susceptible
    pop0[index['S__']] = N - P0
    pop0[index['P__1']] = P0

    soln = solve_ivp(f,
                     [0, days],
                     pop0,
                     method="RK45",
                     t_eval=np.arange(0, days),
                     dense_output=True)

    # print(ints)
    np.savetxt(pathOut + "/ebola_" + name + ".txt", soln.y)
    #print(rec)
    np.savetxt(pathOut + "/ebolaVar_" + name + ".txt", rec,fmt='%.5f')
    return (name)

