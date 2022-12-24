# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:13:23 2020

@author: Kristina B. Helle, Aliou Bouba, Kristan A. Schneider
"""

from parameters import *
from scipy.integrate import solve_ivp
#import numpy as np

# index of compartments (depending on number of Erlang stages per compartment - assuming one general number of Erlang stages everywhere
def indexFunction(n):
    compartments = [('S', 0, '__'),
                    ('E', n, '__'), ('E', n, 't_'), ('E', n, 's_'),
                    ('P', n, '__'), ('P', n, 't_'), ('P', n, 's_'),
                    ('I', n, '_p'), ('I', n, '_h'), ('I', n, '_i'), ('I', n, 'sh'), ('I', n, 'sp'),
                    ('F', 0, '__'),
                    ('B', 0, '_j'), ('B', 0, '_f'),
                    ('R', 0, '__')
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
    # print(index)
    return index


# sum of population in compartment_script1_script1
def popsum(pop, compartment, script1, script2, n, index):
    x = 0
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
        Q = (popsum(pop=pop, compartment='E', script1='t', script2='_', n=n, index=index) + \
             popsum(pop=pop, compartment='P', script1='t', script2='_', n=n, index=index) + \
             popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index))
        if Q <= qmax:
            q = 1
        else:
            q = qmax / Q
    return [q, Q]


# persons to be traced back: potential (c) and actual (C)
# only after t_iso
# depends on number of people entering status Pt or Ii
def c(pop, t, t_iso, cmax, n, index, FT, FP, FE, f_iso):
    if t < t_iso:
        c = 0
        C = 0
    else:
        C = FT * (1 / FE * pop[index['Et_' + str(n)]] + \
                  1 / FP * (pop[index['P__' + str(n)]] + pop[index['Ps_' + str(n)]]) + \
                  1 / FT * (popsum(pop=pop, compartment='P', script1='s', script2='_', n=n, index=index) + \
                            popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index) + \
                            popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index)))

        if C <= cmax:
            c = 1
        else:
            c = cmax / C
    return [c, C]


# force of infection
def la(pop, fiso, f_tb, betaP, betaIp, betaIh, betaF, ph, q, c, fc, n, index):
    l = betaP * ((1 - fiso * f_tb * c) * popsum(pop=pop, compartment='P', script1='_', script2='_', n=n, index=index) + \
                 (1 - f_tb * c) * (popsum(pop=pop, compartment='P', script1='s', script2='_', n=n, index=index) +
                                   (1 - ph) * (1 - q) * popsum(pop=pop, compartment='P', script1='t', script2='_', n=n,
                                                               index=index))) + \
        betaIp * (popsum(pop=pop, compartment='I', script1='_', script2='p', n=n, index=index) + \
                  (1 - f_tb * c) * popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index)) + \
        betaIh * (popsum(pop=pop, compartment='I', script1='_', script2='h', n=n, index=index) +
                  (1 - f_tb * c) * popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index)) + \
        betaIh * (1 - ph) * (1 - q) * (1 - f_tb * c) * \
        popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index) + \
        betaF * pop[index['F__']]

    ls = f_tb * c * (betaP * (fiso * popsum(pop=pop, compartment='P', script1='_', script2='_', n=n, index=index) + \
                              popsum(pop=pop, compartment='P', script1='s', script2='_', n=n, index=index) + \
                              (1 - ph) * (1 - q) * popsum(pop=pop, compartment='P', script1='t', script2='_', n=n,
                                                          index=index)) + \
                     betaIp * popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index) + \
                     betaIh * popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index) + \
                     betaIh * (1 - ph) * (1 - q) * popsum(pop=pop, compartment='I', script1='_', script2='i', n=n,
                                                          index=index))
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
        D=D,
        DT=DT,
        fdead=fdead,
        R0=R0,
        cc=cc,
        P0=P0,  # P(0)
        t_iso=t_iso,
        I_iso=I_iso,
        f_ph0=f_ph0,
        f_ph1=f_ph1,
        d_ph0=d_ph0,
        d_ph1=d_ph1,
        f_tb=f_tb,
        ph=ph,
        qmax=qmax,
        cmax=cmax,
        fc=fc,
        t_vac=t_vac,
        N_vac=N_vac,
        nameIn='',
        method="RK45",
        pathOut='results'):
    D[5] = DT  # trace back time
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
           + method + '_' \
           + nameIn

    print("'" + name + "',")

    # compute values that do not change by time (or population)
    index = indexFunction(n)

    # cD = R0 / (cP * DP + cI * DI + cF * DF)
    u = (1 - d_ph0[0]) * fdead[0] * f_ph0[0] + \
        (1 - d_ph0[1]) * fdead[1] * f_ph0[1]
    cD = R0 / (cc[0] * D[1] + \
               cc[1] * D[3] * f_ph0[0] + \
               cc[2] * D[3] * f_ph0[1] + \
               cc[3] * D[4] * u)
    cD = float(cD)
    betaP = cc[0] * cD
    betaIp = cc[1] * cD
    betaIh = cc[2] * cD
    betaF = cc[3] * cD

    FE = n / D[0]  # epsilon
    FP = n / D[1]  # gamma
    FI = n / D[2]  # delta
    FIp = n / D[3]  # delta at home
    FF = 1 / D[4]  # phi
    FT = 1 / D[5]  # alpha

    ####################################################
    # rec = [[-10000 for i in np.arange(2)] for j in np.arange(days + 1)]
    rec = [[-10000 for i in np.arange(10 + len(index))] for j in np.arange(days + 1)]

    def f(t, pop):
        # Initialize
        out = [0 for i in np.arange(len(index))]

        ## Values that change by time or population but are not differential equations
        # compute once per time and then use constants

        # are countermeasures in place? after time t_iso or after number of cases in I_h + I_i >= I_iso
        t_iso_ = t_iso
        if (t < t_iso):
            if popsum(pop=pop, compartment='I', script1='_', script2='h', n=n, index=index) + \
                    popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index) + \
                    popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index) + \
                    pop[index['F__']] + pop[index['B_j']] + pop[index['B_f']] + pop[index['R__']] \
                    >= I_iso:
                t_iso_ = t - 1

        fdead_i_ = fdead[2]
        fdead_h_ = fdead[1]
        fdead_p_ = fdead[0]

        fc_ = fct(t=t, t_iso=t_iso_, fc=fc)
        f_phi_ = f_phi(t=t, t_iso=t_iso_, f_ph0=f_ph0, f_ph1=f_ph1)
        d_ph_ = d_ph(t=t, t_iso=t_iso_, d_ph0=d_ph0, d_ph1=d_ph1)
        qq = q(pop=pop, t=t, t_iso=t_iso_, qmax=qmax, n=n, index=index)
        q_ = qq[0]
        ccc = c(pop=pop, t=t, t_iso=t_iso_, cmax=cmax, n=n, index=index, FT=FT, FP=FP, FE=FE, f_iso=f_phi_[2])
        c_ = ccc[0]
        la__ = la(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph,
                  q=q_, c=c_, fc=fc_, n=n, index=index)
        la_ = la__[0]  # lambda
        ls_ = la__[1]  # lambda^*
        lt_ = la_ + ls_  # lt(la=la_, ls=ls_)

        vac_ = vac(pop=pop, index=index, t=t, t_vac=t_vac, N_vac=N_vac)
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
        out[index['Isp1']] = dIsp1(pop=pop, fp=f_phi_[0], FP=FP, FT=FT, FI=FI, n=n, index=index)
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
        out[index['R__']] = dR(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, fdead_i=fdead_i_, FI=FI, FIp=FIp, n=n,
                               index=index, vac=vac_)

        out[index['F__']] = dF(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, FI=FI, FIp=FIp, FF=FF, n=n, index=index,
                               d_h=d_ph_[1], d_p=d_ph_[0])

        out[index['B_f']] = dB_f(pop=pop, FF=FF, index=index)

        out[index['B_j']] = dB_j(pop=pop, fdead_i=fdead_i_, FI=FI, FIp=FIp, n=n, index=index, d_h=d_ph_[1],
                                 d_p=d_ph_[0], fdead_h=fdead_h_, fdead_p=fdead_p_)

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
    '''
    pop0[index['E__1']] = P0/10 # 
    pop0[index['E__2']] = P0/10
    pop0[index['E__3']] = P0/10
    pop0[index['E__4']] = P0/10
    pop0[index['E__5']] = P0/10
    pop0[index['E__6']] = P0/10
    pop0[index['E__7']] = P0/10
    pop0[index['E__8']] = P0/10
    pop0[index['E__9']] = P0/10
    pop0[index['E__10']] = P0/10
    '''
    soln = solve_ivp(f,
                     [0, days],
                     pop0,
                     method=method,
                     t_eval=np.arange(0, days),
                     dense_output=True,
                     max_step=1)

    # print(ints)
    np.savetxt(pathOut + "/ebola_" + name + ".txt", soln.y)
    # print(rec)
    np.savetxt(pathOut + "/ebolaVar_" + name + ".txt", rec, fmt='%.5f')
    return (name)


############
# Plot     #
############
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 20:21:52 2020

@author: helle
"""

import matplotlib.pyplot as plt
import numpy as np
import csv as csv
from matplotlib.lines import Line2D
# import parameters_1

# from index import *
from functions2 import indexFunction

# CD6700
colsA0 = ["#000000", "#801980", "#0073B3", "#59B3E6", "#009980", "#E69900", "#CC6600", "#CD6700", "#0073B3", "grey"]
colsA1 = ["#801980", "#0073B3", "#009980", "#E69900", "#CC6600", "#CD6700", "#0073B3"]
colsA0 = ["#000000", "#0073B3", "#009980", "#E69900", "#CC6600", "#CD6700", "#0073B3", "grey"]
colsA = ["#000000", "#801980", "#59B3E6", "#009980", "#E69900", "#CC6600", "#CD6700", "#0073B3", "grey"]
colsA = ["#000000", "#801980", "#0073B3", "#009980", "#E69900", "#CC6600", "#CD6700", "#0073B3"]


#colsA = ["#000000", "#801980", "#59B3E6", "#009980", "#E69900",   "#0073B3", "#CC6600","#CD6700","grey"]
lstsB = ['-', (0, (1, 1)), (0, (1, 2)), (0, (2, 1)),
         (0, (2, 2)), (0, (3, 1)), (0, (3, 2)), (0, (2, 1, 1, 1)), (0, (2, 2, 1, 2)), (0, (3, 1, 1, 1)),
         (0, (3, 2, 1, 2)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 2, 1, 2, 1, 2)), (0, (3, 1, 1, 1, 1, 1, 1, 1)),
         (0, (3, 2, 1, 2, 1, 2, 1, 2))]
lstsC = ['-', (0, (1, 1)), (0, (3, 1)), (0, (1, 3)), (0, (3, 3)),
         (0, (1, 1, 1, 3)), (0, (1, 1, 3, 1)), (0, (1, 1, 3, 3)), (0, (1, 3, 3, 1)), (0, (1, 3, 3, 3)),
         (0, (3, 1, 3, 3))]
lstsD = ['-', '-', (0, (1, 1)), (0, (3, 1)), (0, (1, 3)), (0, (1, 1, 1, 3)), (0, (1, 1, 3, 1)), (0, (1, 1, 3, 3)),
         '-', ]

# lstsA = ['-', (0, (5, 1)), (0, (3, 1)), (0, (3, 3)), (0, (1, 3)), (0, (1, 1, 1, 3)), (0, (1, 1, 3, 1)), (0, (1, 1, 3, 3)), (0, (1, 1)),]

# 11;1113, 1311;1131, 3111;1133, 3311;13;3113;31;1333,3313;3133,3331;33


pathIn = 'results'
#pathOut = 'plots'
pathOut = 'C:/Users/helle/PycharmProjects/Ebola2/plots'

def popsum2d(pops, n):
    index = indexFunction(n)
    indexComp = [1, 1 + n, 1 + 2 * n, 1 + 3 * n,
                 1 + 4 * n, 1 + 5 * n, 1 + 6 * n,
                 1 + 7 * n, 1 + 8 * n, 1 + 9 * n,
                 1 + 10 * n, 1 + 11 * n]

    popSum = np.zeros((16, pops.shape[1]))
    popSum[0] = pops[index['S__']]
    popSum[1] = np.sum(pops[indexComp[0]:indexComp[1]], axis=0)  # E__
    popSum[2] = np.sum(pops[indexComp[1]:indexComp[2]], axis=0)  # Et_
    popSum[3] = np.sum(pops[indexComp[2]:indexComp[3]], axis=0)  # Es_
    popSum[4] = np.sum(pops[indexComp[3]:indexComp[4]], axis=0)  # P__
    popSum[5] = np.sum(pops[indexComp[4]:indexComp[5]], axis=0)  # Pt_
    popSum[6] = np.sum(pops[indexComp[5]:indexComp[6]], axis=0)  # Ps_
    popSum[7] = np.sum(pops[indexComp[6]:indexComp[7]], axis=0)  # I_p
    popSum[8] = np.sum(pops[indexComp[7]:indexComp[8]], axis=0)  # I_h
    popSum[9] = np.sum(pops[indexComp[8]:indexComp[9]], axis=0)  # I_i
    popSum[10] = np.sum(pops[indexComp[9]:indexComp[10]], axis=0)  # Ish
    popSum[11] = np.sum(pops[indexComp[10]:indexComp[11]], axis=0)  # Isp
    popSum[12] = pops[index['F__']]
    popSum[13] = pops[index['B_j']]
    popSum[14] = pops[index['B_f']]
    popSum[15] = pops[index['R__']]
    return popSum


def getQ(pathIn, name, i, days):
    path = pathIn + '/ebolaVar_' + name + '.txt'
    q_ = np.loadtxt(path)
    q__ = np.delete(q_, np.where(q_ == [-1.00000000e+04, -1.00000000e+04]), axis=0)
    q___ = np.interp(x=np.arange(days), xp=q__[:, 0], fp=q__[:, i])
    return q___


def plotEbolaParameters(names, savename, n, days, pathIn=pathIn, pathOut=pathOut, nplots=10, col=colsA):
    index = indexFunction(n)
    par_ji = getQ(pathIn=pathIn, name=names[1], i=1, days=days)
    par = [[-10000 for i in np.arange(2)] for j in np.arange(days)]
    for j in range(0, nplots):
        par[j] = np.empty(shape=[len(names), days])  # s[1]])
        for i in range(0, len(names)):
            par[j][i] = getQ(pathIn=pathIn, name=names[i], i=j, days=days)

    fig = plt.figure()
    fig.set_size_inches(12, 12)

    p = fig.add_subplot(331)
    for i in range(0, len(names)):
        p.plot(par[0][i], color=col[i])

    p = fig.add_subplot(332)
    for i in range(0, len(names)):
        p.plot(par[1][i], color=col[i])

    p = fig.add_subplot(333)
    for i in range(0, len(names)):
        p.plot(par[2][i], color=col[i])

    p = fig.add_subplot(334)
    for i in range(0, len(names)):
        p.plot(par[3][i], color=col[i])

    p = fig.add_subplot(335)
    for i in range(0, len(names)):
        p.plot(par[4][i], color=col[i])

    p = fig.add_subplot(336)
    for i in range(0, len(names)):
        p.plot(par[5][i], color=col[i])

    p = fig.add_subplot(337)
    for i in range(0, len(names)):
        p.plot(par[6][i], color=col[i])

    p = fig.add_subplot(338)
    for i in range(0, len(names)):
        p.plot(par[7][i], color=col[i])

    p = fig.add_subplot(339)
    for i in range(0, len(names)):
        p.plot(par[8][i], color=col[i])
    # plt.show()

legendsA = [' ', 'never traced back', 'not yet traced back', 'diagnosed or traced back',
            'diagnosed or traced back, not in ward', 'diagnosed or traced back, in ward', 'unsafely buried',
            'safely buried', 'total buried']

# day=0 plots full scenario
def plotEbolaAll(names: object, savename: object, lab: object, n, pathIn: object = pathIn, pathOut: object = pathOut,
                 col: object = colsA, lst=lstsD, leg=legendsA,
                 q_max: object = False,
                 tb: object = False,
                 sf: object = False,
                 legendout: object = True, legendlong=False, days=0) -> object:
    # -------- layout settings
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['mathtext.default'] = 'regular'  # new
    # legendout: legend of everything outside plots (if False: each plot has extra legend)
    if legendout:
        plt.rcParams['legend.fontsize'] = 12
        legend_lines = []
        legend_text = []

    '''
    # dependent on tb etc. choose the correct linestyles and legends (some have to be the same)
    if sf == False:
        lst[5] = lst[0]
        lst[6] = lst[0]
    if q_max == False:
        lst[4] = lst[0]
    if q_max == True:
        leg[3] = leg[5]
    if tb == False:
        lst[1] = lst[0]
        lst[2] = lst[0]
        lst[3] = lst[0]
    '''

    # -------- load data
    # determine number of days
    popSum0 = popsum2d(pops=np.loadtxt(pathIn + '/ebola_' + names[0] + '.txt'), n=n)
    s = np.shape(popSum0)

    # if q_max: generate factor of persons that do not fit into the wards, to mulitply with
    if q_max == True:
        q_ = np.empty(shape=[len(names), s[1]])
        for i in range(0, len(names)):
            q_[i] = getQ(pathIn=pathIn, name=names[i], i=1, days=s[1])

    # load simulation
    popSum = np.empty(shape=[len(names), s[0], s[1]])

    popSum[0] = popSum0
    for i in range(0, len(names)):
        pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
        popSum[i] = popsum2d(pops_i, n=n)

    if days != 0:
        popSum1 = np.empty(shape=[len(names), s[0], days])
        for i in range(0, len(names)):
            for j in range(0, n):
                popSum1[i][j] = popSum[i][j][0:days]
        popSum = popSum1

    # ------------- plot  (9 subplots)
    fig = plt.figure()
    fig.set_size_inches(12, 12)

    # Susceptible
    p = fig.add_subplot(331)
    for i in range(0, len(names)):
        # generate legend of this subplot
        p.plot(popSum[i][0], label=lab[i], color=col[i], )
        # generate legend for general plot (colors)
        legend_lines = legend_lines + [Line2D([0], [0], lw=1, color=col[i])]
        legend_text = legend_text + [lab[i]]
    legend_lines = legend_lines + [Line2D([0], [0], lw=1, color='white')]
    legend_text = legend_text + [' ']
    if legendout == False:
        p.legend()
    p.set_ylabel('Susceptible ind.')
    p.set_title(label='A.', loc='left')
    p.set_ylim(bottom=0)

    # Latent
    p = fig.add_subplot(332)
    # generate legend of this subplot (linestyles)

    if tb == False:
        p.plot(popSum[0][1], color=col[0], linestyle=lst[1])
    if tb == True:
        p.plot(popSum[0][1], color=col[0], linestyle=lst[1], label=leg[1])  # 'never traced back'
        p.plot(popSum[0][2], color=col[0], linestyle=lst[2], label=leg[2])  # 'traced back'
        p.plot(popSum[0][3], color=col[0], linestyle=lst[3], label=leg[3])  # 'not yet traced back'

    for i in range(0, len(names)):
        p.plot(popSum[i][1], label=lab[i], color=col[i], linestyle=lst[1])
        if tb == True:
            p.plot(popSum[i][2], color=col[i], linestyle=lst[2])
            p.plot(popSum[i][3], color=col[i], linestyle=lst[3])
    if legendout == False:
        p.legend()
    p.set_ylabel('Latent ind.')
    p.set_title(label='B.', loc='left')
    p.set_ylim(bottom=0)

    # Podromal
    p = fig.add_subplot(333)
    if tb == False:
        p.plot(popSum[0][4], color=col[0])
    if tb == True:
        p.plot(popSum[0][4], color=col[0], linestyle=lst[1], label=leg[1])  # 'never traced back'
        p.plot(popSum[0][6], color=col[0], linestyle=lst[2], label=leg[2])  # 'not yet traced back'
        # generate general legend (linestyles)
        if legendlong == True:
            legend_lines = legend_lines \
                           + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[1])] \
                           + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[2])]
            legend_text = legend_text + [leg[1]] + [leg[2]]
        if q_max == False:
            p.plot(popSum[0][5], color=col[0], linestyle=lst[3], label=leg[3])  # 'traced back'
            if legendlong == True:
                legend_lines = legend_lines \
                               + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[3])]
                legend_text = legend_text + [leg[3]]
        if q_max == True:
            p.plot(np.multiply(popSum[0][5], q_[0]), color=col[0], linestyle=lst[3],
                   label=leg[3])  # 'traced back, in ward'
            p.plot(np.multiply(popSum[0][5], np.multiply(q_[0], -1) + 1), color=col[0], linestyle=lst[4],
                   label=leg[4])  # 'traced back, not in ward'
            if legendlong == True:
                legend_lines = legend_lines \
                               + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[3])] \
                               + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[4])]
                legend_text = legend_text + [leg[3]] + [leg[4]]
    for i in range(0, len(names)):
        if tb == False:
            p.plot(popSum[i][4], label=lab[i], color=col[i])
        if tb == True:
            p.plot(popSum[i][4], color=col[i], linestyle=lst[1], label=lab[i])
            p.plot(popSum[i][6], color=col[i], linestyle=lst[2])
            if q_max == False:
                p.plot(popSum[i][5], color=col[i], linestyle=lst[3])
            if q_max == True:
                p.plot(np.multiply(popSum[i][5], q_[i]), color=col[i], linestyle=lst[3])
                p.plot(np.multiply(popSum[i][5], np.multiply(q_[i], -1) + 1), color=col[i], linestyle=lst[4])

    if legendout == False:
        p.legend()
    # else:
    # p.legend(bbox_to_anchor=(-3, -3, 3.5, 0.5), loc="upper left", mode="expand", ncol=3)
    p.set_ylabel('Podromal ind.')
    p.set_title(label='C.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected at home
    p = fig.add_subplot(334)
    if tb == False:
        p.plot(popSum[0][7], label=lab[i], color=col[i])
    if tb == True:
        p.plot(popSum[0][7], color=col[0], linestyle=lst[1], label=leg[1])  # 'never traced back'
        p.plot(popSum[0][11], color=col[0], linestyle=lst[2], label=leg[2])  # 'not yet traced back'
    for i in range(0, len(names)):
        p.plot(popSum[i][7], label=lab[i], color=col[i], linestyle=lst[1])
        if tb == True:
            p.plot(popSum[i][11], color=col[i], linestyle=lst[2])
    if legendout == False:
        p.legend()
    p.set_ylabel('Fully inf. ind. at home')
    p.set_title(label='D.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in hospital
    p = fig.add_subplot(335)
    if tb == False:
        p.plot(popSum[0][8], color=col[0], linestyle=lst[1])
    if tb == True:
        p.plot(popSum[0][8], color=col[0], linestyle=lst[1], label=leg[1])  # 'never traced back'
        p.plot(popSum[0][10], color=col[0], linestyle=lst[2], label=leg[2])  # 'not yet traced back'
    for i in range(0, len(names)):
        p.plot(popSum[i][8], label=lab[i], color=col[i], linestyle=lst[1])
        if tb == True:
            p.plot(popSum[i][10], color=col[i], linestyle=lst[2])
    if legendout == False:
        p.legend()
    p.set_ylabel('Fully inf. ind. in hospital')
    p.set_title(label='E.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in isolation
    p = fig.add_subplot(336)
    if q_max == False:
        p.plot(popSum[0][9], color=col[0], linestyle=lst[3], label=leg[3])
    if q_max == True:
        p.plot(np.multiply(popSum[0][9], q_[0]), color=col[0], linestyle=lst[3], label=leg[3])  # 'in ward'
        p.plot(np.multiply(popSum[0][9], np.multiply(q_[0], -1) + 1), color=col[0], linestyle=lst[4],
               label=leg[4])  # 'not in ward'

    for i in range(0, len(names)):
        if q_max == False:
            p.plot(popSum[i][9], label=lab[i], color=col[i], linestyle=lst[3])
        if q_max == True:
            p.plot(np.multiply(popSum[i][9], q_[i]), color=col[i], linestyle=lst[3])
            p.plot(np.multiply(popSum[i][9], np.multiply(q_[i], -1) + 1), color=col[i], linestyle=lst[4])
    if legendout == False:
        p.legend()
    p.set_ylabel('Fully inf. ind. in isolation')
    p.set_title(label='F.', loc='left')
    p.set_ylim(bottom=0)

    # Unsafe funerals
    p = fig.add_subplot(337)
    for i in range(0, len(names)):
        p.plot(popSum[i][12], label=lab[i], color=col[i], linestyle=lst[5])
    if legendout == False:
        p.legend()
    p.set_ylabel('Unsafe funerals')
    p.set_title(label='G.', loc='left')
    p.set_ylim(bottom=0)

    # Buried
    p = fig.add_subplot(338)
    if sf == False:
        p.plot(popSum[0][14], color=col[0], linestyle=lst[5])
    if sf == True:
        p.plot(popSum[0][13], color=col[0], linestyle=lst[6], label=leg[7])  # 'safely')
        #    p.plot(popSum[0][14], color=col[0], linestyle=lst[5], label=leg[6])  # 'unsafely')
        p.plot(popSum[0][13] + popSum[0][14], color=col[0], linestyle=lst[8], label=leg[8])  # 'total')
        if legendlong == True:
            legend_lines = legend_lines \
                           + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[6])] \
                           + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[8])]
            #                           + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[5])]
            #           legend_text = legend_text + [leg[7]] + [leg[6]]+ [leg[8]]
            legend_text = legend_text + [leg[7]] + [leg[8]]
    for i in range(0, len(names)):
        if sf == False:
            p.plot(popSum[i][14], label=lab[i], color=col[i], linestyle=lst[5])
        if sf == True:
            p.plot(popSum[i][13], color=col[i], linestyle=lst[6])
            p.plot(popSum[i][13] + popSum[i][14], color=col[i], linestyle=lst[8])  # 'total')
    if legendout == False:
        p.legend()
    p.set_ylabel('Buried ind.')
    p.set_title(label='H.', loc='left')
    p.set_ylim(bottom=0)

    # Recovered
    p = fig.add_subplot(339)
    for i in range(0, len(names)):
        p.plot(popSum[i][15], label=lab[i], color=col[i])
    if legendout == False:
        p.legend()
    p.set_ylabel('Recovered ind.')
    p.set_title(label='I.', loc='left')
    p.set_ylim(bottom=0)

    if legendout:
        # plot general legend
        if legendlong == True:
            fig.legend(legend_lines, legend_text, bbox_to_anchor=(0.1, 0.09), loc="upper left",
                       ncol=6, mode="expand")  # mode="expand", ncol=)
        else:
            fig.legend(legend_lines, legend_text, bbox_to_anchor=(0.1, 0.09), loc="upper left",
                       ncol=6, mode="expand")  # mode="expand", ncol=)
        # bbox_to_anchor=(-3, -3, 3.5, 0.5)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
    # plt.tight_layout()
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.savefig(pathOut + '/Ebola_all_' + savename + '.pdf', dpi=100)
    plt.savefig(fname = pathOut + '/Ebola_all_' + savename, dpi=100)
    plt.show()
    '''
    plt.plot(popSum[0][0] + popSum[0][15], label='healthy (S+R)', color=col[0], linestyle='-')
    plt.plot(np.sum(popSum[0][1:12], axis=0), label='infected (E+P+I)', color=col[0], linestyle='--')
    plt.plot(np.sum(popSum[0][12:15], axis=0), label='dead (D+B)', color=col[0], linestyle=':')
    for i in range(0, len(names)):
        # plt.plot(popSum[i][0] + popSum[i][15], label=lab[i],color=col[i], linestyle='-')
        # print('healthy ' + lab[i] + ': '+ str(round((popSum[i][0] + popSum[i][15])[-1])))
        plt.plot(np.sum(popSum[i][1:12], axis=0), color=col[i], linestyle='--')
        print('infected ' + lab[i] + ': ' + str(round(np.sum(popSum[i][1:12], axis=0)[-1])))
        plt.plot(np.sum(popSum[i][12:15], axis=0), color=col[i], linestyle=':')
        print('dead ' + lab[i] + ': ' + str(round(np.sum(popSum[i][12:15], axis=0)[-1])))
        print('---------')
        plt.ylabel('Individuals')
    plt.legend()
    # plt.title(savename)
    plt.savefig(pathOut + '/Ebola_mix_' + savename + '.pdf', dpi=100)
    plt.show()
    '''

    result = [[-10000 for i in np.arange(3 * 16 + 3)] for j in np.arange(len(names) + 1)]
    result[0] = ['S__final', 'S__max', 'S__whenmax',
                 'E__final', 'E__max', 'E__whenmax',
                 'Et_final', 'Et_max', 'Et_whenmax',
                 'Es_final', 'Es_max', 'Es_whenmax',
                 'P__final', 'P__max', 'P__whenmax',
                 'Pt_final', 'Pt_max', 'Pt_whenmax',
                 'Ps_final', 'Ps_max', 'Ps_whenmax',
                 'I_pfinal', 'I_pmax', 'I_pwhenmax',
                 'I_hfinal', 'I_hmax', 'I_hwhenmax',
                 'I_ifinal', 'I_imax', 'I_iwhenmax',
                 'Ishfinal', 'Ishmax', 'Ishwhenmax',
                 'Ispfinal', 'Ispmax', 'Ispwhenmax',
                 'F__final', 'F__max', 'F__whenmax',
                 'B_jfinal', 'B_jmax', 'B_jwhenmax',
                 'B_ffinal', 'B_fmax', 'B_fwhenmax',
                 'R__final', 'R__max', 'R__whenmax',
                 'Iallfinal', 'Iallmax', 'Iallwhenmax']
    for i in np.arange(len(names)):
        result_i = []
        for j in np.arange(16):
            final_ij = popSum[i][j][-1]
            max_ij = max(popSum[i][j])
            maxwhen_ij = np.argmax(popSum[i][j])
            result_i = result_i + [final_ij] + [max_ij] + [maxwhen_ij]
        I_i = np.sum(popSum[i][1:11], axis=0)
        result_i = result_i + [I_i[-1]] + \
                   [max(I_i)] + \
                   [np.argmax(I_i)]

        print(result_i)
        result[i + 1] = result_i
    with open(pathOut + "/ebolaFinal_" + savename + ".csv", "w+") as my_csv:  # writing the file as my_csv
        csvWriter = csv.writer(my_csv, delimiter=',')  # using the csv module to write the file
        csvWriter.writerows(result)
    # np.savetxt(pathOut + "/ebolaFinal_" + savename + ".txt", result, fmt='%.5f')

    return (savename)



legendsA = ['total',
            'traced back at some time (B-E) / safely buried (H)',
            'diagnosed or traced back already',
            'diagnosed or traced back, in ward']

lstsA = ['-', (0, (3, 1)), (0, (1, 1)), (0, (1, 1, 3, 1))]


def plotEbolaAllCumulative(names: object, savename: object, lab: object, n, pathIn: object = pathIn, pathOut: object = pathOut,
                 col: object = colsA, lst=lstsA, leg=legendsA,
                 #            q_max: object = False,
                 tb: object = False,
                 sf: object = False,
                 legendout: object = True, legendlong=False, days=0, ncolsLeg = 6) -> object:
    # -------- layout settings
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['mathtext.default'] = 'regular'  # new
    # legendout: legend of everything outside plots (if False: each plot has extra legend)
    if legendout:
        plt.rcParams['legend.fontsize'] = 12
    legend_lines = []
    legend_text = []

    # -------- load data
    # determine number of days
    popSum0 = popsum2d(pops=np.loadtxt(pathIn + '/ebola_' + names[0] + '.txt'), n=n)
    s = np.shape(popSum0)

    # if q_max: generate factor of persons that do not fit into the wards, to mulitply with
    # if q_max == True:
    #    q_ = np.empty(shape=[len(names), s[1]])
    #    for i in range(0, len(names)):
    #        q_[i] = getQ(pathIn=pathIn, name=names[i], i=1, days=s[1])

    # load simulation
    popSum = np.empty(shape=[len(names), s[0], s[1]])

    popSum[0] = popSum0
    for i in range(0, len(names)):
        pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
        popSum[i] = popsum2d(pops_i, n=n)

    if days != 0:
        popSum1 = np.empty(shape=[len(names), s[0], days])
        for i in range(0, len(names)):
            for j in range(0, n):
                popSum1[i][j] = popSum[i][j][0:days]
        popSum = popSum1

    # ------------- plot  (9 subplots)
    fig = plt.figure()
    fig.set_size_inches(12, 12)

    # Susceptible
    p = fig.add_subplot(331)
    for i in range(0, len(names)):
        # generate legend of this subplot
        p.plot(popSum[i][0], label=lab[i], color=col[i], )
        # generate legend for general plot (colors)
        legend_lines = legend_lines + [Line2D([0], [0], lw=1, color=col[i])]
        legend_text = legend_text + [lab[i]]
    legend_lines = legend_lines + [Line2D([0], [0], lw=1, color='white')]
    legend_text = legend_text + [' ']
    if legendout == False:
        p.legend()
    p.set_ylabel('Susceptible ind.')
    p.set_title(label='A.', loc='left')
    p.set_ylim(bottom=0)

    # Latent
    p = fig.add_subplot(332)
    # generate legend of this subplot (linestyles)

    p.plot(popSum[0][1] + popSum[0][2] + popSum[0][3], color=col[0], linestyle=lst[0], label=leg[0])  # total
    if tb == True:
        p.plot(popSum[0][2] + popSum[0][3], color=col[0], linestyle=lst[1],
               label=leg[1])  # 'will be traced back at some time'
        p.plot(popSum[0][2], color=col[0], linestyle=lst[2], label=leg[2])  # 'traced back already'

    for i in range(0, len(names)):
        p.plot(popSum[i][1] + popSum[i][2] + popSum[i][3], color=col[i], linestyle=lst[0], label=lab[i])
        if tb == True:
            p.plot(popSum[i][2] + popSum[i][3], color=col[i], linestyle=lst[1])
            p.plot(popSum[i][2], color=col[i], linestyle=lst[2])
    if legendout == False:
        p.legend()
    p.set_ylabel('Latent ind.')
    p.set_title(label='B.', loc='left')
    p.set_ylim(bottom=0)

    # Podromal
    p = fig.add_subplot(333)

    p.plot(popSum[0][4] + popSum[0][5] + popSum[0][6], color=col[0], linestyle=lst[0])  # total
    if tb == True:
        p.plot(popSum[0][5] + popSum[0][6], color=col[0], linestyle=lst[1])  # 'will be traced back at some time'
        p.plot(popSum[0][5], color=col[0], linestyle=lst[2])  # traced back already
    #    if q_max == True:
    #        p.plot(np.multiply(popSum[0][5], q_[0]), color=col[0], linestyle=lst[3], label=leg[3])# 'traced back, in ward'

    for i in range(0, len(names)):
        p.plot(popSum[i][4] + popSum[i][5] + popSum[i][6], color=col[i], linestyle=lst[0])
        if tb == True:
            p.plot(popSum[i][5] + popSum[i][6], color=col[i], linestyle=lst[1])  # 'will be traced back at some time'
            p.plot(popSum[i][5], color=col[i], linestyle=lst[2])  # traced back already
    #        if q_max == True:
    #            p.plot(np.multiply(popSum[i][5], q_[0]), color=col[i], linestyle=lst[3], label=leg[3])  # 'traced back, in ward'
    if legendout == False:
        p.legend()
    # else:
    # p.legend(bbox_to_anchor=(-3, -3, 3.5, 0.5), loc="upper left", mode="expand", ncol=3)
    p.set_ylabel('Podromal ind.')
    p.set_title(label='C.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected at home
    p = fig.add_subplot(334)
    p.plot(popSum[0][7] + popSum[0][11], color=col[0], linestyle=lst[0])
    if tb == True:
        p.plot(popSum[0][11], color=col[0], linestyle=lst[1])  # 'not yet traced back'
    for i in range(0, len(names)):
        p.plot(popSum[i][7] + popSum[i][11], color=col[i], linestyle=lst[0])
        if tb == True:
            p.plot(popSum[i][11], color=col[i], linestyle=lst[1])  # 'not yet traced back'
    if legendout == False:
        p.legend()
    p.set_ylabel('Fully inf. ind. at home')
    p.set_title(label='D.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in hospital
    p = fig.add_subplot(335)
    p.plot(popSum[0][8] + popSum[0][10], color=col[0], linestyle=lst[0])
    if tb == True:
        p.plot(popSum[0][10], color=col[0], linestyle=lst[1])  # 'not yet traced back'
    for i in range(0, len(names)):
        p.plot(popSum[i][8] + popSum[i][10], label=lab[i], color=col[i], linestyle=lst[0])
        if tb == True:
            p.plot(popSum[i][10], color=col[i], linestyle=lst[1])
    if legendout == False:
        p.legend()
    p.set_ylabel('Fully inf. ind. in hospital')
    p.set_title(label='E.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in isolation
    p = fig.add_subplot(336)
    # if q_max == False:
    p.plot(popSum[0][9], color=col[0], linestyle=lst[2], label=leg[2])
    # if q_max == True:
    #    p.plot(np.multiply(popSum[0][9], q_[0]), color=col[0], linestyle=lst[3], label=leg[3])  # 'in ward'
    for i in range(0, len(names)):
        p.plot(popSum[i][9], color=col[i], linestyle=lst[2])
    #    if q_max == True:
    #        p.plot(np.multiply(popSum[i][9], q_[0]), color=col[i], linestyle=lst[3])  # 'in ward'
    if legendout == False:
        p.legend()
    p.set_ylabel('Fully inf. ind. in isolation')
    p.set_title(label='F.', loc='left')
    p.set_ylim(bottom=0)

    # Unsafe funerals
    p = fig.add_subplot(337)
    for i in range(0, len(names)):
        p.plot(popSum[i][12], label=lab[i], color=col[i], linestyle=lst[0])
    if legendout == False:
        p.legend()
    p.set_ylabel('Unsafe funerals')
    p.set_title(label='G.', loc='left')
    p.set_ylim(bottom=0)

    # Buried
    p = fig.add_subplot(338)
    p.plot(popSum[0][13] + popSum[0][14], color=col[0], linestyle=lst[0], label=leg[0])
    if sf == True:
        p.plot(popSum[0][13], color=col[0], linestyle=lst[1], label=leg[1])  # 'safely')
    for i in range(0, len(names)):
        p.plot(popSum[i][13] + popSum[i][14], color=col[i], linestyle=lst[0])
        if sf == True:
            p.plot(popSum[i][13], color=col[i], linestyle=lst[1])  # 'safely')
    if legendout == False:
        p.legend()
    p.set_ylabel('Buried ind.')
    p.set_title(label='H.', loc='left')
    p.set_ylim(bottom=0)

    # Recovered
    p = fig.add_subplot(339)
    for i in range(0, len(names)):
        p.plot(popSum[i][15], label=lab[i], color=col[i], linestyle=lst[0])
    if legendout == False:
        p.legend()
    p.set_ylabel('Recovered ind.')
    p.set_title(label='I.', loc='left')
    p.set_ylim(bottom=0)

    # --------- Legend
    if legendout:
        if legendlong == True:
            legend_lines = legend_lines \
                           + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[0])]
            legend_text = legend_text + [leg[0]]
            if tb == True:
                legend_lines = legend_lines \
                               + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[1])] \
                               + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[2])]
                legend_text = legend_text + [leg[1]] + [leg[2]]
        #        if q_max == True:
        #            legend_lines = legend_lines \
        #                           + [Line2D([0], [0], lw=1, color=col[0], linestyle=lst[3])]
        #            legend_text = legend_text + [leg[3]]

        # plot general legend
        if legendlong == True:
            fig.legend(legend_lines, legend_text, bbox_to_anchor=(0.09, 0.09), loc="upper left",
                       ncol=ncolsLeg)  # mode="expand", ncol=)
        else:
            fig.legend(legend_lines, legend_text, bbox_to_anchor=(0.09, 0.09), loc="upper left",
                       ncol=ncolsLeg)  # mode="expand", ncol=)
        # bbox_to_anchor=(-3, -3, 3.5, 0.5)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
    # plt.tight_layout()
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    #plt.savefig(pathOut + '/Ebola_all_' + savename + '.jpg', dpi=100)
    plt.savefig(pathOut + '/Ebola_all_' + savename, dpi=100)
    plt.show()
    '''
    plt.plot(popSum[0][0] + popSum[0][15], label='healthy (S+R)', color=col[0], linestyle='-')
    plt.plot(np.sum(popSum[0][1:12], axis=0), label='infected (E+P+I)', color=col[0], linestyle='--')
    plt.plot(np.sum(popSum[0][12:15], axis=0), label='dead (D+B)', color=col[0], linestyle=':')
    for i in range(0, len(names)):
        # plt.plot(popSum[i][0] + popSum[i][15], label=lab[i],color=col[i], linestyle='-')
        # print('healthy ' + lab[i] + ': '+ str(round((popSum[i][0] + popSum[i][15])[-1])))
        plt.plot(np.sum(popSum[i][1:12], axis=0), color=col[i], linestyle='--')
        print('infected ' + lab[i] + ': ' + str(round(np.sum(popSum[i][1:12], axis=0)[-1])))
        plt.plot(np.sum(popSum[i][12:15], axis=0), color=col[i], linestyle=':')
        print('dead ' + lab[i] + ': ' + str(round(np.sum(popSum[i][12:15], axis=0)[-1])))
        print('---------')
        plt.ylabel('Individuals')
    plt.legend()
    # plt.title(savename)
    plt.savefig(pathOut + '/Ebola_mix_' + savename + '.pdf', dpi=100)
    plt.show()
    '''

    result = [[-10000 for i in np.arange(3 * 16 + 3)] for j in np.arange(len(names) + 1)]
    result[0] = ['S__final', 'S__max', 'S__whenmax',
                 'E__final', 'E__max', 'E__whenmax',
                 'Et_final', 'Et_max', 'Et_whenmax',
                 'Es_final', 'Es_max', 'Es_whenmax',
                 'P__final', 'P__max', 'P__whenmax',
                 'Pt_final', 'Pt_max', 'Pt_whenmax',
                 'Ps_final', 'Ps_max', 'Ps_whenmax',
                 'I_pfinal', 'I_pmax', 'I_pwhenmax',
                 'I_hfinal', 'I_hmax', 'I_hwhenmax',
                 'I_ifinal', 'I_imax', 'I_iwhenmax',
                 'Ishfinal', 'Ishmax', 'Ishwhenmax',
                 'Ispfinal', 'Ispmax', 'Ispwhenmax',
                 'F__final', 'F__max', 'F__whenmax',
                 'B_jfinal', 'B_jmax', 'B_jwhenmax',
                 'B_ffinal', 'B_fmax', 'B_fwhenmax',
                 'R__final', 'R__max', 'R__whenmax',
                 'Iallfinal', 'Iallmax', 'Iallwhenmax']
    for i in np.arange(len(names)):
        result_i = []
        for j in np.arange(16):
            final_ij = popSum[i][j][-1]
            max_ij = max(popSum[i][j])
            maxwhen_ij = np.argmax(popSum[i][j])
            result_i = result_i + [final_ij] + [max_ij] + [maxwhen_ij]
        I_i = np.sum(popSum[i][1:11], axis=0)
        result_i = result_i + [I_i[-1]] + \
                   [max(I_i)] + \
                   [np.argmax(I_i)]

        print(result_i)
        result[i + 1] = result_i
    with open(pathOut + "/ebolaFinal_" + savename + ".csv", "w+") as my_csv:  # writing the file as my_csv
        csvWriter = csv.writer(my_csv, delimiter=',')  # using the csv module to write the file
        csvWriter.writerows(result)
    # np.savetxt(pathOut + "/ebolaFinal_" + savename + ".txt", result, fmt='%.5f')

    return (savename)


def plotEbolaScenarios(namesAllScenarios, lab, savename, n=16, pathIn=pathIn, pathOut=pathOut, col=colsA):
    fig = plt.figure()
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['axes.labelsize'] = 12

    # M sets of scenarios
    M = len(namesAllScenarios)
    # fig.set_size_inches(12, M * 4)
    fig.set_size_inches(12, 12)

    alphabet = []
    for letter in range(65, 91):
        alphabet.append(chr(letter))

    popSums = []
    for m in range(0, M):
        mm = 3 + (m - 1) * 3

        names = namesAllScenarios[m]
        popSum0 = popsum2d(pops=np.loadtxt(pathIn + '/ebola_' + names[0] + '.txt'), n=n)
        s = np.shape(popSum0)
        popSum = np.empty(shape=[len(names), s[0], s[1]])
        popSum[0] = popSum0

        for i in range(0, len(names)):
            pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
            popSum[i] = popsum2d(pops_i, n=n)

            # Result: final number of dead
            # print(names[i])
            print(np.sum(popSum[i][12:15], axis=0)[-1])

        # Susceptible
        p = fig.add_subplot(M, 3, mm + 1)

        for i in range(0, len(names)):
            p.plot(popSum[i][0], label=r'$f_{iso} = $' + lab[i], color=col[i], linestyle='-')
            # p.plot(popSum[i][0], label=lab[i], color=col[i], linestyle='-')
        if m == 0:
            # p.legend(bbox_to_anchor=(-0.05, -4.9, 3.9, 1), loc="upper center", mode="expand", ncol=5)
            p.legend(bbox_to_anchor=(-0.05, -4.9, 3.9, 1), loc="upper center", mode="expand", ncol=5)

        p.set_ylabel('Susceptible ind.')
        p.set_title(label=alphabet[mm] + '.', loc='left')
        p.set_ylim(bottom=0)

        # Infected
        p = fig.add_subplot(M, 3, mm + 2)
        for i in range(0, len(names)):
            p.plot(np.sum(popSum[i][1:12], axis=0), label=lab[i], color=col[i], linestyle='-')
        p.set_ylabel('Infected ind.')
        p.set_title(label=alphabet[mm + 1] + '.', loc='left')
        p.set_ylim(bottom=0)

        # Dead
        p = fig.add_subplot(M, 3, mm + 3)
        for i in range(0, len(names)):
            p.plot(np.sum(popSum[i][12:15], axis=0), label=lab[i], color=col[i], linestyle='-')
        p.set_ylabel('Dead ind.')
        p.set_title(label=alphabet[mm + 2] + '.', loc='left')
        p.set_ylim(bottom=0)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.25)

    plt.savefig(pathOut + '/Ebola_scenarios_' + savename + '.pdf', dpi=100, bbox_inches='tight')
    plt.show()

    return (savename)
