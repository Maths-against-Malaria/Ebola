# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:13:23 2020

@author: Kristina B. Helle, Aliou Bouba, Kristan A. Schneider
"""

from parameters import *
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import csv as csv

# ================ model =================================
# --------------- help functions for model ---------------
# index of compartments
# (depending on number of Erlang stages per compartment - assuming one general number of Erlang stages everywhere: n)
def indexFunction(n):
    # stages (no Erlang substages)
    # number at position [1] is number of Erlang stages
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
    # Erlang substages
    for i in compartments:
        notation = i[0] + i[2]
        # with no Erlang stages
        if i[1] == 0:
            index[notation] = ind
            ind += 1
        # with Erlang stages
        else:
            # generate sub-stages
            for k in np.arange(1, i[1] + 1):
                notation_k = notation + str(k)
                index[notation_k] = ind
                ind += 1
    return index


# sum of population in compartment: compartment_script1_script2
def popsum(pop, compartment, script1, script2, n, index):
    x = 0
    # for compartments with only 1 stage
    if compartment in {'S', 'F', 'B', 'R'}:
        x = x + pop[index[compartment + script1 + script2]]
    # for compartments with n Erlang-stages
    else:
        # sum up population in all related Erlang-stages
        for k in range(1, n + 1):
            x = x + pop[index[compartment + script1 + script2 + str(k)]]
    return x


# general contact reduction
def fc_t(t, t_iso, fc):
    if t < t_iso:
        fc_ = 1
    # (countermeasures being in place after t_iso)
    else:
        fc_ = fc
    return fc_


# fraction of safe funeral at home (d_p), in hospital (d_h)
def d_ph(t, t_iso, d_ph0, d_ph1):
    if t <= t_iso:
        d_p = d_ph0[0]
        d_h = d_ph0[1]
    # countermeasures being in place after t_iso
    else:
        d_p = d_ph1[0]
        d_h = d_ph1[1]
    x = [d_p, d_h]
    return x


# fraction at home (f_p), in hospital (f_h), in isolation (f_i)
def f_phi(t, t_iso, f_ph0, f_ph1):
    if t <= t_iso:
        f_p = f_ph0[0]
        f_h = f_ph0[1]
    # countermeasures being in place after t_iso
    else:
        f_p = f_ph1[0]
        f_h = f_ph1[1]
    f_i = 1 - (f_p + f_h)
    x = [f_p, f_h, f_i]
    return x


# persons to be in quarantine: potential (q) and actual (qq)
def q__(pop, t, t_iso, qmax, n, index):
    if t < t_iso:
        q = 0
        qq = 0
    # countermeasures being in place after t_iso
    else:
        # all in compartments Et_, Pt_, I_i
        qq = (popsum(pop=pop, compartment='E', script1='t', script2='_', n=n, index=index) +
              popsum(pop=pop, compartment='P', script1='t', script2='_', n=n, index=index) +
              popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index))
        # limited by capacity qmax
        if qq <= qmax:
            q = 1
        else:
            q = qmax / qq
    return [q, qq]


# persons to be traced back: potential (c) and actual (cc)
def c(pop, t, t_iso, cmax, n, index, FT, FP, FE, f_iso):
    if t < t_iso:
        c = 0
        cc = 0
    # countermeasures being in place after t_iso
    else:
        # potential
        # number of people entering status Pt_ or I_i
        # from last Erlang stage of Et_->Pt_, f_iso * ( P__->I_i, Ps_->I_i)
        # by rate 1/FT from all stages of Ps_-> Pt_, Isp->I_i, Ish->I_i
        # all staying in trace back for the average duration FT
        cc = FT * (1 / FE * pop[index['Et_' + str(n)]] + \
                  1 / FP * f_iso * (pop[index['P__' + str(n)]] + pop[index['Ps_' + str(n)]]) + \
                  1 / FT * (popsum(pop=pop, compartment='P', script1='s', script2='_', n=n, index=index) + \
                            popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index) + \
                            popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index)))
        # limited by cmax
        if cc <= cmax:
            c = 1
        else:
            c = cmax / cc
    return [c, cc]


# force of infection
# leading to infections subject to trace back (ls) or not (la)
# (1 - fiso * f_tb * c): not diagnosed, not found by trace back
def la(pop, fiso, f_tb, betaP, betaIp, betaIh, betaF, ph, q, c, fc, n, index):
    # not subject to trace back (l)
    # all modified by general contact reduction fc
    # P modified by betaP; I modified by betaIp, betaIh respectively; F modified by betaF
    # P__: except those who are diagnosed (fiso), and traced back (f_tb) within capacity (c)
    # Ps_: except those who traced back (f_tb) within capacity (c)
    # Pt_: only those exceeding isolation capacity (1-q), with reduced force (1-ph)
    # I_p all, Isp except those traced back (f_tb) within capacity (c)
    # I_h all, Ish except those traced back (f_tb) within capacity (c)
    # I_i only those exceeding isolation capacity (1-q), with reduced force (1-ph), except those traced back (f_tb) within capacity (c)
    # F all
    la = fc * (betaP * ((1 - fiso * f_tb * c) * popsum(pop=pop, compartment='P', script1='_', script2='_', n=n, index=index) + \
                 (1 - f_tb * c) * (popsum(pop=pop, compartment='P', script1='s', script2='_', n=n, index=index) +
                                   (1 - ph) * (1 - q) * popsum(pop=pop, compartment='P', script1='t', script2='_', n=n,
                                                               index=index))) + \
        betaIp * (popsum(pop=pop, compartment='I', script1='_', script2='p', n=n, index=index) + \
                  (1 - f_tb * c) * popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index)) + \
        betaIh * (popsum(pop=pop, compartment='I', script1='_', script2='h', n=n, index=index) +
                  (1 - f_tb * c) * popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index)) + \
        betaIh * (1 - ph) * (1 - q) * (1 - f_tb * c) * \
        popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index) + \
        betaF * pop[index['F__']])
    # subject to trace back (la)
    # all modified by general contact reduction fc, only those that are traced back (f_tb) within capacity (c)
    # P modified by betaP; I modified by betaIp, betaIh respectively
    # P__: those who are diagnosed (fiso)
    # Ps_: all
    # Pt_: only those exceeding isolation capacity (1-q), with reduced force (1-ph)
    # Isp: all
    # Ish: all
    # I_i only those exceeding isolation capacity (1-q), with reduced force (1-ph)
    ls = fc * f_tb * c * (betaP * (fiso * popsum(pop=pop, compartment='P', script1='_', script2='_', n=n, index=index) + \
                              popsum(pop=pop, compartment='P', script1='s', script2='_', n=n, index=index) + \
                              (1 - ph) * (1 - q) * popsum(pop=pop, compartment='P', script1='t', script2='_', n=n,
                                                          index=index)) + \
                     betaIp * popsum(pop=pop, compartment='I', script1='s', script2='p', n=n, index=index) + \
                     betaIh * popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index) + \
                     betaIh * (1 - ph) * (1 - q) * popsum(pop=pop, compartment='I', script1='_', script2='i', n=n,
                                                          index=index))
    return [la, ls]


# number of vaccinations
def vac(pop, index, t, t_vac, N_vac):
    if t <= t_vac:
        x = 0
    # being in place after t_vac
    # only applied to individuals in S__
    # limited by N_vac
    if t > t_vac:
        x = min(pop[index['S__']], N_vac)
    return x


# -------------- differential equations ------------

# suspectible
#           reduced by total force of infection
#           reduced by vaccination
def dS(pop, lt, N, index, vac):
    x = - (lt / N) * pop[index['S__']] - vac
    return x

# latent
#   not subject to trace back
#     1st Erlang stage
#           increased by infections (not subject to trace back)
#           decreased by rate
def dE__1(pop, la, N, FE, index):
    x = (la / N) * pop[index['S__']] \
        - FE * pop[index['E__1']]
    return x

#     2nd - last Erlang stage: 2 <= j <= NE
#           increased and decreased by rate
def dE__k(pop, j, FE, index):
    x = FE * pop[index['E__' + str(j - 1)]] \
        - FE * pop[index['E__' + str(j)]]
    return x

#   subject to trace back later
#     1st Erlang stage
#           increased by infections (subject to trace back)
#           decreased by rate
#           decreased by trace back
def dEs_1(pop, ls, N, FT, FE, index):
    x = (ls / N) * pop[index['S__']] \
        - (FT + FE) * pop[index['Es_1']]
    return x

#     2nd - last Erlang stage: 2 <= j <= NE
#           increased and decreased by rate
#           decreased by trace back
def dEs_k(pop, j, FE, FT, index):
    x = FE * pop[index['Es_' + str(j - 1)]] \
        - (FT + FE) * pop[index['Es_' + str(j)]]
    return x

#   found by trace back
#     1st Erlang stage
#           decreased by rate
#           increased by trace back
def dEt_1(pop, FT, FE, index):
    x = FT * pop[index['Es_1']] \
        - FE * pop[index['Et_1']]
    return x

#     2nd - last Erlang stage: 2 <= j <= NE
#           increased and decreased by rate
#           increased by trace back
def dEt_k(pop, j, FT, FE, index):
    x = FT * pop[index['Es_' + str(j)]] \
        + FE * pop[index['Et_' + str(j - 1)]] \
        - FE * pop[index['Et_' + str(j)]]
    return x

# podromal
#   not subject to trace back
#     1st Erlang stage
#           increased and decreased by rate
def dP__1(pop, FE, FP, n, index):
    x = FE * pop[index['E__' + str(n)]] \
        - FP * pop[index['P__1']]
    return x

#     2nd - last Erlang stage: 2 <= j <= NP
#           increased and decreased by rate
def dP__k(pop, j, FP, index):
    x = FP * pop[index['P__' + str(j - 1)]] \
        - FP * pop[index['P__' + str(j)]]
    return x

#   subject to trace back later
#     1st Erlang stage
#           increased and decreased by rate
#           decreased by trace back
def dPs_1(pop, FE, FT, FP, n, index):
    x = FE * pop[index['Es_' + str(n)]] \
        - (FT + FP) * pop[index['Ps_1']]
    return x

#     2nd - last Erlang stage: 2 <= j <= NP
#           increased and decreased by rate
#           decreased by trace back
def dPs_k(pop, j, FP, FT, index):
    x = FP * pop[index['Ps_' + str(j - 1)]] \
        - (FT + FP) * pop[index['Ps_' + str(j)]]
    return x

#   found by trace back
#     1st Erlang stage
#           increased and decreased by rate
#           decreased by trace back
def dPt_1(pop, FT, FE, FP, n, index):
    x = FT * pop[index['Ps_1']] \
        + FE * pop[index['Et_' + str(n)]] \
        - FP * pop[index['Pt_1']]
    return x

#     2nd - last Erlang stage: 2 <= j <= NP
#           increased and decreased by rate
#           increased by trace back
def dPt_k(pop, j, FT, FP, index):
    x = FT * pop[index['Ps_' + str(j)]] \
        + FP * pop[index['Pt_' + str(j - 1)]] \
        - FP * pop[index['Pt_' + str(j)]]
    return x

# fully infectious
#   not subject to trace back
#     at home
#       1st Erlang stage
#           increased by rate, fraction getting this treatment; decreased by rate
def dI_p1(pop, fp, FP, FI, n, index):
    x = fp * FP * pop[index['P__' + str(n)]] \
        - FI * pop[index['I_p1']]
    return x

#       2nd - last Erlang stage: 2 <= j <= NIp
#           increased and decreased by rate
def dI_pk(pop, j, FI, index):
    x = FI * pop[index['I_p' + str(j - 1)]] \
        - FI * pop[index['I_p' + str(j)]]
    return x

#     in hospital
#       1st Erlang stage
#           increased by rate, fraction getting this treatment; decreased by rate
def dI_h1(pop, fh, FP, FI, n, index):
    x = fh * FP * pop[index['P__' + str(n)]] \
        - FI * pop[index['I_h1']]
    return x

#       2nd - last Erlang stage: 2 <= j <= NIh
#           increased and decreased by rate
def dI_hk(pop, j, FI, index):
    x = FI * pop[index['I_h' + str(j - 1)]] \
        - FI * pop[index['I_h' + str(j)]]
    return x

#    subject to trace back later
#     at home
#       1st Erlang stage
#           increased by rate, fraction getting this treatment; decreased by rate
#           decreased by trace back
def dIsp1(pop, fp, FP, FT, FIp, n, index):
    x = fp * FP * pop[index['Ps_' + str(n)]] \
        - (FT + FIp) * pop[index['Isp1']]
    return x

#       2nd - second last Erlang stage: 2 <= j <= NIp-1
#           increased and decreased by rate
#           decreased by trace back
def dIspk(pop, j, FIp, FT, index):
    x = FIp * pop[index['Isp' + str(j - 1)]] \
        - (FT + FIp) * pop[index['Isp' + str(j)]]
    return x

#       last Erlang stage: NIp
#           increased by rate
#           decreased by trace back
def dIspNI(pop, FIp, FT, n, index):
    x = FIp * pop[index['Isp' + str(n - 1)]] \
        - FT * pop[index['Isp' + str(n)]]
    return x

#     in hospital
#       1st Erlang stage
#           increased by rate, fraction getting this treatment; decreased by rate
#           decreased by trace back
def dIsh1(pop, fh, FP, FT, FI, n, index):
    x = fh * FP * pop[index['Ps_' + str(n)]] \
        - (FT + FI) * pop[index['Ish1']]
    return x

#       2nd - second last Erlang stage: 2 <= j <= NIh-1
#           increased and decreased by rate
#           decreased by trace back
def dIshk(pop, j, FI, FT, index):
    x = FI * pop[index['Ish' + str(j - 1)]] \
        - (FT + FI) * pop[index['Ish' + str(j)]]
    return x

#       last Erlang stage: NIh
#           increased by rate
#           decreased by trace back
def dIshNI(pop, FI, FT, n, index):
    x = FI * pop[index['Ish' + str(n - 1)]] \
        - FT * pop[index['Ish' + str(n)]]
    return x

#     in isolation
#       1st Erlang stage
#           by P__, Ps_: increased by rate, fraction getting this treatment
#           by Pt_: increased by rate
#           by Isp1, ISh1: increased by trace back
#           decreased by rate
def dI_i1(pop, fi, FP, FT, FI, n, index):
    x = fi * FP * pop[index['P__' + str(n)]] \
        + fi * FP * pop[index['Ps_' + str(n)]] \
        + FP * pop[index['Pt_' + str(n)]] \
        + FT * pop[index['Isp1']] \
        + FT * pop[index['Ish1']] \
        - FI * pop[index['I_i1']]
    return x

#       2nd - last Erlang stage: 2 <= j <= NIi
#           increased and decreased by rate
#           increased by trace back
def dI_ik(pop, j, FI, FT, index):
    x = FI * pop[index['I_i' + str(j - 1)]] \
        + FT * pop[index['Isp' + str(j)]] \
        + FT * pop[index['Ish' + str(j)]] \
        - FI * pop[index['I_i' + str(j)]]
    return x

# recovered
# increased by rates, those who do not die
# increased by vaccinated
def dR(pop, fdead_p, fdead_h, fdead_i, FI, FIp, n, index, vac):
    x = FIp * (1 - fdead_p) * pop[index['I_p' + str(n)]] \
        + FI * ((1 - fdead_h) * pop[index['I_h' + str(n)]] \
                + (1 - fdead_i) * pop[index['I_i' + str(n)]]) \
        + vac
    return x

# dead and awaiting unsafe funeral
# increased by rates, those who die, and do not receive safe funeral
# decreased by rate
def dF(pop, fdead_p, fdead_h, FI, FIp, FF, n, index, d_h, d_p):
    x = FIp * (1 - d_p) * fdead_p * pop[index['I_p' + str(n)]] \
        + FI * (1 - d_h) * fdead_h * pop[index['I_h' + str(n)]] \
        - FF * pop[index['F__']]
    return x

# dead, unsafely buried
# increased by rate
def dB_f(pop, FF, index):
    x = FF * pop[index['F__']]
    return x

# dead, safely buried
# increased by rates, those who die, and receive safe funeral
def dB_j(pop, fdead_i, FI, FIp, n, index, d_h, d_p, fdead_h, fdead_p):
    x = FIp * d_p * fdead_p * pop[index['I_p' + str(n)]] +\
         FI * (fdead_i * pop[index['I_i' + str(n)]]
                + d_h * fdead_h * pop[index['I_h' + str(n)]])
    return x


# ---------------- model ------------------------------
# SEIR model (function solving ODE inside)
# for parameters, see parameters file
# pathResult: where to safe result files to, this folder has to exist
def modelEbola(
        days = days,
        n = n,
        N = N,
        D = D,
        DT = DT,
        fdead = fdead,
        R0 = R0,
        cc = cc,
        P0 = P0,  # P(0)
        t_iso = t_iso,
        I_iso = I_iso,
        f_ph0 = f_ph0,
        f_ph1 = f_ph1,
        d_ph0 = d_ph0,
        d_ph1 = d_ph1,
        f_tb = f_tb,
        ph = ph,
        qmax = qmax,
        cmax = cmax,
        fc = fc,
        t_vac = t_vac,
        N_vac = N_vac,
        pathResult=pathResult):
    D[5] = DT  # trace back time
    # ---------- create file name of all parameters -----------------
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
           + str(N_vac) + '_'


    # generate string with name that can be used to call result file for plotting
    print("'" + name + "',")

    # -------------------- values that do not change by time (or population) -----------------
    index = indexFunction(n)

    # contact rates for different stages
    #   common denominator for cD = R0 / (cP * DP + cI * DI + cF * DF)
    u = (1 - d_ph0[0]) * fdead[0] * f_ph0[0] + \
        (1 - d_ph0[1]) * fdead[1] * f_ph0[1]

    #   modified reproduction number for different stages
    cD = R0 / (cc[0] * D[1] +
               cc[1] * D[3] * f_ph0[0] +
               cc[2] * D[3] * f_ph0[1] +
               cc[3] * D[4] * u)
    cD = float(cD)

    # adapted contact rates
    betaP = cc[0] * cD
    betaIp = cc[1] * cD
    betaIh = cc[2] * cD
    betaF = cc[3] * cD

    # rates of disease progression
    FE = n / D[0]  # epsilon
    FP = n / D[1]  # gamma
    FI = n / D[2]  # delta in hospital and isolation
    FIp = n / D[3]  # delta at home
    FF = 1 / D[4]  # phi
    FT = 1 / D[5]  # alpha

    # --------------- time dependent values --------------------------
    # initialize object to save parameters and population in all compartments, at all times (days)
    rec = [[-10000 for i in np.arange(10 + len(index))] for j in np.arange(days + 1)]

    # function of population per compartment at all times
    def f(t, pop):
        # Initialize population in all compartments, at current time
        out = [0 for i in np.arange(len(index))]

        # Values that change by time or population but are not differential equations
        #   compute once per time and then use constants

        # are countermeasures in place at this time?
        t_iso_ = t_iso
        #   after time t_iso
        if (t < t_iso):
            #  or after number of cases exceeds I_iso
            if popsum(pop=pop, compartment='I', script1='_', script2='h', n=n, index=index) + \
                    popsum(pop=pop, compartment='I', script1='_', script2='i', n=n, index=index) + \
                    popsum(pop=pop, compartment='I', script1='s', script2='h', n=n, index=index) + \
                    pop[index['F__']] + pop[index['B_j']] + pop[index['B_f']] + pop[index['R__']] \
                    >= I_iso:
                t_iso_ = t - 1

        # adjust all parameters that depend on countermeasures being in place
        # death rates
        fdead_i_ = fdead[2]
        fdead_h_ = fdead[1]
        fdead_p_ = fdead[0]

        # general reduction of contacts
        fc_ = fc_t(t=t, t_iso=t_iso_, fc=fc)

        # fraction receiving different treatment
        f_phi_ = f_phi(t=t, t_iso=t_iso_, f_ph0=f_ph0, f_ph1=f_ph1)

        # fraction receiving safe funeral
        d_ph_ = d_ph(t=t, t_iso=t_iso_, d_ph0=d_ph0, d_ph1=d_ph1)

        # check, if isolation capacity is exceeded (to adjust force of infection)
        qq = q__(pop=pop, t=t, t_iso=t_iso_, qmax=qmax, n=n, index=index)
        q_ = qq[0]

        # check, if trace back capacity is exceeded  (to adjust force of infection)
        ccc = c(pop=pop, t=t, t_iso=t_iso_, cmax=cmax, n=n, index=index, FT=FT, FP=FP, FE=FE, f_iso=f_phi_[2])
        c_ = ccc[0]

        # force of infection
        la__ = la(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph,
                  q=q_, c=c_, fc=fc_, n=n, index=index)
        # not subject to trace back (lambda)
        la_ = la__[0]
        # subject to trace back (lambda^*)
        ls_ = la__[1]
        # total
        lt_ = la_ + ls_

        # number of current vaccinations
        vac_ = vac(pop=pop, index=index, t=t, t_vac=t_vac, N_vac=N_vac)


        # keep record of current parameters and population
        rec[int(t)] = [t] + [t_iso_] + [q_] + [fc_] + f_phi_ + [c_] + la__ + [vac_] + np.ndarray.tolist(pop)

        # ------------ update differential equations for compartments ----------
        # Susceptible
        out[index['S__']] = dS(pop=pop, lt=lt_, N=N, index=index, vac=vac_)

        # Latent
        #   not subject to trace back
        out[index['E__1']] = dE__1(pop=pop, la=la_, N=N, FE=FE, index=index)
        for i in range(2, n + 1):
            out[index['E__' + str(i)]] = dE__k(pop=pop, j=i, FE=FE, index=index)

        #   traced back later
        out[index['Es_1']] = dEs_1(pop=pop, ls=ls_, N=N, FT=FT, FE=FE, index=index)
        for i in range(2, n + 1):
            out[index['Es_' + str(i)]] = dEs_k(pop=pop, j=i, FE=FE, FT=FT, index=index)

        #   found by trace back
        out[index['Et_1']] = dEt_1(pop=pop, FT=FT, FE=FE, index=index)
        for i in range(2, n + 1):
            out[index['Et_' + str(i)]] = dEt_k(pop=pop, j=i, FE=FE, FT=FT, index=index)

        # Podromal
        #   not subject to trace back
        out[index['P__1']] = dP__1(pop=pop, FE=FE, FP=FP, n=n, index=index)
        for i in range(2, n + 1):
            out[index['P__' + str(i)]] = dP__k(pop=pop, j=i, FP=FP, index=index)

        #   traced back later
        out[index['Ps_1']] = dPs_1(pop=pop, FE=FE, FT=FT, FP=FP, n=n, index=index)
        for i in range(2, n + 1):
            out[index['Ps_' + str(i)]] = dPs_k(pop=pop, j=i, FP=FP, FT=FT, index=index)

        #   found by trace back
        out[index['Pt_1']] = dPt_1(pop=pop, FT=FT, FE=FE, FP=FP, n=n, index=index)
        for i in range(2, n + 1):
            out[index['Pt_' + str(i)]] = dPt_k(pop=pop, j=i, FP=FP, FT=FT, index=index)

        # Fully Infectious
        #   not subject to trace back
        #     at home
        out[index['I_p1']] = dI_p1(pop=pop, fp=f_phi_[0], FP=FP, FI=FIp, n=n, index=index)
        for i in range(2, n + 1):
            out[index['I_p' + str(i)]] = dI_pk(pop=pop, j=i, FI=FIp, index=index)

        #     in hospital
        out[index['I_h1']] = dI_h1(pop=pop, fh=f_phi_[1], FP=FP, FI=FI, n=n, index=index)
        for i in range(2, n + 1):
            out[index['I_h' + str(i)]] = dI_hk(pop=pop, j=i, FI=FI, index=index)

        #   traced back later
        #     at home
        out[index['Isp1']] = dIsp1(pop=pop, fp=f_phi_[0], FP=FP, FT=FT, FIp=FIp, n=n, index=index)
        for i in range(2, n):
            out[index['Isp' + str(i)]] = dIspk(pop=pop, j=i, FIp=FIp, FT=FT, index=index)
        if n > 1:
            out[index['Isp' + str(n)]] = dIspNI(pop=pop, FIp=FIp, FT=FT, n=n, index=index)

        #     in hospital
        out[index['Ish1']] = dIsh1(pop=pop, fh=f_phi_[1], FP=FP, FT=FT, FI=FI, n=n, index=index)
        for i in range(2, n):
            out[index['Ish' + str(i)]] = dIshk(pop=pop, j=i, FI=FI, FT=FT, index=index)
        if n > 1:
            out[index['Ish' + str(n)]] = dIshNI(pop=pop, FI=FI, FT=FT, n=n, index=index)

        #   in isolation
        out[index['I_i1']] = dI_i1(pop=pop, fi=f_phi_[2], FP=FP, FT=FT, FI=FI, n=n, index=index)
        for i in range(2, n + 1):
            out[index['I_i' + str(i)]] = dI_ik(pop=pop, j=i, FI=FI, FT=FT, index=index)

        # Recovered
        out[index['R__']] = dR(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, fdead_i=fdead_i_, FI=FI, FIp=FIp, n=n,
                               index=index, vac=vac_)

        # Dead, awaiting unsafe funeral
        out[index['F__']] = dF(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, FI=FI, FIp=FIp, FF=FF, n=n, index=index,
                               d_h=d_ph_[1], d_p=d_ph_[0])

        # Unsafely buried
        out[index['B_f']] = dB_f(pop=pop, FF=FF, index=index)

        # Safely buried
        out[index['B_j']] = dB_j(pop=pop, fdead_i=fdead_i_, FI=FI, FIp=FIp, n=n, index=index, d_h=d_ph_[1],
                                 d_p=d_ph_[0], fdead_h=fdead_h_, fdead_p=fdead_p_)

        # return population of all compartments
        return out

    #------------- Solve ODEs -------------------

    # Initial values
    pop0 = [0 for i in np.arange(len(index))]

    # All individuals are susceptible
    pop0[index['S__']] = N - P0
    # except the initial set of infected
    pop0[index['P__1']] = P0

    # solve ODE
    #   function f takes population as input and returns new one as output
    #   run in (max) daily steps for duration (days)
    soln = solve_ivp(f,
                     [0, days],
                     pop0,
                     t_eval=np.arange(0, days),
                     dense_output=True,
                     max_step=1)

    # save population at all compartments and time steps
    np.savetxt(pathResult + "/ebola_" + name + ".txt", soln.y)
    # save variables at all time steps
    np.savetxt(pathResult + "/ebolaVar_" + name + ".txt", rec, fmt='%.5f')
    return (name)
