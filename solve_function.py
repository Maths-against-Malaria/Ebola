# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 15:03:23 2020

@author: helle
"""

# This is the basic model as a function
# It imports parameters (l 15) to organise the simulation of different scenarios (each to be saved with own name).
# It saves the result to a textfile (last lines) [and some more things], that can be used for plotting etc.

import numpy as np
from scipy.integrate import solve_ivp

from parameters_1 import *  # basic scenario (all parameters)
from index import *
from functions import *
from differential_equations import *


def modelEbola(
        days=days,
        NE=NE,
        NP=NP,
        NIp=NIp,
        NIh=NIh,
        NIi=NIi,
        N=N,
        DE=DE,
        DP=DP,
        DI=DI,
        DF=DF,
        DT=DT,
        fdead_p=fdead_p,
        fdead_h=fdead_h,
        fdead_i=fdead_i,
        R0=R0,
        cP=cP,
        cIp=cIp,
        cIh=cIh,
        cF=cF,
        P0=P0,  # P(0)
        t_iso=t_iso,
        f_p1=f_p1,
        f_h1=f_h1,
        f_p2=f_p2,
        f_h2=f_h2,
        k=k,
        f_tb=f_tb,
        ph=ph,
        qmax=qmax,
        f_di_h=f_di_h,
        f_di_p=f_di_p,
        t_vac=t_vac,
        Nvac=Nvac,
        nameIn='',
        pathOut='results'):
    name = str(days) + '_' \
           + str(NE) + '_' \
           + str(NP) + '_' \
           + str(NIp) + '_' \
           + str(NIp) + '_' \
           + str(NIh) + '_' \
           + str(NIi) + '_' \
           + str(N) + '_' \
           + str(DE) + '_' \
           + str(DP) + '_' \
           + str(DI) + '_' \
           + str(DF) + '_' \
           + str(DT) + '_' \
           + str(fdead_p) + '_' \
           + str(fdead_h) + '_' \
           + str(fdead_i) + '_' \
           + str(R0) + '_' \
           + str(cP) + '_' \
           + str(cIp) + '_' \
           + str(cIh) + '_' \
           + str(cF) + '_' \
           + str(P0) + '_' \
           + str(t_iso) + '_' \
           + str(f_p1) + '_' \
           + str(f_h1) + '_' \
           + str(f_p2) + '_' \
           + str(f_h2) + '_' \
           + str(k) + '_' \
           + str(f_tb) + '_' \
           + str(ph) + '_' \
           + str(qmax) + '_' \
           + str(f_di_h) + '_' \
           + str(f_di_p) + '_' \
           + str(t_vac) + '_' \
           + str(Nvac) + '_' \
           + nameIn

    print(name)

    # compute values that do not change by time (or population)
    Nerls = [NE, NP, NIp, NIh, NIi]
    index = indexFunction(Nerls)

    cD = R0 / (cP * DP + (cIh + cIp) * 0.5 * DI + cF * DF)
    cD = float(cD)

    betaP = cP * cD
    betaIp = cIp * cD
    betaIh = cIh * cD
    betaF = cF * cD

    ####################################################
    rec = [[-10000 for i in np.arange(2)] for j in np.arange(days + 1)]

    def f(t, pop):
        # Initialize
        out = [0 for i in np.arange(len(index))]

        ## Values that change by time or population but are not differential equations
        # compute once per time and then use constants
        f_phi_ = f_phi(t=t, k=k, t_iso=t_iso, f_p1=f_p1, f_h1=f_h1, f_p2=f_p2, f_h2=f_h2)
        q_ = q(pop=pop, t=t, t_iso=t_iso, qmax=qmax, Nerls=Nerls, index=index)
        la__ = la(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph,
                  q=q_, Nerls=Nerls, index=index)
        # la__ = la_Aliou(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph, q=q_, Nerls=Nerls, index=index)
        # la__ = la_aliou2(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph,q=q_, Nerls=Nerls, index=index)

        # la__ = la(f_phi_[2], f_tb, pop, betaP, betaIp, betaIh, betaF)
        la_ = la__[0]
        ls_ = la__[1]  # ls(f_phi_[2], pop, betaP, betaIp, betaIh)
        lt_ = lt(la=la_, ls=ls_)

        vac_ = vac(pop=pop, index=index, t=t,t_vac=t_vac, Nvac=Nvac)

        rec[int(t)] = [t, q_]
        #rec[j] = [t, q_]
        #j = j + 1
        # trec.append(t)
        # rec.append([t] + f_phi_+ [q_]+ la__ + [lt_])
        # ints.append(int(t))

        ## Differential equations for compartments
        # Susceptible
        out[index['S__']] = dS(pop=pop, lt=lt_, N=N, index=index, vac=vac_)

        # Latent
        out[index['E__1']] = dE__1(pop=pop, la=la_, N=N, DE=DE, index=index)

        for i in range(2, NE + 1):
            out[index['E__' + str(i)]] = dE__k(pop=pop, j=i, DE=DE, index=index)

        out[index['Es_1']] = dEs_1(pop=pop, ls=ls_, N=N, DT=DT, DE=DE, index=index)

        for i in range(2, NE + 1):
            out[index['Es_' + str(i)]] = dEs_k(pop=pop, j=i, DE=DE, DT=DT, index=index)

        out[index['Et_1']] = dEt_1(pop=pop, DT=DT, DE=DE, index=index)

        for i in range(2, NE + 1):
            out[index['Et_' + str(i)]] = dEt_k(pop=pop, j=i, DE=DE, DT=DT, index=index)

        # Podromal
        out[index['P__1']] = dP__1(pop=pop, DE=DE, DP=DP, NE=NE, index=index)

        for i in range(2, NP + 1):
            out[index['P__' + str(i)]] = dP__k(pop=pop, j=i, DP=DP, index=index)

        out[index['Ps_1']] = dPs_1(pop=pop, DE=DE, DT=DT, DP=DP, NE=NE, index=index)

        for i in range(2, NP + 1):
            out[index['Ps_' + str(i)]] = dPs_k(pop=pop, j=i, DP=DP, DT=DT, index=index)

        out[index['Pt_1']] = dPt_1(pop=pop, DT=DT, DE=DE, DP=DP, NE=NE, index=index)

        for i in range(2, NP + 1):
            out[index['Pt_' + str(i)]] = dPt_k(pop=pop, j=i, DP=DP, DT=DT, index=index)

        # Fully Infectious
        # I_p home        
        out[index['I_p1']] = dI_p1(pop=pop, fp=f_phi_[0], DP=DP, DI=DI, NP=NP, index=index)

        for i in range(2, NIp + 1):
            out[index['I_p' + str(i)]] = dI_pk(pop=pop, j=i, DI=DI, index=index)

        # I_h hosp
        out[index['I_h1']] = dI_h1(pop=pop, fh=f_phi_[1], DP=DP, DI=DI, NP=NP, index=index)

        for i in range(2, NIh + 1):
            out[index['I_h' + str(i)]] = dI_hk(pop=pop, j=i, DI=DI, index=index)

        # Isp * home
        out[index['Isp1']] = dIsp1(pop=pop, fp=f_phi_[0], DP=DP, DT=DT, NP=NP, index=index)

        for i in range(2, NIp):
            out[index['Isp' + str(i)]] = dIspk(pop=pop, j=i, DI=DI, DT=DT, index=index)
        if NIp > 1:
            out[index['Isp' + str(NIp)]] = dIspNI(pop=pop, DI=DI, DT=DT, NIp=NIp, index=index)

        # Ish * hosp
        out[index['Ish1']] = dIsh1(pop=pop, fh=f_phi_[1], DP=DP, DT=DT, DI=DI, NP=NP, index=index)

        for i in range(2, NIh):
            out[index['Ish' + str(i)]] = dIshk(pop=pop, j=i, DI=DI, DT=DT, index=index)
        if NIh > 1:
            out[index['Ish' + str(NIh)]] = dIshNI(pop=pop, DI=DI, DT=DT, NIh=NIh, index=index)

        # I_i iso 
        out[index['I_i1']] = dI_i1(pop=pop, fi=f_phi_[2], DP=DP, DT=DT, DI=DI, NP=NP, index=index)

        for i in range(2, NIi + 1):
            out[index['I_i' + str(i)]] = dI_ik(pop=pop, j=i, DI=DI, DT=DT, index=index)

        # Recovered, Dead
        out[index['R__']] = dR(pop=pop, fdead_p=fdead_p, fdead_h=fdead_h, fdead_i=fdead_i, DI=DI, NIp=NIp, NIh=NIh,
                               NIi=NIi, index=index, vac=vac_)

        out[index['F__']] = dF(pop=pop, fdead_p=fdead_p, fdead_h=fdead_h, DI=DI, DF=DF, NIp=NIp, NIh=NIh, index=index, f_di_h=f_di_h, f_di_p=f_di_p)

        out[index['B_f']] = dB_f(pop=pop, DF=DF, index=index)

        out[index['B_j']] = dB_j(pop=pop, fdead_i=fdead_i, DI=DI, NIi=NIi, index=index, f_di_h=f_di_h, f_di_p=f_di_p)

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
    np.savetxt(pathOut + "/ebolaVar_" + name + ".txt", rec)
    # np.savetxt(pathOut + "/ebolaInt_" + name + ".txt", trec)
    #    np.savetxt("20201010_R0_25/superinfection_parameters_" + name + ".txt", rec)
    return (name)
