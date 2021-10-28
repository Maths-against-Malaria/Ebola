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

from parameters_3 import *  # basic scenario (all parameters)
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
        m = m,
        R0=R0,
        cP=cP,
        cI=cI,
        cIh=cIh,
        cF=cF,
        P0=P0,  # P(0)
        t_iso=t_iso,
        I_iso = I_iso,
        f_p1=f_p1,
        f_h1=f_h1,
        f_p2=f_p2,
        f_h2=f_h2,
        k=k,
        d_p1=d_p1,
        d_h1=d_h1,
        d_p2=d_p2,
        d_h2=d_h2,
        l=l,
        f_tb=f_tb,
        ph=ph,
        qmax=qmax,
        cmax=cmax,
        fc = fc,
        t_vac=t_vac,
        Nvac=Nvac,
        nameIn='',
        pathOut='results'):
    name = str(days) + '_' \
           + str(NE) + '_' \
           + str(NP) + '_' \
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
           + str(m) + '_' \
           + str(R0) + '_' \
           + str(cP) + '_' \
           + str(cI) + '_' \
           + str(cIh) + '_' \
           + str(cF) + '_' \
           + str(P0) + '_' \
           + str(t_iso) + '_' \
           + str(I_iso) + '_' \
           + str(f_p1) + '_' \
           + str(f_h1) + '_' \
           + str(f_p2) + '_' \
           + str(f_h2) + '_' \
           + str(k) + '_' \
           + str(l) + '_' \
           + str(f_tb) + '_' \
           + str(ph) + '_' \
           + str(qmax) + '_' \
           + str(cmax) + '_' \
           + str(fc)  + '_' \
           + str(t_vac) + '_' \
           + str(Nvac) + '_' \
           + nameIn

    print("'" + name + "',")
    #     + str(d_p1) + '_' \
    #     + str(d_h1) + '_' \
    #     + str(d_p2) + '_' \
    #     + str(d_h2) + '_' \
    # compute values that do not change by time (or population)
    Nerls = [NE, NP, NIp, NIh, NIi]
    index = indexFunction(Nerls)
    #print(index)

    cD = R0 / (cP * DP + cI * DI + cF * DF)
    cD = float(cD)

    betaP = cP * cD
    betaIp = cI * cD
    betaIh = cIh * cD
    betaF = cF * cD

    FE = NE/DE # epsilon
    FP = NP/DP # phi
    FI = NIi/DI # gamma
    FF = 1/FI # phi2
    FT = 1/DI # alpha

    ####################################################
    #rec = [[-10000 for i in np.arange(2)] for j in np.arange(days + 1)]
    rec = [[-10000 for i in np.arange(10 + len(index))] for j in np.arange(days + 1)]

    def f(t, pop):
        # Initialize
        out = [0 for i in np.arange(len(index))]

        ## Values that change by time or population but are not differential equations
        # compute once per time and then use constants
        t_iso_ = t_iso
        if (t < t_iso):
            if popsum(pop=pop, compartment ='I', script1='_', script2='h', Nerls=Nerls, index=index) + \
                    popsum(pop=pop, compartment='I', script1='_', script2='i', Nerls=Nerls, index=index) >= I_iso:
                t_iso_ = t-1
        fdead_i_ = fdead_i[m]
        fdead_h_ = fdead_h[m]
        fdead_p_ = fdead_p[m]
        fc_ = fct(t=t, t_iso=t_iso_, fc=fc)
        f_phi_ = f_phi(t=t, k=k, t_iso=t_iso_, f_p1=f_p1, f_h1=f_h1, f_p2=f_p2, f_h2=f_h2)
        d_ph_ = d_ph(t=t, l=l, t_iso=t_iso_, d_p1=d_p1, d_h1=d_h1, d_p2=d_p2, d_h2=d_h2)
        q_ = q(pop=pop, t=t, t_iso=t_iso_, qmax=qmax, Nerls=Nerls, index=index)
        cc = c(pop=pop, t=t, t_iso=t_iso_, cmax=cmax, Nerls=Nerls, index=index, FT=FT, FP=FP, NP = NP, f_iso=f_phi_[2])
        c_ = cc[0]
        la__ = la(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph,
                  q=q_, c = c_, fc = fc_, Nerls=Nerls, index=index)
        # la__ = la_Aliou(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph, q=q_, Nerls=Nerls, index=index)
        # la__ = la_aliou2(pop=pop, fiso=f_phi_[2], f_tb=f_tb, betaP=betaP, betaIp=betaIp, betaIh=betaIh, betaF=betaF, ph=ph,q=q_, Nerls=Nerls, index=index)

        # la__ = la(f_phi_[2], f_tb, pop, betaP, betaIp, betaIh, betaF)
        la_ = la__[0]
        ls_ = la__[1]  # ls(f_phi_[2], pop, betaP, betaIp, betaIh)
        lt_ = lt(la=la_, ls=ls_)
        lt_ = lt(la=la_, ls=ls_)

        vac_ = vac(pop=pop, index=index, t=t,t_vac=t_vac, Nvac=Nvac)
        #rec[int(t)] = [t, q_]
        rec[int(t)] = [t] + [q_] + [fc_] + f_phi_ + [c_] + la__ + [vac_] + np.ndarray.tolist(pop)
        #rec[j] = [t, q_]
        #j = j + 1
        # trec.append(t)
        # rec.append([t] + f_phi_+ [q_]+ la__ + [lt_])
        # ints.append(int(t))

        ## Differential equations for compartments
        # Susceptible
        out[index['S__']] = dS(pop=pop, lt=lt_, N=N, index=index, vac=vac_)

        # Latent
        out[index['E__1']] = dE__1(pop=pop, la=la_, N=N, FE=FE, index=index)

        for i in range(2, NE + 1):
            out[index['E__' + str(i)]] = dE__k(pop=pop, j=i, FE=FE, index=index)

        out[index['Es_1']] = dEs_1(pop=pop, ls=ls_, N=N, FT=FT, FE=FE, index=index)

        for i in range(2, NE + 1):
            out[index['Es_' + str(i)]] = dEs_k(pop=pop, j=i, FE=FE, FT=FT, index=index)

        out[index['Et_1']] = dEt_1(pop=pop, FT=FT, FE=FE, index=index)

        for i in range(2, NE + 1):
            out[index['Et_' + str(i)]] = dEt_k(pop=pop, j=i, FE=FE, FT=FT, index=index)

        # Podromal
        out[index['P__1']] = dP__1(pop=pop, FE=FE, FP=FP, NE=NE, index=index)

        for i in range(2, NP + 1):
            out[index['P__' + str(i)]] = dP__k(pop=pop, j=i, FP=FP, index=index)

        out[index['Ps_1']] = dPs_1(pop=pop, FE=FE, FT=FT, FP=FP, NE=NE, index=index)

        for i in range(2, NP + 1):
            out[index['Ps_' + str(i)]] = dPs_k(pop=pop, j=i, FP=FP, FT=FT, index=index)

        out[index['Pt_1']] = dPt_1(pop=pop, FT=FT, FE=FE, FP=FP, NE=NE, index=index)

        for i in range(2, NP + 1):
            out[index['Pt_' + str(i)]] = dPt_k(pop=pop, j=i, FP=FP, FT=FT, index=index)

        # Fully Infectious
        # I_p home        
        out[index['I_p1']] = dI_p1(pop=pop, fp=f_phi_[0], FP=FP, FI=FI, NP=NP, index=index)

        for i in range(2, NIp + 1):
            out[index['I_p' + str(i)]] = dI_pk(pop=pop, j=i, FI=FI, index=index)

        # I_h hosp
        out[index['I_h1']] = dI_h1(pop=pop, fh=f_phi_[1], FP=FP, FI=FI, NP=NP, index=index)

        for i in range(2, NIh + 1):
            out[index['I_h' + str(i)]] = dI_hk(pop=pop, j=i, FI=FI, index=index)

        # Isp * home
        out[index['Isp1']] = dIsp1(pop=pop, fp=f_phi_[0], FP=FP, FT=FT, FI = FI, NP=NP, index=index)

        for i in range(2, NIp):
            out[index['Isp' + str(i)]] = dIspk(pop=pop, j=i, FI=FI, FT=FT, index=index)
        if NIp > 1:
            out[index['Isp' + str(NIp)]] = dIspNI(pop=pop, FI=FI, FT=FT, NIp=NIp, index=index)

        # Ish * hosp
        out[index['Ish1']] = dIsh1(pop=pop, fh=f_phi_[1], FP=FP, FT=FT, FI=FI, NP=NP, index=index)

        for i in range(2, NIh):
            out[index['Ish' + str(i)]] = dIshk(pop=pop, j=i, FI=FI, FT=FT, index=index)
        if NIh > 1:
            out[index['Ish' + str(NIh)]] = dIshNI(pop=pop, FI=FI, FT=FT, NIh=NIh, index=index)

        # I_i iso 
        out[index['I_i1']] = dI_i1(pop=pop, fi=f_phi_[2], FP=FP, FT=FT, FI=FI, NP=NP, index=index)

        for i in range(2, NIi + 1):
            out[index['I_i' + str(i)]] = dI_ik(pop=pop, j=i, FI=FI, FT=FT, index=index)

        # Recovered, Dead
        out[index['R__']] = dR(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, fdead_i=fdead_i_, FI=FI, NIp=NIp, NIh=NIh,
                               NIi=NIi, index=index, vac=vac_)

        out[index['F__']] = dF(pop=pop, fdead_p=fdead_p_, fdead_h=fdead_h_, FI=FI, FF=FF, NIp=NIp, NIh=NIh, index=index, d_h=d_ph_[1], d_p=d_ph_[0])

        out[index['B_f']] = dB_f(pop=pop, FF=FF, index=index)

        out[index['B_j']] = dB_j(pop=pop, fdead_i=fdead_i_, FI=FI, NIi=NIi, NIp=NIp, NIh = NIh, index=index, d_h=d_ph_[1], d_p=d_ph_[0], fdead_h=fdead_h_, fdead_p=fdead_p_)

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
    # np.savetxt(pathOut + "/ebolaInt_" + name + ".txt", trec)
    #    np.savetxt("20201010_R0_25/superinfection_parameters_" + name + ".txt", rec)
    return (name)
