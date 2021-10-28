# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 20:21:52 2020

@author: helle
"""

import matplotlib.pyplot as plt
import numpy as np
#import parameters_1

from index import *


#CD6700
colsA = ["#000000", "#801980", "#59B3E6", "#009980", "#E69900", "#CC6600", "#CD6700", "#0073B3"]
lab = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0']
#index = indexFunction(Nerls)
#pathIn = 'C:/Users/helle/PycharmProjects/ebola'
pathIn = 'results'
pathOut= 'plots'

def popsum2d(pops, Nerls):
    index = indexFunction(Nerls)
    indexComp = [1, 1 + Nerls[0], 1 + 2*Nerls[0], 1 + 3*Nerls[0],
                    1 + 3*Nerls[0]+Nerls[1], 1 + 3*Nerls[0]+2*Nerls[1], 1 + 3*Nerls[0]+3*Nerls[1],
                    1 + 3*Nerls[0]+3*Nerls[1]+Nerls[2], 1 + 3*Nerls[0]+3*Nerls[1]+Nerls[2]+Nerls[3], 1 + 3*Nerls[0]+3*Nerls[1]+Nerls[2]+Nerls[3]+Nerls[4],
                    1 + 3 * Nerls[0] + 3 * Nerls[1] + Nerls[2] + 2*Nerls[3] + Nerls[4], 1 + 3 * Nerls[0] + 3 * Nerls[1] + 2*Nerls[2] + 2*Nerls[3] + Nerls[4]]

    popSum = np.zeros((16,pops.shape[1]))
    popSum[0] = pops[index['S__']]
    popSum[1] = np.sum(pops[indexComp[0]:indexComp[1]], axis = 0)            # E__
    popSum[2] = np.sum(pops[indexComp[1]:indexComp[2]], axis = 0)       # Et_
    popSum[3] = np.sum(pops[indexComp[2]:indexComp[3]], axis = 0)     # Es_
    popSum[4] = np.sum(pops[indexComp[3]:indexComp[4]], axis = 0)     # P__
    popSum[5] = np.sum(pops[indexComp[4]:indexComp[5]], axis = 0)     # Pt_
    popSum[6] = np.sum(pops[indexComp[5]:indexComp[6]], axis = 0)     # Ps_
    popSum[7] = np.sum(pops[indexComp[6]:indexComp[7]], axis = 0)     # I_p
    popSum[8] = np.sum(pops[indexComp[7]:indexComp[8]], axis = 0)     # I_h
    popSum[9] = np.sum(pops[indexComp[8]:indexComp[9]], axis = 0)     # I_i
    popSum[10] = np.sum(pops[indexComp[9]:indexComp[10]], axis = 0)    # Ish
    popSum[11] = np.sum(pops[indexComp[10]:indexComp[11]], axis = 0)   # Isp
    popSum[12] = pops[index['F__']]
    popSum[13] = pops[index['B_j']]
    popSum[14] = pops[index['B_f']]
    popSum[15] = pops[index['R__']]
    return popSum


def getQ(pathIn, name, i, days):
    path = pathIn + '/ebolaVar_'+ name + '.txt'
    q_  = np.loadtxt(path)
    q__ = np.delete(q_, np.where(q_ == [-1.00000000e+04, -1.00000000e+04]), axis=0)
    q___ = np.interp(x=np.arange(days), xp=q__[:, 0], fp=q__[:, i])
    return q___

def plotEbolaParameters(names, savename, Nerls, days, pathIn=pathIn, pathOut = pathOut, n=10, col=colsA):
    index = indexFunction(Nerls)
    par_ji = getQ(pathIn=pathIn, name=names[1], i=1)
    par = [[-10000 for i in np.arange(2)] for j in np.arange(days)]
    for j in range(0, n):
        par[j] = np.empty(shape = [len(names), days])#s[1]])
        for i in range(0, len(names)):
            par[j][i] = getQ(pathIn=pathIn, name = names[i], i = j, days = days)

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
    plt.show()

def plotEbolaAll(names, savename, Nerls, pathIn=pathIn, pathOut = pathOut, col=colsA, lab=lab, q_max=False, tb = True, legendout=True, days=1000):
    if legendout:
        plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['axes.labelsize'] = 12
    index = indexFunction(Nerls)
    popSum0 = popsum2d(pops=np.loadtxt(pathIn + '/ebola_' + names[0] + '.txt'), Nerls=Nerls)
    s = np.shape(popSum0)
    s1 = int(max(s[1], days))
    print(s)
    if q_max == True:
        q_ = np.empty(shape = [len(names),s[1]])
        for i in range(0, len(names)):
            q_[i] = getQ(pathIn=pathIn, name = names[i], i = 1, days = days)

    #popSum = np.empty(shape = [len(names),s[0], s[1]])
    #popSum[0] = popSum0
    popSum = np.empty(shape=[len(names), s[0], s[1]])
    popSum[0] = popSum0
    for i in range(0,len(names)):
        pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
        popSum[i] = popsum2d(pops_i, Nerls=Nerls)


    fig = plt.figure()
    fig.set_size_inches(12, 12)

    # Susceptible
    p = fig.add_subplot(331)
    for i in range(0, len(names)):
        l4=l3=l2=l1=p.plot(popSum[i][0], label=lab[i], color=col[i])
    if legendout == False:
        #p.legend(bbox_to_anchor=(-0.2, 1.1, 3.5, -0.2), loc="lower left", mode="expand", ncol=len(lab))
        #p.legend(bbox_to_anchor=(-0.2, 1.1, 2.3, -0.2), loc="lower left", mode="expand", ncol=3)
        p.legend()
    p.set_ylabel('Susceptible ind.')
    p.set_title(label='A.', loc='left')
    p.set_ylim(bottom=0)

    # Latent
    p = fig.add_subplot(332)
    p.plot(popSum[0][1], color=col[0], linestyle='-')
    if tb == True:
        p.plot(popSum[0][1], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][3], color=col[0], linestyle='--', label='not yet traced back')
        p.plot(popSum[0][2], color=col[0], linestyle=':', label='traced back')
    if legendout == False:
        p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][1], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][2], color=col[i], linestyle=':')
            p.plot(popSum[i][3], color=col[i], linestyle='--')
    p.set_ylabel('Latent ind.')
    p.set_title(label='B.', loc='left')
    p.set_ylim(bottom=0)

    # Podromal
    p = fig.add_subplot(333)
    if tb == True:
        if legendout:
            p.plot(popSum[0][4], color=col[0], linestyle='-', label='never traced back (and unsafe funeral)')
        else:
            p.plot(popSum[0][4], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][6], color=col[0], linestyle='--', label='not yet traced back')
        if q_max == True:
            if legendout:
                l3=p.plot(np.multiply(popSum[0][5], q_[0]), color=col[0], linestyle=':', label='traced back, in ward (and safe funeral)')
            else:
                p.plot(np.multiply(popSum[0][5], q_[0]), color=col[0], linestyle=':', label='traced back, in ward')
            p.plot(np.multiply(popSum[0][5], np.multiply(q_[0],-1)+1), color=col[0], linestyle='-.', label='traced back, not in ward')
            for i in range(0, len(names)):
                p.plot(np.multiply(popSum[i][5], q_[i]), color=col[i], linestyle=':')
                p.plot(np.multiply(popSum[i][5], np.multiply(q_[i],-1)+1), color=col[i], linestyle='-.')
        if q_max == False:
            if legendout:
                p.plot(popSum[0][5], color=col[0], linestyle=':', label='traced back (and safe funeral)')
            else:
                p.plot(popSum[0][5], color=col[0], linestyle=':', label='traced back')
    if tb == False:
        if legendout:
            p.plot(popSum[0][4], color=col[0], linestyle='-', label='unsafe funeral')
            p.plot(popSum[0][4], color=col[0], linestyle=':', label='safe funeral')
            p.plot(popSum[0][4], color='w', linestyle=':', label=' ')
    if legendout:
        for i in range(0, len(names)):
            p.plot(popSum[i][4], label=lab[i], color=col[i], linestyle='-')
        #p.legend((l2,l3,l4),bbox_to_anchor=(-2.6, 1.3, 1.5, 0), loc="lower left", mode="expand", ncol=3)
        #p.legend(bbox_to_anchor=(-0.4, 1.1, 1.5, 0), loc="lower left", mode="expand", ncol=1)
        p.legend(bbox_to_anchor=(-3, -3, 3.5, 0.5), loc="upper left", mode="expand", ncol=3)
                 #handles=reversed(p.legend().legendHandles), labels=lab + ['a', 'b', 'c'])
        #print('a')
    else:
        p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][4], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][5], color=col[i], linestyle=':')
            p.plot(popSum[i][6], color=col[i], linestyle='--')
    p.set_ylabel('Podromal ind.')
    p.set_title(label='C.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected at home
    p = fig.add_subplot(334)
    if tb == True:
        p.plot(popSum[0][7], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][11], color=col[0], linestyle='--', label='not yet traced back')
        if legendout == False:
            p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][7], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][11], color=col[i], linestyle='--')
    p.set_ylabel('Fully inf. ind. at home')
    p.set_title(label='D.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in hospital
    p = fig.add_subplot(335)
    if tb == True:
        p.plot(popSum[0][8], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][10], color=col[0], linestyle='--', label='not yet traced back')
        if legendout == False:
            p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][8], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][10], color=col[i], linestyle=':')
    p.set_ylabel('Fully inf. ind. in hospital')
    p.set_title(label='E.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in isolation
    p = fig.add_subplot(336)
    if q_max == True:
        p.plot(np.multiply(popSum[0][9], q_[0]), color=col[0], linestyle='-', label='in ward')
        p.plot(np.multiply(popSum[0][9], np.multiply(q_[0], -1) + 1), color=col[0], linestyle='-.',
                 label='not in ward')
        if legendout == False:
            p.legend()
        for i in range(0, len(names)):
            p.plot(np.multiply(popSum[i][9], q_[i]), color=col[i], linestyle='-')
            p.plot(np.multiply(popSum[i][9], np.multiply(q_[i], -1) + 1), color=col[i], linestyle='-.')
    if q_max == False:
        #p.plot(popSum[0][9], color=col[0], linestyle='-', label = 'isolated in ward')
        #p.legend()
        for i in range(0, len(names)):
            p.plot(popSum[i][9], label=lab[i], color=col[i], linestyle='-')
    #p.ylim([-5, np.ndarray.max(popSum[:, 9]) * 1.05])
    p.set_ylabel('Fully inf. ind. in isolation')
    p.set_title(label='F.', loc='left')
    p.set_ylim(bottom=0)

    # Unsafe funerals
    p = fig.add_subplot(337)
    for i in range(0, len(names)):
        p.plot(popSum[i][12], label=lab[i], color=col[i], linestyle='-')
    p.set_ylabel('Unsafe funerals')
    p.set_title(label='G.', loc='left')
    p.set_ylim(bottom=0)

    # Buried
    p = fig.add_subplot(338)
    p.plot(popSum[0][13], color=col[0], linestyle=':', label='safely')
    p.plot(popSum[0][14], color=col[0], linestyle='-', label='unsafely')
    #p.legend()
    if legendout == False:
        p.legend()#bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand", ncol=2)
    for i in range(0, len(names)):
        p.plot(popSum[i][13], color=col[i], linestyle=':')
        p.plot(popSum[i][14], label=lab[i], color=col[i], linestyle='-')
    p.set_ylabel('Buried ind.')
    p.set_title(label='H.', loc='left')
    p.set_ylim(bottom=0)

    # Recovered
    p = fig.add_subplot(339)
    for i in range(0, len(names)):
        p.plot(popSum[i][15], label=lab[i], color=col[i], linestyle='-')
    p.set_ylabel('Recovered ind')
    p.set_title(label='I.', loc='left')
    p.set_ylim(bottom=0)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
    #plt.tight_layout()
    #ig = plt.figure(figsize=(15, 13))
    #fig.legend(loc="upper center", bbox_to_anchor=(0.5,1))#, mode="expand", ncol=3)
    #plt.legend(bbox_to_anchor=(0.5, 1), loc="upper center")#, mode="expand", ncol=3),
    #p.legend(bbox_to_anchor=(-2.5, 1.1, 3.5, 0.5), loc="upper left", mode="expand", ncol=3),

    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.savefig(pathOut + '/Ebola_all_' + savename + '.pdf', dpi=100)
    plt.show()


    plt.plot(popSum[0][0] + popSum[0][15], label='healthy (S+R)', color=col[0], linestyle='-')
    plt.plot(np.sum(popSum[0][1:12], axis=0), label='infected (E+P+I)', color=col[0], linestyle='--')
    plt.plot(np.sum(popSum[0][12:15], axis=0), label='dead (D+B)', color=col[0], linestyle=':')
    for i in range(0, len(names)):
        plt.plot(popSum[i][0] + popSum[i][15], label=lab[i],color=col[i], linestyle='-')
        print('healthy ' + lab[i] + ': '+ str(round((popSum[i][0] + popSum[i][15])[-1])))
        plt.plot(np.sum(popSum[i][1:12], axis=0), color=col[i], linestyle='--')
        print('infected ' + lab[i] + ': ' + str(round(np.sum(popSum[i][1:12], axis=0)[-1])))
        plt.plot(np.sum(popSum[i][12:15], axis=0), color=col[i], linestyle=':')
        print('dead ' + lab[i] + ': ' + str(round(np.sum(popSum[i][12:15], axis=0)[-1])))
        print('---------')
        plt.ylabel('Individuals')
    plt.legend()
    #plt.title(savename)
    plt.savefig(pathOut + '/Ebola_mix_' + savename + '.pdf', dpi=100)
    plt.show()
    return(savename)




def plotEbolaAll(names, savename, Nerls, pathIn=pathIn, pathOut = pathOut, col=colsA, lab=lab, q_max=False, tb = True, legendout=True, days=1000):
    if legendout:
        plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['axes.labelsize'] = 12
    index = indexFunction(Nerls)
    popSum0 = popsum2d(pops=np.loadtxt(pathIn + '/ebola_' + names[0] + '.txt'), Nerls=Nerls)
    s = np.shape(popSum0)
    s1 = int(max(s[1], days))
    print(s)
    if q_max == True:
        q_ = np.empty(shape = [len(names),s[1]])
        for i in range(0, len(names)):
            q_[i] = getQ(pathIn=pathIn, name = names[i], i = 1, days = days)

    #popSum = np.empty(shape = [len(names),s[0], s[1]])
    #popSum[0] = popSum0
    popSum = np.empty(shape=[len(names), s[0], s[1]])
    popSum[0] = popSum0
    for i in range(0,len(names)):
        pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
        popSum[i] = popsum2d(pops_i, Nerls=Nerls)


    fig = plt.figure()
    fig.set_size_inches(12, 12)

    # Susceptible
    p = fig.add_subplot(331)
    for i in range(0, len(names)):
        l4=l3=l2=l1=p.plot(popSum[i][0], label=lab[i], color=col[i])
    if legendout == False:
        #p.legend(bbox_to_anchor=(-0.2, 1.1, 3.5, -0.2), loc="lower left", mode="expand", ncol=len(lab))
        #p.legend(bbox_to_anchor=(-0.2, 1.1, 2.3, -0.2), loc="lower left", mode="expand", ncol=3)
        p.legend()
    p.set_ylabel('Susceptible ind.')
    p.set_title(label='A.', loc='left')
    p.set_ylim(bottom=0)

    # Latent
    p = fig.add_subplot(332)
    p.plot(popSum[0][1], color=col[0], linestyle='-')
    if tb == True:
        p.plot(popSum[0][1], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][3], color=col[0], linestyle='--', label='not yet traced back')
        p.plot(popSum[0][2], color=col[0], linestyle=':', label='traced back')
    if legendout == False:
        p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][1], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][2], color=col[i], linestyle=':')
            p.plot(popSum[i][3], color=col[i], linestyle='--')
    p.set_ylabel('Latent ind.')
    p.set_title(label='B.', loc='left')
    p.set_ylim(bottom=0)

    # Podromal
    p = fig.add_subplot(333)
    if tb == True:
        if legendout:
            p.plot(popSum[0][4], color=col[0], linestyle='-', label='never traced back (and unsafe funeral)')
        else:
            p.plot(popSum[0][4], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][6], color=col[0], linestyle='--', label='not yet traced back')
        if q_max == True:
            if legendout:
                l3=p.plot(np.multiply(popSum[0][5], q_[0]), color=col[0], linestyle=':', label='traced back, in ward (and safe funeral)')
            else:
                p.plot(np.multiply(popSum[0][5], q_[0]), color=col[0], linestyle=':', label='traced back, in ward')
            p.plot(np.multiply(popSum[0][5], np.multiply(q_[0],-1)+1), color=col[0], linestyle='-.', label='traced back, not in ward')
            for i in range(0, len(names)):
                p.plot(np.multiply(popSum[i][5], q_[i]), color=col[i], linestyle=':')
                p.plot(np.multiply(popSum[i][5], np.multiply(q_[i],-1)+1), color=col[i], linestyle='-.')
        if q_max == False:
            if legendout:
                p.plot(popSum[0][5], color=col[0], linestyle=':', label='traced back (and safe funeral)')
            else:
                p.plot(popSum[0][5], color=col[0], linestyle=':', label='traced back')
    if tb == False:
        if legendout:
            p.plot(popSum[0][4], color=col[0], linestyle='-', label='unsafe funeral')
            p.plot(popSum[0][4], color=col[0], linestyle=':', label='safe funeral')
            p.plot(popSum[0][4], color='w', linestyle=':', label=' ')
    if legendout:
        for i in range(0, len(names)):
            p.plot(popSum[i][4], label=lab[i], color=col[i], linestyle='-')
        #p.legend((l2,l3,l4),bbox_to_anchor=(-2.6, 1.3, 1.5, 0), loc="lower left", mode="expand", ncol=3)
        #p.legend(bbox_to_anchor=(-0.4, 1.1, 1.5, 0), loc="lower left", mode="expand", ncol=1)
        p.legend(bbox_to_anchor=(-3, -3, 3.5, 0.5), loc="upper left", mode="expand", ncol=3)
                 #handles=reversed(p.legend().legendHandles), labels=lab + ['a', 'b', 'c'])
        #print('a')
    else:
        p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][4], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][5], color=col[i], linestyle=':')
            p.plot(popSum[i][6], color=col[i], linestyle='--')
    p.set_ylabel('Podromal ind.')
    p.set_title(label='C.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected at home
    p = fig.add_subplot(334)
    if tb == True:
        p.plot(popSum[0][7], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][11], color=col[0], linestyle='--', label='not yet traced back')
        if legendout == False:
            p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][7], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][11], color=col[i], linestyle='--')
    p.set_ylabel('Fully inf. ind. at home')
    p.set_title(label='D.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in hospital
    p = fig.add_subplot(335)
    if tb == True:
        p.plot(popSum[0][8], color=col[0], linestyle='-', label='never traced back')
        p.plot(popSum[0][10], color=col[0], linestyle='--', label='not yet traced back')
        if legendout == False:
            p.legend()
    for i in range(0, len(names)):
        p.plot(popSum[i][8], label=lab[i], color=col[i], linestyle='-')
        if tb == True:
            p.plot(popSum[i][10], color=col[i], linestyle=':')
    p.set_ylabel('Fully inf. ind. in hospital')
    p.set_title(label='E.', loc='left')
    p.set_ylim(bottom=0)

    # Fully infected in isolation
    p = fig.add_subplot(336)
    if q_max == True:
        p.plot(np.multiply(popSum[0][9], q_[0]), color=col[0], linestyle='-', label='in ward')
        p.plot(np.multiply(popSum[0][9], np.multiply(q_[0], -1) + 1), color=col[0], linestyle='-.',
                 label='not in ward')
        if legendout == False:
            p.legend()
        for i in range(0, len(names)):
            p.plot(np.multiply(popSum[i][9], q_[i]), color=col[i], linestyle='-')
            p.plot(np.multiply(popSum[i][9], np.multiply(q_[i], -1) + 1), color=col[i], linestyle='-.')
    if q_max == False:
        #p.plot(popSum[0][9], color=col[0], linestyle='-', label = 'isolated in ward')
        #p.legend()
        for i in range(0, len(names)):
            p.plot(popSum[i][9], label=lab[i], color=col[i], linestyle='-')
    #p.ylim([-5, np.ndarray.max(popSum[:, 9]) * 1.05])
    p.set_ylabel('Fully inf. ind. in isolation')
    p.set_title(label='F.', loc='left')
    p.set_ylim(bottom=0)

    # Unsafe funerals
    p = fig.add_subplot(337)
    for i in range(0, len(names)):
        p.plot(popSum[i][12], label=lab[i], color=col[i], linestyle='-')
    p.set_ylabel('Unsafe funerals')
    p.set_title(label='G.', loc='left')
    p.set_ylim(bottom=0)

    # Buried
    p = fig.add_subplot(338)
    p.plot(popSum[0][13], color=col[0], linestyle=':', label='safely')
    p.plot(popSum[0][14], color=col[0], linestyle='-', label='unsafely')
    #p.legend()
    if legendout == False:
        p.legend()#bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand", ncol=2)
    for i in range(0, len(names)):
        p.plot(popSum[i][13], color=col[i], linestyle=':')
        p.plot(popSum[i][14], label=lab[i], color=col[i], linestyle='-')
    p.set_ylabel('Buried ind.')
    p.set_title(label='H.', loc='left')
    p.set_ylim(bottom=0)

    # Recovered
    p = fig.add_subplot(339)
    for i in range(0, len(names)):
        p.plot(popSum[i][15], label=lab[i], color=col[i], linestyle='-')
    p.set_ylabel('Recovered ind')
    p.set_title(label='I.', loc='left')
    p.set_ylim(bottom=0)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
    #plt.tight_layout()
    #ig = plt.figure(figsize=(15, 13))
    #fig.legend(loc="upper center", bbox_to_anchor=(0.5,1))#, mode="expand", ncol=3)
    #plt.legend(bbox_to_anchor=(0.5, 1), loc="upper center")#, mode="expand", ncol=3),
    #p.legend(bbox_to_anchor=(-2.5, 1.1, 3.5, 0.5), loc="upper left", mode="expand", ncol=3),

    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.savefig(pathOut + '/Ebola_all_' + savename + '.pdf', dpi=100)
    plt.show()


    plt.plot(popSum[0][0] + popSum[0][15], label='healthy (S+R)', color=col[0], linestyle='-')
    plt.plot(np.sum(popSum[0][1:12], axis=0), label='infected (E+P+I)', color=col[0], linestyle='--')
    plt.plot(np.sum(popSum[0][12:15], axis=0), label='dead (D+B)', color=col[0], linestyle=':')
    for i in range(0, len(names)):
        plt.plot(popSum[i][0] + popSum[i][15], label=lab[i],color=col[i], linestyle='-')
        print('healthy ' + lab[i] + ': '+ str(round((popSum[i][0] + popSum[i][15])[-1])))
        plt.plot(np.sum(popSum[i][1:12], axis=0), color=col[i], linestyle='--')
        print('infected ' + lab[i] + ': ' + str(round(np.sum(popSum[i][1:12], axis=0)[-1])))
        plt.plot(np.sum(popSum[i][12:15], axis=0), color=col[i], linestyle=':')
        print('dead ' + lab[i] + ': ' + str(round(np.sum(popSum[i][12:15], axis=0)[-1])))
        print('---------')
        plt.ylabel('Individuals')
    plt.legend()
    #plt.title(savename)
    plt.savefig(pathOut + '/Ebola_mix_' + savename + '.pdf', dpi=100)
    plt.show()
    return(savename)


def plotEbolaScenarios(namesAllScenarios, savename, Nerls, pathIn=pathIn, pathOut=pathOut, col=colsA, lab=lab):
    fig = plt.figure()
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['axes.labelsize'] = 12

    # M sets of scenarios
    M = len(namesAllScenarios)
    #fig.set_size_inches(12, M * 4)
    fig.set_size_inches(12, 12)

    alphabet = []
    for letter in range(65, 91):
        alphabet.append(chr(letter))

    for m in range(0, M):
        mm = 3 + (m - 1) * 3

        names = namesAllScenarios[m]
        popSum0 = popsum2d(pops=np.loadtxt(pathIn + '/ebola_' + names[0] + '.txt'), Nerls=Nerls)
        s = np.shape(popSum0)
        popSum = np.empty(shape=[len(names), s[0], s[1]])
        popSum[0] = popSum0

        for i in range(0, len(names)):
            pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
            popSum[i] = popsum2d(pops_i, Nerls=Nerls)

        # Susceptible
        p = fig.add_subplot(M, 3, mm + 1)

        for i in range(0, len(names)):
            p.plot(popSum[i][0], label= r'$f_{iso} = $' + lab[i], color=col[i], linestyle='-')
        if m == 0:
            #p.legend(bbox_to_anchor=(-0.05, -4.9, 3.9, 1), loc="upper center", mode="expand", ncol=5)
            p.legend(bbox_to_anchor=(-0.05, -4.9, 3.9, 1), loc="upper center", mode="expand", ncol=5)

        p.set_ylabel('Susceptible ind.')
        p.set_title(label=alphabet[mm] + '.', loc='left')
        p.set_ylim(bottom=0)

        # Infected
        p = fig.add_subplot(M, 3, mm + 2)
        for i in range(0, len(names)):
            p.plot(np.sum(popSum[i][1:12], axis=0), label=lab[i], color=col[i], linestyle='-')
        p.set_ylabel('Infected ind.')
        p.set_title(label=alphabet[mm+1] + '.', loc='left')
        p.set_ylim(bottom=0)

        # Dead
        p = fig.add_subplot(M, 3, mm + 3)
        for i in range(0, len(names)):
            p.plot(np.sum(popSum[i][12:15], axis=0), label=lab[i], color=col[i], linestyle='-')
        p.set_ylabel('Dead ind.')
        p.set_title(label=alphabet[mm+2] + '.', loc='left')
        p.set_ylim(bottom=0)


    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.25)

    plt.savefig(pathOut + '/Ebola_scenarios_' + savename + '.pdf', dpi=100,bbox_inches='tight')
    plt.show()
    return (savename)
