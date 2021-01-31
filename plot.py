# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 20:21:52 2020

@author: helle
"""

import matplotlib.pyplot as plt
import numpy as np


from index import *

Nerls=[NE, NP, NIp, NIh, NIi]

def popsum2d(pops, Nerls=Nerls):
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

colsA = ["#000000", "#59B3E6", "#009980", "#E69900", "#CC6600", "#0073B3", 'red']
lab = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0']
index = indexFunction(Nerls)
path = 'C:/Users/helle/PycharmProjects/ebola'

def plotEbola(names, savename, path=path, col=colsA, lab=lab, type='all', Nerls=Nerls, legendtitle=''):
    index = indexFunction(Nerls)
    popSum0 = popsum2d(pops=np.loadtxt(path +'/ebola_' + names[0] + '.txt'), Nerls=Nerls)
    s = np.shape(popSum0)
    popSum = np.empty(shape = [len(names),s[0], s[1]])
    popSum[0] = popSum0
    for i in range(1,len(names)):
        pops_i = np.loadtxt(path +'/ebola_' + names[i] + '.txt')
        #plt.plot(np.sum(pops_i, axis=0), color=col[i])
        popSum[i] = popsum2d(pops_i, Nerls=Nerls)
    #plt.show()
    # total population (for test)
    #for i in range(0, len(names)):
    #    plt.plot(np.sum(popSum[i], axis =0), color=col[i])
    #plt.show()
    # Susceptible
    if type == 'all' or type == 'S':
        for i in range(0,len(names)):
            plt.plot(popSum[i][0], label = lab[i], color=col[i])
        plt.legend(title=legendtitle)
        plt.ylabel('Susceptible ind.')
        plt.savefig (path + '/Ebola_S_' + savename + '.pdf', dpi=100)
        plt.show()

    # Latent
    if type == 'all' or type == 'E':
        plt.plot(popSum[0][1], color=col[0], linestyle='-')
        #plt.plot(popSum[0][1], color=col[0], linestyle='-', label='never traced back')
        #plt.plot(popSum[0][3], color=col[0], linestyle='--', label='not yet traced back')
        #plt.plot(popSum[0][2], color=col[0], linestyle=':', label='traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][1], label=lab[i], color=col[i], linestyle='-')
            #plt.plot(popSum[i][2], color=col[i], linestyle=':')
            #plt.plot(popSum[i][3], color=col[i], linestyle='--')
        plt.legend(title=legendtitle)
        plt.ylabel('Latent ind.')
        plt.savefig(path + '/Ebola_E_' + savename + '.pdf', dpi=100)
        plt.show()
    # Podromal
    if type == 'all' or type == 'P':
        plt.plot(popSum[0][4], color=col[0], linestyle='-')
        #plt.plot(popSum[0][4], color=col[0], linestyle='-', label='never traced back')
        #plt.plot(popSum[0][6], color=col[0], linestyle='--', label='not yet traced back')
        #plt.plot(popSum[0][5], color=col[0], linestyle=':', label='traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][4], label=lab[i], color=col[i], linestyle='-')
            #plt.plot(popSum[i][5], color=col[i], linestyle=':')
            #plt.plot(popSum[i][6], color=col[i], linestyle='--')
        plt.legend(title=legendtitle)
        plt.ylabel('Podromal ind.')
        plt.savefig(path + '/Ebola_P_' + savename + '.pdf', dpi=100)
        plt.show()
    # Fully infected at home
    if type == 'all' or type == 'Ip':
        plt.plot(popSum[0][7], color=col[0], linestyle='-')
        #plt.plot(popSum[0][7], color=col[0], linestyle='-', label='never traced back')
        #plt.plot(popSum[0][11], color=col[0], linestyle='--', label='not yet traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][7], label=lab[i], color=col[i], linestyle='-')
            #plt.plot(popSum[i][11], color=col[i], linestyle=':')
        plt.legend(title=legendtitle)
        plt.ylabel('Fully inf. ind. at home')
        plt.savefig(path + '/Ebola_Ip_' + savename + '.pdf', dpi=100)
        plt.show()
    # Fully infected in hospital
    if type == 'all' or type == 'Ih':
        plt.plot(popSum[0][8], color=col[0], linestyle='-')
        #plt.plot(popSum[0][8], color=col[0], linestyle='-', label='never traced back')
        #plt.plot(popSum[0][10], color=col[0], linestyle='--', label='not yet traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][8], label=lab[i], color=col[i], linestyle='-')
            #plt.plot(popSum[i][10], color=col[i], linestyle=':')
        plt.legend(title=legendtitle)
        plt.ylabel('Fully inf. ind. in hospital')
        plt.savefig(path + '/Ebola_Ih_' + savename + '.pdf', dpi=100)
        plt.show()
    # Fully infected in isolation
    if type == 'all' or type == 'Ii':
        for i in range(0, len(names)):
            #plt.plot(popSum[i][9], label=lab[i], color=col[i], linestyle=':')
            plt.plot(popSum[i][9], label=lab[i], color=col[i], linestyle='-')
        plt.legend(title=legendtitle)
        plt.ylabel('Fully inf. ind. in isolation')
        plt.savefig(path + '/Ebola_Ii_' + savename + '.pdf', dpi=100)
        plt.show()
    # Buried
    if type == 'all' or type == 'B':
        plt.plot(popSum[0][14], color=col[0], linestyle='-', label='unsafely')
        plt.plot(popSum[0][13], color=col[0], linestyle='--', label='safely')
        for i in range(0, len(names)):
            plt.plot(popSum[i][13], color=col[i], linestyle='--')
            plt.plot(popSum[i][14], label=lab[i], color=col[i], linestyle='-')
        plt.legend(title=legendtitle)
        plt.ylabel('Buried ind.')
        plt.savefig(path + '/Ebola_B_' + savename + '.pdf', dpi=100)
        plt.show()
    # Unsafe funerals
    if type == 'all' or type == 'F':
        for i in range(0, len(names)):
            plt.plot(popSum[i][12], label=lab[i], color=col[i], linestyle='-')
        plt.legend(title=legendtitle)
        plt.ylabel('Unsafe funerals')
        plt.savefig(path + '/Ebola_F_' + savename + '.pdf', dpi=100)
        plt.show()
    # Recovered
    if type == 'all' or type == 'R':
        for i in range(0, len(names)):
            plt.plot(popSum[i][15], label=lab[i], color=col[i], linestyle='-')
        plt.legend(title=legendtitle)
        plt.ylabel('Recovered ind')
        plt.savefig(path + '/Ebola_R_' + savename + '.pdf', dpi=100)
        plt.show()

    if type == 'all' or type == 'mix':
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
        plt.legend(title=legendtitle)
        plt.ylabel('Individuals')
        plt.savefig(path + '/Ebola_mix_' + savename + '.pdf', dpi=100)
        plt.show()
    return(savename)


