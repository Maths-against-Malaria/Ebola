
def plotEbola(names, savename, pathIn=pathIn, pathOut = pathOut, col=colsA, lab=lab, type='all', q_max=False, Nerls=Nerls, legendtitle=''):

    index = indexFunction(Nerls)
    popSum0 = popsum2d(pops=np.loadtxt(pathIn + '/ebola_' + names[0] + '.txt'), Nerls=Nerls)
    s = np.shape(popSum0)
    if q_max == True:
        q_ = np.empty(shape = [len(names),s[1]])
        for i in range(0, len(names)):
            q_[i] = getQ(pathIn=pathIn, name = names[i])
        if type == 'q':
            for i in range(0,len(names)):
                plt.plot(q_[i],color=col[i], label = lab[i],)
            plt.ylim([-0.1,1.1])
            plt.legend()
            plt.ylabel('q')
            plt.savefig (pathOut + '/Ebola_q_' + savename + '.pdf', dpi=100)
            plt.show()

    popSum = np.empty(shape = [len(names),s[0], s[1]])
    popSum[0] = popSum0
    for i in range(0,len(names)):
        pops_i = np.loadtxt(pathIn + '/ebola_' + names[i] + '.txt')
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
        plt.savefig (pathOut + '/Ebola_S_' + savename + '.pdf', dpi=100)
        plt.show()

    # Latent
    if type == 'all' or type == 'E':
        plt.plot(popSum[0][1], color=col[0], linestyle='-')
        plt.plot(popSum[0][1], color=col[0], linestyle='-', label='never traced back')
        plt.plot(popSum[0][3], color=col[0], linestyle='--', label='not yet traced back')
        plt.plot(popSum[0][2], color=col[0], linestyle=':', label='traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][1], label=lab[i], color=col[i], linestyle='-')
            plt.plot(popSum[i][2], color=col[i], linestyle=':')
            plt.plot(popSum[i][3], color=col[i], linestyle='--')
        plt.legend(title=legendtitle)
        plt.ylabel('Latent ind.')
        plt.savefig(pathOut + '/Ebola_E_' + savename + '.pdf', dpi=100)
        plt.show()
    # Podromal
    if type == 'all' or type == 'P':
        plt.plot(popSum[0][4], color=col[0], linestyle='-')
        plt.plot(popSum[0][4], color=col[0], linestyle='-', label='never traced back')
        plt.plot(popSum[0][6], color=col[0], linestyle='--', label='not yet traced back')
        if q_max == True:
            plt.plot(np.multiply(popSum[0][5], q_[0]), color=col[0], linestyle=':', label='traced back, in ward')
            plt.plot(np.multiply(popSum[0][5], np.multiply(q_[0],-1)+1), color=col[0], linestyle='-.', label='traced back, not in ward')
        else:
            plt.plot(popSum[0][5], color=col[0], linestyle=':', label='traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][4], label=lab[i], color=col[i], linestyle='-')
            plt.plot(popSum[i][6], color=col[i], linestyle='--')
            if q_max == True:
                plt.plot(np.multiply(popSum[i][5], q_[i]), color=col[i], linestyle=':')
                plt.plot(np.multiply(popSum[i][5], np.multiply(q_[i],-1)+1), color=col[i], linestyle='-.')
            else:
                plt.plot(popSum[i][5], color=col[i], linestyle=':')

        plt.legend(title=legendtitle)
        plt.ylabel('Podromal ind.')
        plt.savefig(pathOut + '/Ebola_P_' + savename + '.pdf', dpi=100)
        plt.show()
    # Fully infected at home
    if type == 'all' or type == 'Ip':
        plt.plot(popSum[0][7], color=col[0], linestyle='-')
        plt.plot(popSum[0][7], color=col[0], linestyle='-', label='never traced back')
        plt.plot(popSum[0][11], color=col[0], linestyle='--', label='not yet traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][7], label=lab[i], color=col[i], linestyle='-')
            plt.plot(popSum[i][11], color=col[i], linestyle=':')
        plt.legend(title=legendtitle)
        plt.ylabel('Fully inf. ind. at home')
        plt.savefig(pathOut + '/Ebola_Ip_' + savename + '.pdf', dpi=100)
        plt.show()
    # Fully infected in hospital
    if type == 'all' or type == 'Ih':
        plt.plot(popSum[0][8], color=col[0], linestyle='-')
        plt.plot(popSum[0][8], color=col[0], linestyle='-', label='never traced back')
        plt.plot(popSum[0][10], color=col[0], linestyle='--', label='not yet traced back')
        for i in range(0, len(names)):
            plt.plot(popSum[i][8], label=lab[i], color=col[i], linestyle='-')
            plt.plot(popSum[i][10], color=col[i], linestyle=':')
        plt.legend(title=legendtitle)
        plt.ylabel('Fully inf. ind. in hospital')
        plt.savefig(pathOut + '/Ebola_Ih_' + savename + '.pdf', dpi=100)
        plt.show()
    # Fully infected in isolation
    if type == 'all' or type == 'Ii':
        if q_max == True:
            plt.plot(np.multiply(popSum[0][9], q_[0]), color=col[0], linestyle='-', label='in ward')
            plt.plot(np.multiply(popSum[0][9], np.multiply(q_[0],-1)+1), color=col[0], linestyle='-.', label='traced back, not in ward')
            for i in range(0, len(names)):
                plt.plot(np.multiply(popSum[i][9], q_[i]), color=col[i], linestyle='-')
                plt.plot(np.multiply(popSum[i][9], np.multiply(q_[i], -1) + 1), color=col[i], linestyle='-.')
        if q_max == False:
            plt.plot(popSum[0][9], label=lab[0], color=col[0], linestyle='-')
            for i in range(0, len(names)):
                plt.plot(popSum[i][9], label=lab[i], color=col[i], linestyle='-')
        plt.ylim([-5,np.ndarray.max(popSum[:,9])*1.05])
        plt.legend(title=legendtitle)
        plt.ylabel('Fully inf. ind. in isolation')
        plt.savefig(pathOut + '/Ebola_Ii_' + savename + '.pdf', dpi=100)
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
        plt.savefig(pathOut + '/Ebola_B_' + savename + '.pdf', dpi=100)
        plt.show()
    # Unsafe funerals
    if type == 'all' or type == 'F':
        for i in range(0, len(names)):
            plt.plot(popSum[i][12], label=lab[i], color=col[i], linestyle='-')
        plt.legend(title=legendtitle)
        plt.ylabel('Unsafe funerals')
        plt.savefig(pathOut + '/Ebola_F_' + savename + '.pdf', dpi=100)
        plt.show()
    # Recovered
    if type == 'all' or type == 'R':
        for i in range(0, len(names)):
            plt.plot(popSum[i][15], label=lab[i], color=col[i], linestyle='-')
        plt.legend(title=legendtitle)
        plt.ylabel('Recovered ind')
        plt.savefig(pathOut + '/Ebola_R_' + savename + '.pdf', dpi=100)
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
        plt.savefig(pathOut + '/Ebola_mix_' + savename + '.pdf', dpi=100)
        plt.show()
    return(savename)


plt.rcParams['legend.framealpha'] = 0
