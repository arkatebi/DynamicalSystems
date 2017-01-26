#/usr/bin/env python
import auxiliary_functions as aux
import PyDSTool as dst
from PyDSTool import common as cmn
import numpy as np
from matplotlib import pyplot as plt
import sys

# Simulation name
DSargs = cmn.args(name='Gene Interaction')

# Parameters
DSargs.pars = { 'gX': 5.0e1,
                'gY': 5.0e1,
                'X0': 1.0e2,
                'Y0': 1.0e2,
                'nX': 3,
                'nY': 3,
                'lX': 0.1,
                'lY': 0.1,
                'kX': 0.10e0, 
                'kY': 0.1e0}                

# auxiliary helper function(s) -- function name: ([func signature], definition)
DSargs.fnspecs={'HS': (['A','A0','nA','lamb'], 
                  'lamb + (1.0-lamb)/(1.0 + (A/A0)**nA)')} 

# rhs of the differential equation, including dummy variable
DSargs.varspecs={'X': 'gX*HS(Y,Y0,nY,lY) - kX*X', 
                 'Y' : 'gY*HS(X,X0,nX,lX) - kY*Y'} 
# initial conditions
DSargs.ics = {'X': 10, 'Y': 10}

DSargs.xdomain = {'X': [0, 1.0e+4], 'Y':[0, 1.0e+4]}
DSargs.tdomain = [0,100]    # set the range of integration.
ode = dst.Generator.Vode_ODEsystem(DSargs) # instance of the 'Generator' class

traj = ode.compute('polarization')  # integrate ODE
pts = traj.sample(dt=0.01)          # Data for plotting

def t_dynamics_X(): 
    # PyPlot commands
    #plt.plot(pts['X'], pts['Y'])
    plt.plot(pts['t'], pts['X'])
    #plt.plot(pts['t'], pts['Y'])

    plt.xlabel('t')       # Axes labels
    plt.ylabel('X')       # ...
    #plt.xlim([0,7000])   
    plt.ylim([0,200])     # Range of the y axis
    plt.title(ode.name)   # Figure title from model name
    plt.show()
    plt.figure()

def t_dynamics_Y(): 
    # PyPlot commands
    #plt.plot(pts['X'], pts['Y'])
    #plt.plot(pts['t'], pts['X'])
    plt.plot(pts['t'], pts['Y'])

    plt.xlabel('t')       # Axes labels
    plt.ylabel('Y')       # ...
    #plt.xlim([0,7000])   
    plt.ylim([0,200])     # Range of the y axis
    plt.title(ode.name)   # Figure title from model name
    plt.show()
    plt.figure()

def t_dynamics_XY(): 
    # PyPlot commands
    plt.plot(pts['X'], pts['Y'])
    #plt.plot(pts['t'], pts['X'])
    #plt.plot(pts['t'], pts['Y'])

    plt.xlabel('X')      # Axes labels
    plt.ylabel('Y')      # ...
    #plt.xlim([0,7000])
    plt.ylim([0,800])    # Range of the y axis
    plt.title(ode.name)  # Figure title from model name
    plt.show()

def t_dynamics_multi_ICs_X():   
    plt.ylim([0,200])
    plt.hold(True) # Sequences of plot commands will not clear existing figures
    for i, x0 in enumerate(np.linspace(-20,10,30)):
        ode.set(ics = { 'X': x0 } )    # Initial condition
        # Trajectories are called pol0, pol1, ...
        # sample them on the fly to create Pointset tmp
        tmp = ode.compute('pol%3i' % i).sample()
        plt.plot(tmp['t'], tmp['X'])
    plt.xlabel('time')
    plt.ylabel('X')
    plt.title(ode.name + ' multi ICs')
    plt.show()

def t_dynamics_multi_ICs_Y():   
    plt.ylim([0,200])
    plt.hold(True) # Sequences of plot commands will not clear existing figures
    for i, y0 in enumerate(np.linspace(-20,10,30)):
        ode.set(ics = { 'Y': y0 } )    # Initial condition
        # Trajectories are called pol0, pol1, ...
        # sample them on the fly to create Pointset tmp
        tmp = ode.compute('pol%3i' % i).sample()
        plt.plot(tmp['t'], tmp['Y'])
    plt.xlabel('time')
    plt.ylabel('Y')
    plt.title(ode.name + ' multi ICs')
    plt.show()

def t_dynamics_multi_ICs_XY_old():   
    plt.figure()
    plt.ylim([0,900]) 
    plt.hold(True) # Sequences of plot commands will not clear existing figures
    for i, v0 in enumerate(np.linspace(-20,10,30)):
        ode.set(ics = { 'X': v0, 'Y': v0 } )     # Initial condition
        # Trajectories are called pol0, pol1, ...
        # sample them on the fly to create Pointset tmp
        tmp = ode.compute('pol%3i' % i).sample()
        plt.plot(tmp['X'], tmp['Y'])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(ode.name + ' multi ICs XY')
    plt.show()

def t_dynamics_multi_ICs_X():
    plt.figure()
    plt.ylim([0,900])
    plt.hold(True) # Sequences of plot commands will not clear existing figures
    for i, x0 in enumerate(np.linspace(-20,10,30)):
        ode.set(ics = { 'X': x0 } )     # Initial condition
        # Trajectories are called pol0, pol1, ...
        # sample them on the fly to create Pointset tmp
        tmp = ode.compute('pol%3i' % i).sample()
        plt.plot(tmp['t'], tmp['X'])
    plt.xlabel('time')
    plt.ylabel('X')
    plt.title(ode.name + ' multi ICs X')
    plt.show()

def t_dynamics_multi_ICs_Y():
    plt.figure()
    plt.ylim([0,900])
    plt.hold(True) # Sequences of plot commands will not clear existing figures
    for i, y0 in enumerate(np.linspace(-20,10,30)):
        ode.set(ics = { 'Y': y0 } )   # Initial condition
        # Trajectories are called pol0, pol1, ...
        # sample them on the fly to create Pointset tmp
        tmp = ode.compute('pol%3i' % i).sample()
        plt.plot(tmp['t'], tmp['Y'])
    plt.xlabel('time')
    plt.ylabel('Y')
    plt.title(ode.name + ' multi ICs Y')
    plt.show()

def t_dynamics_multi_ICs_XY():   
    plt.figure()
    plt.ylim([0,900])
    plt.hold(True) # Sequences of plot commands will not clear existing figures
    for i, x0 in enumerate(np.linspace(1,1000,4)):
        for i, y0 in enumerate(np.linspace(1,1000,4)):
            ode.set(ics = { 'X': x0, 'Y': y0 } )    # Initial condition
            # Trajectories are called pol0, pol1, ...
            # sample them on the fly to create Pointset tmp
            tmp = ode.compute('pol%3i' % i).sample()
            plt.plot(tmp['X'], tmp['Y'])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(ode.name + ' multi ICs XY')
    plt.show()

def getBifDiagrams():
    ode.set(pars={'gX':5.0e1, 'gY':5.0e1,
                  'X0':1.0e2, 'Y0': 1.0e2,
                  'nX':3, 'nY':3,
                  'lX':0.1,'lY':0.1,
                  'kX':0.10e0, 'kY':0.1e0})
    ode.set(ics = {'X': 100, 'Y': 50})
    freepar='gX'
    fp=aux.fast_fixedpoint(ode)
    print(fp.values())
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=4, 
                          LocBifPoints=['Lp','B'], bif_startpoint=50, 
                          maxstep=1e+1, minstep=0.01, step=0.1, 
                          silence=True, fs=[4,4], ics=[fp], xlim=[0,300], 
                          ylim=[0,800], fontsize=10)
    sys.exit(0)
    freepar='gY'
    fp=aux.fast_fixedpoint(ode)
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=4, 
                          LocBifPoints=['Lp','B'], bif_startpoint=10, 
                          maxstep=1e+4, minstep=1e-1, step=5e+2, 
                          silence=True, fs=[6,5], ics=[fp], xlim=[0,1000], 
                          fontsize=10)
 
def getNullClines(): 
    ode.set(pars={'gX':5.0e1, 'gY':5.0e1,
                  'X0':1.0e2, 'Y0':1.0e2,
                  'nX':3, 'nY':3,
                  'lX':0.1,'lY':0.1,
                  'kX':0.10e0, 'kY':0.1e0})
    #ode.set(ics = {'X': 30, 'Y': 1})
 
    from PyDSTool.Toolbox import phaseplane as pp
    vlim = {'X': [1, 1000], 'Y': [1, 1000]}
    fp = aux.eliminate_redundants(pp.find_fixedpoints(ode, n=2, maxsearch=1e+4,
                                                     eps=1e-12),
                                                     4)
    stab = aux.stability(fp, ode)
    
    for i in range(len(fp)):
        print(stab[i], fp[i])
    nfp=0
    aux.nullclines(['X','Y'], DSargs, stab, fp, nfp=nfp, vlim=vlim,
                   maxpoints=[1000,1000],
                   xticks=[0, 100, 200, 300, 400, 500, 600, 700, 800],
                   yticks=[0, 100, 200, 300, 400, 500, 600, 700, 800],
                   step=0.01, minstep=0.001, maxstep=10, fs=[4,4], 
                   fontsize=8, silence=False)

if __name__ == '__main__': 
    #t_dynamics_X()
    #t_dynamics_Y()
    #t_dynamics_XY()

    #t_dynamics_multi_ICs_X()
    #t_dynamics_multi_ICs_Y()
    #t_dynamics_multi_ICs_XY()

    #getBifDiagrams()
    getNullClines()
