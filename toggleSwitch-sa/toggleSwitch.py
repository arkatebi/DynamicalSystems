#/usr/bin/env python
import auxiliary_functions as aux
import PyDSTool as dst
from PyDSTool import common as cmn
import numpy as np
from matplotlib import pyplot as plt
import sys

def defineSystem(): 
    # Create an object of args class from common module 
    DSargs = cmn.args(name='Genetic Toggle Switch with SA')

    # Initialize the DSargs object with parameters
    #DSargs.pars = aux.parameter_set_1()
    #DSargs.pars = aux.parameter_set_2()
    #DSargs.pars = aux.parameter_set_3() # gives interesting intersection
    DSargs.pars = aux.parameter_set_4()

    # obtain the differential equations:
    DSargs.varspecs = aux.equations() 
    # obtain the auxiliary functions:
    DSargs.fnspecs = aux.functions()

   # Set initial conditions:
    #DSargs.ics = {'X': 10, 'Y': 10}
    DSargs.ics = {'X': 10, 'Y': 10}

    DSargs.xdomain = {'X': [0, 1.0e+4], 'Y':[0, 1.0e+4]}

    # Set the range of integration:
    DSargs.tdomain = [0,100]
    return DSargs

def t_dynamics_X(pts): 
    # PyPlot commands
    plt.plot(pts['t'], pts['X'])

    plt.xlabel('t')       # Axes labels
    plt.ylabel('X')       # ...
    #plt.xlim([0,7000])   
    plt.ylim([0,600])     # Range of the y axis
    plt.title(ode.name + ': time dynamics')   # Figure title from model name
    plt.show()
    plt.figure()

def t_dynamics_Y(pts): 
    # PyPlot commands
    plt.plot(pts['t'], pts['Y'])

    plt.xlabel('t')       # Axes labels
    plt.ylabel('Y')       # ...
    #plt.xlim([0,7000])   
    plt.ylim([0,600])     # Range of the y axis
    plt.title(ode.name + ': time dynamics')   # Figure title from model name
    plt.show()
    plt.figure()

def t_dynamics_XY(pts): 
    # PyPlot commands
    plt.plot(pts['X'], pts['Y'])

    plt.xlabel('X')      # Axes labels
    plt.ylabel('Y')      # ...
    #plt.xlim([0,7000])
    plt.ylim([0,800])    # Range of the y axis
    plt.title(ode.name + ': time dynamics')   # Figure title from model name
    plt.show()

def t_dynamics_multi_ICs_X(ode):   
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
    plt.title(ode.name + ': multi ICs')
    plt.show()

def t_dynamics_multi_ICs_Y(ode):   
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

def t_dynamics_multi_ICs_X(ode):
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

def t_dynamics_multi_ICs_Y(ode):
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

def t_dynamics_multi_ICs_XY(ode):   
    plt.figure()
    #plt.ylim([0,600])
    plt.ylim([0,1200])
    # Sequences of plot commands will not clear existing figures:
    plt.hold(True) 
    for i, x0 in enumerate(np.linspace(1,1000,4)):
        for i, y0 in enumerate(np.linspace(1,1000,4)):
            # Reset the initial conditions in the Vode_ODEsystem object ode:
            ode.set(ics = { 'X': x0, 'Y': y0 } )    
            # Trajectories are called pol0, pol1, ...
            # Sample them on the fly to create tmp, a Pointset object: 
            tmp = ode.compute('pol%3i' % i).sample()
            plt.plot(tmp['X'], tmp['Y'])
    plt.xlabel('X')
    plt.ylabel('Y')
    #plt.title(ode.name + ': multi ICs for both')
    plt.show()
    #plt.savefig('./figures/parSet-1_tdynamics.pdf')

def getBifDiagrams(ode):
    freepar='gX'
    fp=aux.fast_fixedpoint(ode)
    print(fp.values())
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=1, 
                          LocBifPoints=['LP','B'], bif_startpoint=50, 
                          maxstep=1e+1, minstep=0.01, step=0.1,
                          silence=True, fs=[6,6], ics=[fp], 
                          xlim=[0,600], ylim=[0,6000], fontsize=10)
    freepar='gY'
    fp=aux.fast_fixedpoint(ode)
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=1, 
                          LocBifPoints=['LP','B'], bif_startpoint=50, 
                          maxstep=1e+1, minstep=1e-2, step=1e-1, 
                          silence=True, fs=[6,6], ics=[fp], 
                          xlim=[0,600], ylim=[0,6000], fontsize=10)
    sys.exit(0)
 
    freepar='lX'
    fp=aux.fast_fixedpoint(ode)
    print(fp.values())
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=1, 
                          LocBifPoints=['LP','B'], bif_startpoint=50, 
                          maxstep=1e+1, minstep=0.01, step=0.1,
                          silence=True, fs=[6,6], ics=[fp], 
                          xlim=[0,200], ylim=[0,700], fontsize=10)
    freepar='lY'
    fp=aux.fast_fixedpoint(ode)
    print(fp.values())
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=1, 
                          LocBifPoints=['LP','B'], bif_startpoint=50, 
                          maxstep=1e+1, minstep=0.01, step=0.1,
                          silence=True, fs=[6,6], ics=[fp], 
                          xlim=[0,200], ylim=[0,700], fontsize=10)
    sys.exit(0)
    freepar='kX'
    fp=aux.fast_fixedpoint(ode)
    print(fp.values())
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=1, 
                          LocBifPoints=['LP','B'], bif_startpoint=50, 
                          maxstep=1e+1, minstep=0.01, step=0.1,
                          silence=True, fs=[6,6], ics=[fp], 
                          xlim=[0,200], ylim=[0,700], fontsize=10)
    freepar='kY'
    fp=aux.fast_fixedpoint(ode)
    print(fp.values())
    aux.plot_continuation(ode, freepar, keys=['X','Y'], ncol=2, nrow=1, 
                          LocBifPoints=['LP','B'], bif_startpoint=50, 
                          maxstep=1e+1, minstep=0.01, step=0.1,
                          silence=True, fs=[6,6], ics=[fp], 
                          xlim=[0,200], ylim=[0,700], fontsize=10)


def getNullClines(DSargs, ode): 
    from PyDSTool.Toolbox import phaseplane as pp
    vlim = {'X': [1, 1400], 'Y': [1, 400]}
    fp = aux.eliminate_redundants(pp.find_fixedpoints(ode, n=2, maxsearch=1e+4,
                                                     eps=1e-12),
                                                     4)
    stab = aux.stability(fp, ode)
     
    for i in range(len(fp)):
        print(stab[i], fp[i])
    nfp=0
    aux.nullclines(['X','Y'], DSargs, stab, fp, nfp=nfp, vlim=vlim,
                   maxpoints=[400,400],
                  #xticks=[0, 200, 400, 600, 800, 1000, 1200, 1400],
                   xticks=[0, 200, 400],
                  #yticks=[0, 100, 200, 300, 400, 500, 600, 700,800, 900, 1000],
                   yticks=[0, 100, 200, 300, 400],
                   step=0.01, minstep=0.001, maxstep=10, fs=[4,4], 
                   fontsize=8, silence=False)

if __name__ == '__main__': 
    DSargs = defineSystem()
    # Obtain a Vode_ODEsystem object: 
    # (similar to VODE from SciPy) 
    ode = dst.Generator.Vode_ODEsystem(DSargs) 
    # Obtain a Trajectory object (integrate ODE): 
    traj = ode.compute('polarization')  
    # Collect data points as a Pointset object:
    pts = traj.sample(dt=0.01)          
     
    #t_dynamics_X(pts)
    #t_dynamics_Y(pts)
    #t_dynamics_XY(pts)

    #t_dynamics_multi_ICs_X(ode)
    #t_dynamics_multi_ICs_Y(ode)
    #t_dynamics_multi_ICs_XY(ode)

    #getBifDiagrams(ode)
    getNullClines(DSargs, ode)
