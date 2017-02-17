#/usr/bin/env python
'''
    Mutually Inhibitory Toggle Switch with protein burst.
    X ---| Y; Y ---| X
    DEs for the toggle switch system:
    dX/dt = gX*HS(Y) - kX*X
    dY/dt = gY*HS(X) - kY*Y

    The equivalent chemical system: 
    R1: phi -> X, with reaction (production) rate gX*HS(Y) 
    R2: X -> phi, with reaction (degradation) rate kX
    R3: phi -> Y with reaction (production) rate gY*HS(X) 
    R4: Y -> phi, with reaction (degradation) rate kY
    
    System variables: X, Y, Px, and Py
    Dummy variable: phi
    System parameters: gX, gY, kX, kY 
    Auxiliary function: HS    
'''
import sys
import random
from collections import OrderedDict
import auxiliary_functions as aux

#-----------------------------------------------------------------------------#
def defineSystem():
    '''
    Create dictionaries to define the system
    '''
    #set parameters:
    pars = OrderedDict() 
    pars = aux.parameter_burst_1()
    #set initial conditions for the system variables: 
    vars=OrderedDict()
    vars['X'] = 100 
    vars['Y'] = 100
    vars['Px'] = 0 
    vars['Py'] = 1
    #set simulation time: 
    tmax=5.0e+6
    return (pars,vars,tmax)

#-----------------------------------------------------------------------------#
def calculate_propensities(pars, vars):
    #create an empty list to store propensities:
    pros = list()

    #reactions of X system: Promoter P on
    #R1:
    pros.append(pars.get('KoffX'))
    #R2:
    pros.append(pars.get('gXon'))
    #R3:
    pros.append(pars.get('kX')*vars.get('X'))

    #reactions of X system: Promoter P off
    #R1:
    pros.append(pars.get('KonX')*(vars.get('Y')**pars.get('nY')))
    #R2:
    pros.append(pars.get('gXoff'))
    #R3:
    pros.append(pars.get('kX')*vars.get('X'))

    #reactions of Y system: Promoter P on
    #R1:
    pros.append(pars.get('KoffY'))
    #R2:
    pros.append(pars.get('gYon'))
    #R3:
    pros.append(pars.get('kY')*vars.get('Y'))

    #reactions Y system: Promoter P off
    #R1:
    pros.append(pars.get('KonY')*(vars.get('X')**pars.get('nX')))
    #R2:
    pros.append(pars.get('gYoff'))
    #R3:
    pros.append(pars.get('kY')*vars.get('Y'))
    return pros 

#-----------------------------------------------------------------------------#
def updateSystem_PonX(pars, vars, pros):
    '''
    This method updates the system when X promoter P is on.
    '''
    #print(pros)
    #print(len(pros))
    #print('PonX')
    #sys.exit(0)
 
    ptotal = sum(pros) #total propensity
    pr = [x/ptotal for x in pros]
    rn=random.uniform(0,1)
    if rn<pr[0]: 
        vars['Y']+=pars['nY']
        vars['Px']=0
    elif rn <pr[0]+pr[1]:   
        vars['X']+=1
    else:
        vars['X']-=1
    return vars

#-----------------------------------------------------------------------------#
def updateSystem_PoffX(pars, vars, pros):
    '''
    This method updates the system when X promoter P is off.
    '''
    #print(pros)
    #print(len(pros))
    #print('PoffX')
    #sys.exit(0)
 
    ptotal = sum(pros) #total propensity
    pr = [x/ptotal for x in pros]
    rn=random.uniform(0,1)
    if rn<pr[0] and vars.get('Y')>=pars.get('nY'):
        vars['Y']-=pars['nY']
        vars['Px']=1 # Promoter switches to ON state
    elif rn <pr[0]+pr[1]:   
        vars['X']+=1
    elif rn <pr[0]+pr[1]+pr[2]:   
        vars['X']-=1
    return vars

#-----------------------------------------------------------------------------#
def updateSystem_PonY(pars, vars, pros):
    '''
    This method updates the system when Y promoter P is on.
    '''
    #print(pros)
    #print(len(pros))
    #print('PonY')
    #sys.exit(0)
 
    ptotal = sum(pros) #total propensity
    pr = [x/ptotal for x in pros]
    rn=random.uniform(0,1)
    if rn<pr[0]: 
        vars['X']+=pars['nX']
        vars['Py']=0
    elif rn <pr[0]+pr[1]:   
        vars['Y']+=1
    else:   
        vars['Y']-=1
    return vars

#-----------------------------------------------------------------------------#
def updateSystem_PoffY(pars, vars, pros):
    '''
    This method updates the system when Y promoter P is off.
    '''
    #print(pros)
    #print(len(pros))
    #print('PoffY')
    #sys.exit(0)
 
    ptotal = sum(pros) #total propensity
    pr = [x/ptotal for x in pros]
    rn=random.uniform(0,1)
    if rn<pr[0] and vars.get('X')>=pars.get('nX'):
        vars['X']-=pars['nX']
        vars['Py']=1 #Promoter switches to ON state
    elif rn<pr[0]+pr[1]:
        vars['Y']+=1
    elif rn<pr[0]+pr[1]+pr[2]:
        vars['Y']-=1
    return vars

#-----------------------------------------------------------------------------#
def updateSystem(pars, vars, pros):
    '''
    This method updates system after performing specific 
    reaction in each subsystem.
    '''
    #print(pros)
    #print(len(pros))
    if vars.get('Px'):
        updateSystem_PonX(pars, vars, pros[0:3])
    elif not vars.get('Px') and vars.get('Y')>=pars.get('nY'):
        updateSystem_PoffX(pars, vars, pros[3:6])
    if vars.get('Py'):
        updateSystem_PonY(pars, vars, pros[6:9])
    elif not vars.get('Py') and vars.get('X')>=pars.get('nX'):
        updateSystem_PoffY(pars, vars, pros[9:])
    return vars

#-----------------------------------------------------------------------------#
def run_simulation(pars, vars, tmax):
    '''
    This method repeatedly performs the following tasks: 
    (1) Calculate propensities for each reaction as a function of the current
        state of the system.
    (2) Stochastically selects the next reaction to occur and update the 
        system after the reaction is performed.
    (3) Stochastically chooses the wait time (dt) for the next reaction to 
        occur.
    '''
    import numpy as np

    #initialize current time:
    tc=0
    #iteration counter:
    factor=1000
    count = 0
    #run the simulation using Gillespie's Direct Method:
    print('tc', '\t', 'X', '\t', 'Y', '\t', 'Px', '\t', 'Py')
    while(tc<tmax):
        #calculate propensities:
        pros = calculate_propensities(pars, vars)
        #print(pros)
        #print(len(pros))
        ptotal=sum(pros) #total propensity
        if(not ptotal):
            continue
        #save configuration at multiple of 'factor' timesteps:
        if (not (count%factor)):
            print(tc, '\t', vars.get('X'), '\t', vars.get('Y'), 
                      '\t', vars.get('Px'), '\t', vars.get('Py'))

        #perform specific reaction based on propensity values:
        vars=updateSystem(pars, vars, pros)
        #if(vars.get('X')<0 or vars.get('Y')<0):
        #    print(tc, '\t', vars.get('X'), '\t', vars.get('Y'), 
        #              '\t', vars.get('Px'), '\t', vars.get('Py'))
        #    return None
        #obtain wait time dt from exponential distribution:
        lambd=ptotal # lambd = 1/mean where mean = 1/Rtotal
        dt=random.expovariate(lambd)
        tc+=dt
        count+=1
    #print(count)
    return None

#-----------------------------------------------------------------------------#
def plot_trajectory(series_cnt, series_X, series_Y):
    import pylab as pl 
    xLen=int(len(series_cnt)/1)
    #pl.plot(series_cnt, series_X, ':k', label='X') # plot X trajectory
    pl.plot(series_cnt[0:xLen], series_X[0:xLen], ':k', label='X') #X trajectory
    #pl.plot(series_cnt, series_Y, '-r', label='Y') # plot Y trajectory
    pl.plot(series_cnt[0:xLen], series_Y[0:xLen], '-r', label='Y') #Y trajectory
    pl.legend()
    pl.show()

#-----------------------------------------------------------------------------#
if __name__=="__main__":
    import time
    (pars,vars,tmax)=defineSystem()
    t_start=time.time()
    run_simulation(pars, vars, tmax)
    t_end=time.time()
    print('simulation time: ', t_end-t_start)
    #plot_trajectory(ser_cnt,ser_X,ser_Y)
