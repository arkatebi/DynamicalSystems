#/usr/bin/env python
'''
    Toggle Switch consisting of two mutually inhibitory genes: 
    X ---| Y; Y ---| X
    DEs for the toggle switch system: 
    dX/dt = gX*HS(Y) - kX*X
    dY/dt = gY*HS(X) - kY*Y

    Stochasticity of the switch is observed by applying Gillespie's 
    Direct Method.
    The equivalent chemical system for the switch: 
    R1: phi -> X, with reaction (production) rate gX*HS(Y) 
    R2: X -> phi, with reaction (degradation) rate kX
    R3: phi -> Y with reaction (production) rate gY*HS(X) 
    R4: Y -> phi, with reaction (degradation) rate kY
    
    System variables: X and Y
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
    pars = aux.parameter_st_1()
    #set initial conditions for the system variables: 
    vars=OrderedDict()
    vars['X'] = 100 
    vars['Y'] = 100

    #set simulation time: 
    tmax=5.0e+6 
    return (pars,vars,tmax)  
    
#-----------------------------------------------------------------------------#
def shifted_hill_impact(X, X0, nX, lamb):
    '''
    This method calculates and returns the effect of shifted Hill function.
    '''
    return lamb + (1.0-lamb)/(1.0+(X/X0)**nX)

#-----------------------------------------------------------------------------#
def calculate_propensities(pars, vars):
    '''
    Calculates propensitity for each reaction.
    '''
    #create an empty list to store propensities:
    pros = list()
    sh_impact_X = shifted_hill_impact(vars.get('X'), 
                                      pars.get('X0'), 
                                      pars.get('nX'), 
                                      pars.get('lX'))
    sh_impact_Y = shifted_hill_impact(vars.get('Y'), 
                                      pars.get('Y0'), 
                                      pars.get('nY'), 
                                      pars.get('lY'))
    #add propensity for reaction 1:
    pros.append(pars.get('gX')*sh_impact_Y)
    #add propensity for reaction 2:
    pros.append(pars.get('kX')*vars.get('X'))
    #add propensity for reaction 3:
    pros.append(pars.get('gY')*sh_impact_X)
    #add propensity for reaction 4:
    pros.append(pars.get('kY')*vars.get('Y'))
    return pros 

#-----------------------------------------------------------------------------#
def perform_sp_reaction(vars, pros):
    '''
    This method stochastically chooses which reaction to occur and updates the 
    system after performing the chosen reaction.
    '''
    ptotal = sum(pros) #total propensity
    pr = [x/ptotal for x in pros]
    rn=random.uniform(0,1)
    #choose reaction based on rn:
    if (rn<=pr[0]):
        vars['X']+=1
    elif (rn<=(pr[0]+pr[1])):
        vars['X']-=1
    elif (rn<=(pr[0]+pr[1]+pr[2])):
        vars['Y']+=1
    else: #(rn>(pr[0]+pr[1]+pr[2]) and pr<=(pr[0]+pr[1]+pr[2]+pr[3])):
        vars['Y']-=1
    return vars 

#-----------------------------------------------------------------------------#
def run_simulation(pars, vars, tmax):
    '''
    This method repeatedly performs the following tasks: 
    (1) Calculate propensities for each reaction as a function of the current
        state of the system.
    (2) Stochastically chooses the wait time (dt) for the next reaction to 
        occur.  
    (3) Stochastically selects the next reaction to occur and update the 
        system after the reaction is performed.
    '''
    import numpy as np

    #initialize current time:
    tc=0
    #iteration counter:
    factor=1000
    count = 0
    #run the simulation using Gillespie's Direct Method:
    while(tc<tmax):
        #calculate propensities:
        pros = calculate_propensities(pars, vars)
        ptotal=sum(pros) #total propensity
        if(not ptotal):
            continue
        #save configuration at multiple of 'factor' timesteps:
        if (not (count%factor)):
            print(tc, ' ', vars.get('X'), ' ', vars.get('Y'))

        #perform specific reaction based on propensity values:
        vars=perform_sp_reaction(vars, pros)
        if(vars.get('X')<0 or vars.get('Y')<0):
            return None
        #obtain wait time dt from exponential distribution:
        lambd=ptotal # lambd = 1/mean where mean = 1/Rtotal
        dt=random.expovariate(lambd)
        tc+=dt
        count+=1
    return None

def run_simulation_old(pars, vars, tmax):
    '''
    This method repeatedly performs the following tasks: 
    (1) Calculate propensities for each reaction as a function of the current
        state of the system.
    (2) Stochastically chooses the wait time (dt) for the next reaction to 
        occur.  
    (3) Stochastically selects the next reaction to occur and update the 
        system after the reaction is performed.
    '''
    import numpy as np

    #initialize current time:
    tc=0
    #iteration counter:
    factor=1000
    count = 0
    #run the simulation using Gillespie's Direct Method:
    while(tc<tmax):
        #calculate propensities:
        pR1,pR2,pR3,pR4=calculate_propensities(pars, vars)
        #pros = calculate_propensities(pars, vars)
        Rtotal= pR1+pR2+pR3+pR4
        #ptotal=sum(pros) #total propensity
        #if(not ptotal):
        if(not Rtotal):
            continue
        #save configuration at multiple of 'factor' timesteps:
        if (not (count%factor)):
            print(tc, ' ', vars.get('X'), ' ', vars.get('Y'))

        #perform specific reaction based on propensity values:
        vars=perform_sp_reaction_old(vars,pR1/Rtotal,pR2/Rtotal,
                                      pR3/Rtotal,pR4/Rtotal)
        #vars=perform_sp_reaction(vars, pros)

        if(vars.get('X')<0 or vars.get('Y')<0):
            return None
        #obtain wait time dt from exponential distribution:
        lambd=Rtotal # lambd = 1/mean where mean = 1/Rtotal
        #lambd=ptotal # lambd = 1/mean where mean = 1/Rtotal
        dt=random.expovariate(lambd)
        tc+=dt
        count+=1
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
