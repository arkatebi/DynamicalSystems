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

#-----------------------------------------------------------------------------#
def shifted_hill_impact(X, X0, nX, lamb):
    '''
    This method calculates and returns the effect of shifted Hill function.
    '''
    return lamb + (1.0-lamb)/(1.0+(X/X0)**nX)

#-----------------------------------------------------------------------------#
def calculate_propensities(X0,Y0,
                           gX,gY,
                           kX,kY,
                           X,Y,
                           nX,nY,
                           lX,lY):
    sh_impact_X = shifted_hill_impact(X, X0, nX, lX)
    sh_impact_Y = shifted_hill_impact(Y, Y0, nY, lY)
    #pR1 = gX*sh_impact_Y*X
    pR1 = gX*sh_impact_Y
    pR2 = kX*X
    #pR3 = gY*sh_impact_X*Y
    pR3 = gY*sh_impact_X
    pR4 = kY*Y
    return pR1,pR2,pR3,pR4

#-----------------------------------------------------------------------------#
def perform_sp_reaction(X, Y, p1, p2, p3, p4):
    '''
    This method stochastically chooses which reaction to occur and updates the 
    system after performing the chosen reaction.
    '''
    pr=random.uniform(0,1)

    if (pr<=p1):
        X+=1
    elif (pr<=(p1+p2)):
        X-=1
    elif (pr<=(p1+p2+p3)):
        Y+=1
    else: #(pr>(p1+p2+p3) and pr<=(p1+p2+p3+p4)):
        Y-=1
    #else:
    #    pass 
    return (X,Y) 

#-----------------------------------------------------------------------------#
def run_simulation(X0,Y0,
                   gX,gY,
                   kX,kY,
                   X,Y,
                   nX,nY,
                   lX,lY,
                   tmax):
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
    cnt=0;
    #initialize variables to save the trajectory: 
    ser_X=[X]
    ser_Y=[Y]
    ser_cnt=[0]
    ser_tm=[tc]
    factor=1000
    count = 0
    #run the simulation using Gillespie's Direct Method:
    while(tc<tmax):
        #calculate propensities:
        pR1,pR2,pR3,pR4=calculate_propensities(X0,Y0,
                                               gX,gY,
                                               kX,kY,
                                               X,Y,
                                               nX,nY,
                                               lX,lY)
        Rtotal= pR1+pR2+pR3+pR4
        if(not Rtotal):
            continue
        #if (int(tc)/1000)
        if (not (count%factor)):
            print(tc, ' ', X, ' ', Y)
        #perform specific reaction based on propensity values:

        (X,Y)=perform_sp_reaction(X,Y,pR1/Rtotal,pR2/Rtotal,
                                          pR3/Rtotal,pR4/Rtotal)
        if(X<0 or Y<0):
            return (ser_cnt,ser_X,ser_Y) 

        #parameter for the exponential distribution:
        lambd=Rtotal # lambd = 1/mean where mean = 1/Rtotal
        dt=random.expovariate(lambd)
        tc+=dt
        count+=1
        #cnt+=1
        #store the state of the system:
        #ser_X.append(X)
        #ser_Y.append(Y)
        #ser_cnt.append(cnt)
    return ser_cnt,ser_X,ser_Y

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
    #set system parameters:
    #X0=100; Y0=100
    #X0=30; Y0=30
    #X0=20; Y0=20
    X0=10; Y0=10
    #gX=50; gY=50; #production rates for X and Y 
    #gX=10; gY=10; #production rates for X and Y 
    gX=5; gY=5; #production rates for X and Y 
    kX=0.1; kY=0.1 #degradation rates for X and Y
    nX=3.0; nY=3.0 #cooperativities of X and Y
    lX=0.1; lY=0.1 #fold changes for the shifted Hill function

    #set initial conditions:
    X=100; Y=100
    #X=30; Y=30
    #set reaction time:
    tmax=5000000
    t_start=time.time()
    #run simulation:
    (ser_cnt,ser_X,ser_Y)=run_simulation(X0,Y0,
                                           gX,gY,
                                           kX,kY,
                                           X,Y,
                                           nX,nY,
                                           lX,lY,
                                           tmax)
    t_end=time.time()
    print('simulation time: ', t_end-t_start)
    #plot_trajectory(ser_cnt,ser_X,ser_Y)
