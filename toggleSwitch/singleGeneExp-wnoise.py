#/usr/bin/env python
''''
Simulation of the birth-death process of a single gene considering the 
intrinsic noise in the birth-death process.
 
ODE: NO intrinsic noise is prsent
dX/dt=gX-kX*X

SDE: intrinsic noise is present
dX/dt=gX-kX*X+g*sqrt(dt)*dW(dt)
where W models the standard Brownian motion
      g models the strength of the intrinsic noise
'''

import sys
import numpy as np
from collections import OrderedDict

#-------------------------------------------------------------------------------
def defineSystem():
    #initialize the system variables:
    #variable sets:
    vars=OrderedDict()
    vars['X']=10

    pars=OrderedDict()
    #production rate:
    pars['gX']=5.0e0 
    #degration rate:
    pars['kX']=1.0e-1 

    #time step size:
    dt=1.0e-3 
    #maximum simulation time:
    tmax=1.0e+3
    return (vars,pars,dt,tmax) 

#-------------------------------------------------------------------------------
def getEffectiveRate(X, gX, kX):
    return gX-kX*X
        
#-------------------------------------------------------------------------------
def getNoiseStrength(X,gX,kX):
    g=np.sqrt(abs(gX+kX*X)/2)
    return g 

#-------------------------------------------------------------------------------
#obtain white noise according to the Wiener process:
#dW = lambda dt: np.random.normal(loc = 0.0, scale = np.sqrt(dt))
dW = lambda dt: np.random.normal(loc = 0.0, scale = 1)

#-------------------------------------------------------------------------------
def run_simulation(vars,pars,dt,tmax):
    factor=1000
    count=0
    tc=0.
    while(tc<=tmax):
        if (not count%factor):
            print(tc,'\t',vars.get('X'))
        #print(tc, '\t', vars.get('X'))
        a=getEffectiveRate(vars.get('X'),pars.get('gX'),pars.get('kX'))
        b=getNoiseStrength(vars.get('X'),pars.get('gX'),pars.get('kX'))
        #update the system:
        vars['X']+=a*dt+b*np.sqrt(dt)*dW(dt)
        tc+=dt
        count+=1

#-------------------------------------------------------------------------------
if __name__=='__main__':
   num_sims = 1
   for i_sim in range(num_sims):
       (vars,pars,dt,tmax)=defineSystem()  
       run_simulation(vars,pars,dt,tmax)
       print('\n')
