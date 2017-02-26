#/usr/bin/env python
''''
Simulation of a toggle switch system consisting of two mutually inhibiting 
transcription factors (TF) considering the existence of intrinsic noise in 
the birth-death process of each TF.
 
SDEs: intrinsic noise is present
dX/dt=gX*HS(Y)-kX*X+noise(X,Y,t) 
dY/dt=gY*HS(X)-kY*Y+noise(X,Y,t)

White noise is introduced in the system according to the Wiener process:
dX/dt=gX*HS(Y)-kX*X+g*sqrt(dt)*dW(dt)
dY/dt=gX*HS(X)-kY*Y+g*sqrt(dt)*dW(dt)
where: 
    HS models the mutual inhibition of the two genes X and Y
    W models the standard Brownian motion 
    g models the noise strength
'''
import sys
import numpy as np
from collections import OrderedDict
import auxiliary_functions as aux

#-----------------------------------------------------------------------------#
def defineSystem():
    #system variables:
    vars = OrderedDict()
    vars['X'] = 100
    vars['Y'] = 100
    
    #system parameters:
    pars = OrderedDict()
    pars = aux.parameter_st_1()

    #time step size:
    dt = 1.0e-3 

    #maximum simulation time:
    tmax = 1.0e+6
    return (vars, pars, dt, tmax) 

#-----------------------------------------------------------------------------#
def shifted_hill_impact(X, X0, nX, lamb):
    '''
    This method calculates and returns the effect of shifted Hill function.
    '''
    return lamb + (1.0-lamb)/(1.0+(X/X0)**nX)

#-----------------------------------------------------------------------------#
def getEffectiveRate(X, gX, kX, sh_impact_Y):
    return gX*sh_impact_Y - kX*X
    #return pars.get('gX') - pars.get('kX')*vars.get('X')
        
#-----------------------------------------------------------------------------#
def getNoiseStrength(X, gX, kX):
    #g = np.sqrt((gX+kX*X)/2)
    g = np.sqrt(abs((gX+kX*X))/2)
    return g 

#-----------------------------------------------------------------------------#
#obtain white noise according to the Wiener process:
#dW = lambda dt: np.random.normal(loc = 0.0, scale = np.sqrt(dt))
dW = lambda dt: np.random.normal(loc = 0.0, scale = 1)

#-----------------------------------------------------------------------------#
def run_simulation(vars, pars, dt, tmax):
    factor = 1000
    count = 0
    tc=0.
    while(tc<=tmax):
        if (not count%factor):
            print(tc, '\t', vars.get('X'), '\t', vars.get('Y'))
        #a = getEffectiveRate(vars, pars)
        #b = getNoiseStrength(vars, pars)
        #approximate the change in X-subsystem:
        sh_impact_Y = shifted_hill_impact(vars.get('Y'), pars.get('Y0'), 
                                          pars.get('nY'), pars.get('lY'))
        a = getEffectiveRate(vars.get('X'), pars.get('gX'), 
                             pars.get('kX'), sh_impact_Y)
        b = getNoiseStrength(vars.get('X'), pars.get('gX'), pars.get('kX'))
        dX = a * dt + b * np.sqrt(dt)*dW(dt)
        #update X-subsytem:
        vars['X'] += dX 

        #approximate the change in Y-subsystem:
        sh_impact_X = shifted_hill_impact(vars.get('X'), 
                                          pars.get('X0'), 
                                          pars.get('nX'), pars.get('lX'))
        a = getEffectiveRate(vars.get('Y'), pars.get('gY'),
                             pars.get('kY'), sh_impact_X)
        b = getNoiseStrength(vars.get('Y'), pars.get('gY'), 
                             pars.get('kY'))
        dY = a * dt + b * np.sqrt(dt)*dW(dt)
        #update Y-subsytem:
        vars['Y'] += dY 
        tc+=dt
        count+=1

#-----------------------------------------------------------------------------#
if __name__=='__main__':
   num_sims = 5
   for i_sim in range(num_sims):
       (vars, pars, dt, tmax) = defineSystem()  
       run_simulation(vars, pars, dt, tmax)
       print('\n')
