# Copyright 2016 Raytheon BBN Technologies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

def rabif(xdata,p0,p1,p2,p3):
    # model A + B * cos(w t + phi)
    return p0 - p1*np.cos(2*np.pi*p2*(xdata - p3))


def analyze_rabi_amp(data):
    
    
    # TODO - if we are trying to fit a curve to f(x) given x, why are we only
    # passing in the f(x)?  If data is not a function of xpts how will this 
    # work.  For now, the xpts vector matches below in the unit test
    numsteps = 40; #should be even
    stepsize = 2/numsteps
    
    xpts = np.arange(-(numsteps/2)*stepsize,-stepsize+stepsize,stepsize)
    xpts = np.append(xpts,np.arange(stepsize,(numsteps/2)*stepsize+stepsize,
    stepsize))

    # use largest FFT frequency component to seed Rabi frequency
    yfft = np.fft.fft(data)
    freqpos = np.argmax(np.abs( yfft[1:np.floor((len(yfft)-1)/2)] ))
    print(freqpos)
    frabi = 0.5*np.maximum(freqpos,1)/xpts[len(xpts)-1];
    print(frabi)
    
    # initial guess for amplitude is max - min
    amp = 0.5*(np.max(data) - np.min(data));
    offset = np.mean(data)
    phase = 0
    
    # check sign of amp
    if data[np.floor((len(yfft)-1)/2)] > offset:
        amp = -amp
    
    beta,_ = curve_fit(rabif,xpts,data,method='lm')
    print('Beta: ',beta)

    #The frequency tells us something about what a pi should calibrate to
    frabi = np.abs(beta[2])
    piAmp = 0.5/frabi;
    
    # The phase tell us somethign about the offset
    offsetPhase = beta[3];
    
    return beta
    
