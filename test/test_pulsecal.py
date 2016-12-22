# Copyright 2016 Raytheon BBN Technologies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0

import unittest
import numpy as np
#import matplotlib.pyplot as plt

import auspex.pulsecal.phase_estimation as pe
import auspex.pulsecal.optimize_amplitude as oa
import auspex.pulsecal.analyze_rabi as ar

def simulate_measurement(amp, target, numPulses):
    
    idealAmp = 0.34
    noiseScale = 0.05
    polarization = 0.99 # residual polarization after each pulse

    # data representing over/under rotation of pi/2 pulse
    # theta = pi/2 * (amp/idealAmp);
    theta = target * (amp/idealAmp)
    ks = [ 2**k for k in range(0,numPulses+1)]
    
    xdata = [ polarization**x * np.sin(x*theta) for x in ks];
    xdata = np.insert(xdata,0,-1.0)
    zdata = [ polarization**x * np.cos(x*theta) for x in ks];
    zdata = np.insert(zdata,0,1.0)
    data = np.array([zdata,xdata]).flatten('F')
    data = np.tile(data,(2,1)).flatten('F')
    
    # add noise
    #data += noiseScale * np.random.randn(len(data));
    vardata = noiseScale**2 * np.ones((len(data,)));
    
    return data, vardata
    

class PhaseEstimateTestCase(unittest.TestCase):

    def test_simulated_measurement(self):
        
        numPulses = 9
        amp = .55
        direction = 'X'
        target = np.pi
        
        # Using the same simulated data as matlab
        data, vardata =  simulate_measurement(amp, target, numPulses)   
        
        # Verify output matches what was previously seen by matlab
        phase, sigma = pe.phase_estimation(data, vardata, verbose=True)
        self.assertAlmostEqual(phase,-1.2012,places=4)
        self.assertAlmostEqual(sigma,0.0245,places=4)
        
class OptimizeAmplitudeTestCase(unittest.TestCase):

    def test_simulated_measurement(self):
        
        numPulses = 9
        amp = .55
        direction = 'X'
        target = np.pi
        
        # NOTE: this function is a place holder to simulate an AWG generating
        # a sequence and a digitizer receiving the sequence.  This function
        # is passed into the optimize_amplitude routine to be able to update
        # the amplitude as part of the optimization loop.
        def update_data(amp):
            data, vardata =  simulate_measurement(amp, target, numPulses)   
            return data,vardata
        
        # Verify output matches what was previously seen by matlab
        amp_out = oa.optimize_amplitude(amp, direction, target, update_data)
        
        # NOTE: expected result is from the same input fed to the matlab
        # routine
        self.assertAlmostEqual(amp_out,0.3400,places=4)

class AnalyzeRabiTestCase(unittest.TestCase):

    def test_analyze_rabi(self):
        
        numsteps = 40; #should be even
        stepsize = 2/numsteps
        
        xpts = np.arange(-(numsteps/2)*stepsize,-stepsize+stepsize,stepsize)
        xpts = np.append(xpts,np.arange(stepsize,(numsteps/2)*stepsize+stepsize,
        stepsize))
        
        amp = .2
        offset = 3.0
        freq = 1.
        phase = np.pi/4
        
        data = np.ndarray(xpts.shape,dtype=np.float)
        for i in range(0,len(xpts)):
            data[i]=ar.rabif(xpts[i],offset,amp,freq,phase)
            #data[i] += np.random.randn()*.01
                    
        piAmp,beta = ar.analyze_rabi_amp( data)
        
        
        fit = np.ndarray(xpts.shape,dtype=np.float)
        for i in range(0,len(xpts)):
            fit[i] = ar.rabif(xpts[i],beta[0],beta[1],beta[2],beta[3])
        
        self.assertAlmostEqual(beta[0],offset,places=3)
        self.assertAlmostEqual(beta[1],amp,places=3)
        self.assertAlmostEqual(beta[2],freq,places=3)
        self.assertAlmostEqual(beta[3],phase,places=3)
        
        '''
        plt.figure()
        plt.plot(xpts,data)
        plt.plot(xpts,fit,'r')
        plt.figure()
        plt.plot(xpts,fit-data)
        plt.show()
        '''
                
        
        
        #self.assertAlmostEqual(amp_out,0.3400,places=4)

if __name__ == '__main__':
    unittest.main()
