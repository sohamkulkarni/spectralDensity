#--- Package import
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
import os



#--- General use functions


def rinmaker(d,trend):
    '''
    Gives Relative Intensity of the time series, 
    'd' is the time series,
    'trend' is the option for detrend
    '''
    d = d/np.mean(d)
    if trend==True:
        d = signal.detrend(d,type='linear')
    return d

def lisa_req(f):
    '''
    Formulas directly from LISA proposal. Give the length and 
        acceleration noise requirement
    
    '''
    l_noise = (1e-12)*np.sqrt(1+(f/2e-3)**(-4))
    a_noise = (3e-15)*np.sqrt(1+(f/4e-4)**(-2))*np.sqrt(1+(f/8e-3)**(-4))
    return a_noise+l_noise

def find_nearest(arr_in, val):
    ''' 
    Finds the index closest to a given 'val' in the array 'arr_in' 
    Needed for 'psdmaker' function.
    
    '''
    array = np.asarray(arr_in)
    idx = np.argmin((np.abs(array - val)))
    return int(idx)

def psdmaker(ts_data, n1, n2, n3, fs):
    '''
    Calculates three PSDs and merges them by stitching together at different 
        frequencies. Number of averages change at 1mHz, 10mHz below.
    'ts_data' is whatever time series
    'n1, n2, n3' are the number of average of each PSD
    'fs' is sampling frequency
    
    '''
    start_psd = time.time()
    n_avg = n1
    f2, Pxx_den2 = signal.welch(ts_data, fs,nperseg = len(ts_data)/n_avg, 
                                detrend = 'linear', return_onesided=True) #-- PSD using the welch periodogram method with 16 averages, and linear detrend
    print("Avg n1 calculation done in:", time.time()-start_psd)
    start_psd = time.time()
    n_avg = n2
    f4, Pxx_den4 = signal.welch(ts_data, fs,nperseg = len(ts_data)/n_avg, 
                                detrend = 'linear', return_onesided=True)
    print("Avg n2 calculation done in:", time.time()-start_psd)
    start_psd = time.time()
    n_avg = n3
    f8, Pxx_den8 = signal.welch(ts_data, fs,nperseg = len(ts_data)/n_avg, 
                                detrend = 'linear', return_onesided=True)
    print("Avg n3 calculation done in:", time.time()-start_psd)
    start_psd = time.time()
    f = np.concatenate((f2[1:find_nearest(f2,1e-3)], 
                        f4[find_nearest(f4,1e-3):find_nearest(f4,1e-2)], 
                        f8[find_nearest(f8,1e-2):]))
    pxx = np.concatenate((Pxx_den2[1:find_nearest(f2,1e-3)], 
                          Pxx_den4[find_nearest(f4,1e-3):find_nearest(f4,1e-2)], 
                          Pxx_den8[find_nearest(f8,1e-2):]))
    print("Concatenation done in:", time.time()-start_psd)
    return f,pxx

def dCleaner(pm_data, trend):
    '''Cleaning data for phasemeter glitching. Also optional detrend of data. '''
    for i in range(len(pm_data)):
        if pm_data[i]>5e8: pm_data[i] = 10003.798e5 - pm_data[i]
    if trend==True: pm_data = signal.detrend(pm_data,type='linear')
    return pm_data



#--Functions for resistence calculation

def R(t,*args):
    i=args[0]
    j=args[1]
    p=[3.9083e-3, -5.775e-7, 0]
    r0 = 1000
    return (r0*(1+p[0]*t+p[1]*t**2+p[2]*t**3*(t-100))-data_temp[i,j])

def rfinder(Vo):
    '''When the entire bridge is used with an external voltage source'''
    Vs = -4
    r = 11e3*(Vs - (2*Vo))/((2*Vo) + Vs)
    return r

def tempfinder(Vo,method):
    '''
    Calculates temperature from the resistence time series
    'Vo' is the time series of the measurement.
    'method' chooses between 4wire and direct voltage
    
    All formulae and constants from the spec sheet of the thermistor.
    '''
    if method=='4wire' : r = rr(Vo)
    if method=='bridge' : r = rfinder(Vo)
    a1,b1,c1,d1 = 0.003354016, 0.000313649, 1.78017e-6, 3.07414e-7
    t_inv = a1 + b1*np.log10(r/10e3) +c1*(np.log10(r/10e3))**2 + d1*(np.log10(r/10e3))**3
    return (1/t_inv) -273.15

def rr(r):
    '''Equivalent resistence calculation for 4-wire calculation'''
#     return (r*33e3 - 242e6)/(11e3-r)
    return (33000*r)/(33000-r)
#     return (-11+((22000*r)/(22000-r)))
