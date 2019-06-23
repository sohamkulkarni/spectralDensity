#--- Package import

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import interpolate
import matplotlib.mlab as mlab
from matplotlib import rc
import matplotlib as mpl
import os

#--- Primary file reading block
'''
'readpath' is the path variable for the files which are streamed.
'readpathSdCard' is the path variable for the files read on the iPad.
'writepath' to save results. All these variables are built here. 
The two 'genpath' variables are one for the folder where the data file
    is streamed. Here, my laptop and the lab computer.

User Inputs:

'date' when the measurment was started.
'fileId' is the timestamp on every data file. Find this from the file.
'skip_header' =8 if you are using an iPad file, else 0 unless you want
    to deliberately not consider a bunch of data at the start
'skip_footer' = 0 unless you want to deliberately not consider
    a bunch of data at the end
'''

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
    n_avg = n1
    f2, Pxx_den2 = signal.welch(ts_data, fs,nperseg = len(ts_data)/n_avg, 
                                detrend = 'linear', return_onesided=True) #-- PSD using the welch periodogram method with 16 averages, and linear detrend

    n_avg = n2
    f4, Pxx_den4 = signal.welch(ts_data, fs,nperseg = len(ts_data)/n_avg, 
                                detrend = 'linear', return_onesided=True)

    n_avg = n3
    f8, Pxx_den8 = signal.welch(ts_data, fs,nperseg = len(ts_data)/n_avg, 
                                detrend = 'linear', return_onesided=True)

    f = np.concatenate((f2[1:find_nearest(f2,1e-3)], 
                        f4[find_nearest(f4,1e-3):find_nearest(f4,1e-2)], 
                        f8[find_nearest(f8,1e-2):]))
    pxx = np.concatenate((Pxx_den2[1:find_nearest(f2,1e-3)], 
                          Pxx_den4[find_nearest(f4,1e-3):find_nearest(f4,1e-2)], 
                          Pxx_den8[find_nearest(f8,1e-2):]))
    
    return f,pxx

def dCleaner(pm_data, trend):
    '''Cleaning data for phasemeter glitching. Also optional detrend of data. '''
    for i in range(len(pm_data)):
        if pm_data[i]>5e8: pm_data[i] = 10003.798e5 - pm_data[i]
    if trend==True: pm_data = signal.detrend(pm_data,type='linear')
    return pm_data


'''Main PSD computation routine'''

def spectral_den(readpath, writepath):
	'''
	This function will read the file from given path, detrend the data,
	plot time series, calculate and then plot the PSD. 
	The data required for all the plots can be returned to main.py 
	if needed for any further analysis.
	'''


	'''Read stream files'''
	data2_time=np.genfromtxt(readpath,usecols=([0]), dtype=(float),delimiter=',',
	                         skip_header=0, skip_footer=0) #-- Column 0==time, column 1==frequency
	data2_freq=np.genfromtxt(readpath,usecols=([2]), dtype=(float),delimiter=',',
	                         skip_header=0, skip_footer=0)
	# data2_I = np.genfromtxt(os.path.join(genpath,date,readpath),usecols=([5]),
	#                         dtype=(float),delimiter=',', skip_header=0, skip_footer=0)
	# data2_Q = np.genfromtxt(os.path.join(genpath,date,readpath),usecols=([6]),
	                          # dtype=(float),delimiter=',', skip_header=0, skip_footer=0)
	duration = str(round(data2_time[-1]/3600,2))
	label = duration+' hour '


	'''Prepare frequency noise data before any processing.'''
	detr2 = dCleaner(data2_freq,True)


	'''Plotting for time series'''
	plt.figure(figsize=(20,20))
	plt.plot(data2_time,detr2,'r.', label=label ) #-- Linear detrended data
	plt.xlabel('Time (Sec)',fontsize=28)
	plt.ylabel('Beat Frequency (Hz)',fontsize=28)
	plt.xticks(fontsize=28, fontstretch = 50)
	plt.yticks(fontsize=28)
	plt.grid(b=True, which='major', color='black', linestyle='-', linewidth=1.2)
	plt.grid(which='minor', color='grey', linestyle='--', linewidth=0.7)
	# plt.legend(fontsize=28)
	plt.title('Time Series Detrend Data',fontsize=28)
	# plt.show()
	plt.savefig(os.path.join(writepath,'detrend.png'))

	plt.figure(figsize=(20,20))
	plt.plot(data2_time,data2_freq,'r.', label=label) #-- Raw data
	plt.xlabel('Time (Sec)',fontsize=28)
	plt.ylabel('Beat Frequency (Hz)',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	# plt.legend(fontsize=28)
	plt.title('Time Series Data',fontsize=28)
	# plt.show()
	plt.savefig(os.path.join(writepath,'time_series.png'))


	'''Main PSD calculation'''
	f, pxx = psdmaker(detr2,1,2,3, 30.5)
	# f1,pxx1 = psdmaker(detr22,1,2,3)
	# f3,pxx3 = psdmaker(detr23,1,2,3)


	'''PSD plot for frequency noise'''
	cav_length = 0.25
	frequency =  2.8195489e+14
	df_to_dl = cav_length/frequency

	plt.figure(figsize=(20,12))
	plt.loglog(f, np.sqrt(pxx)*df_to_dl,'r', label="TTS-LRC-new EOM "+label,
	           linewidth='2.5')
	plt.loglog(f[1:], lisa_req(f[1:]),'b', linewidth=2.5, linestyle='--', 
		       label='LISA requirement')
	# plt.ylim([5e-15, 4e-9])
	plt.xlim([1e-5, 1e0])
	plt.grid(b=True, which='major', linewidth = '1.5',color='black')
	plt.grid(which='minor')
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.legend(fontsize=28)
	plt.title('Cavity Stability Noise Spectrum', fontsize=32)
	#plt.title('MOKU Noise Spectrum', fontsize=28)
	plt.xlabel(r'Frequency [$Hz$]',fontsize=28)
	plt.ylabel(r'Frequency ASD [$m/ \sqrt{Hz}$]',fontsize=28)
	plt.savefig(os.path.join(writepath,'length_noise_tts.png'))
	plt.show()

	

	# return data2_time, detr2, data2_freq, f, pxx 	'''If the spectral density data is needed for further analysis'''
	return 

