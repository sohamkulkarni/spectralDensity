#--- Package import

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os
import time
import pyLISA_telescope as func_library


'''Main PSD computation routine'''

def spectral_den(readpath, writepath):
	'''
	This function will read the file from given path, detrend the data,
	plot time series, calculate and then plot the PSD. 
	The data required for all the plots can be returned to main.py 
	if needed for any further analysis.
	'''

	start_time = time.time()
	abs_start = time.time()
	'''Read stream files'''
	data2_time=np.genfromtxt(readpath,usecols=([0]), dtype=(float),delimiter=',',
	                         skip_header=0, skip_footer=0) #-- Column 0==time, column 1==frequency
	data2_freq=np.genfromtxt(readpath,usecols=([2]), dtype=(float),delimiter=',',
	                         skip_header=0, skip_footer=0)
	# data2_I = np.genfromtxt(os.path.join(genpath,date,readpath),usecols=([5]),
	#                         dtype=(float),delimiter=',', skip_header=0, skip_footer=0)
	# data2_Q = np.genfromtxt(os.path.join(genpath,date,readpath),usecols=([6]),
	                          # dtype=(float),delimiter=',', skip_header=0, skip_footer=0)
	print("Done reading files in:", time.time() - start_time)
	start_time = time.time()
	duration = str(round(data2_time[-1]/3600,2))
	label = duration+' hour '


	'''Prepare frequency noise data before any processing.'''
	detr2 = func_library.dCleaner(data2_freq,True)
	print("Done detrending files in:", time.time()-start_time)
	start_time = time.time()

	'''Plotting for time series'''
	plt.figure(figsize=(20,20))
	plt.plot(data2_time/3600,detr2,'r.', label=label ) #-- Linear detrended data
	plt.xlabel('Time (Hours)',fontsize=28)
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
	plt.plot(data2_time/3600,data2_freq,'r.', label=label) #-- Raw data
	plt.xlabel('Time (Hours)',fontsize=28)
	plt.ylabel('Beat Frequency (Hz)',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.grid(b=True, which='major', color='black', linestyle='-', linewidth=1.2)
	plt.grid(which='minor', color='grey', linestyle='--', linewidth=0.7)
	# plt.legend(fontsize=28)
	plt.title('Time Series Data',fontsize=28)
	# plt.show()
	plt.savefig(os.path.join(writepath,'time_series.png'))
	print('Plotting time series done in:', time.time()-start_time)
	start_time = time.time()


	'''Main PSD calculation'''
	f, pxx = func_library.psdmaker(detr2,1,2,4, 30.5)
	# f1,pxx1 = psdmaker(detr22,1,2,3)
	# f3,pxx3 = psdmaker(detr23,1,2,3)

	print('PSD calculation done in:', time.time()-start_time)
	start_time = time.time()


	'''PSD plot for frequency noise'''
	cav_length = 0.25
	frequency =  2.8195489e+14
	df_to_dl = cav_length/frequency

	plt.figure(figsize=(20,12))
	plt.loglog(f, np.sqrt(pxx),'r', label="SRC-LRC "+label,
	           linewidth='2.5')
	# plt.loglog(f[1:], func_library.lisa_req(f[1:]),'b', linewidth=2.5,
	# 		   linestyle='--', label='LISA requirement')
	# plt.ylim([5e-15, 4e-9])
	plt.xlim([1e-5, 1e0])
	plt.grid(b=True, which='major', linewidth = '1.5',color='black')
	plt.grid(which='minor')
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.legend(fontsize=28)
	plt.title('Cavity Length Noise Spectrum', fontsize=32)
	plt.xlabel(r'Frequency [$Hz$]',fontsize=28)
	# plt.ylabel(r'Length ASD [$m/ \sqrt{Hz}$]',fontsize=28)
	plt.ylabel(r'Frequency ASD [$Hz/ \sqrt{Hz}$]',fontsize=28)
	# plt.savefig(os.path.join(writepath,'length_noise.png'))
	plt.savefig(os.path.join(writepath,'frequency_noise.png'))
	plt.show()

	print('Plotting Spectrum done in:', time.time()-start_time)
	print('spectral_density done in:', time.time()-abs_start)

	

	# return data2_time, detr2, data2_freq, f, pxx 	'''If the spectral density data is needed for further analysis'''
	return 

