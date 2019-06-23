#--- Package import
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import interpolate
import matplotlib.mlab as mlab
from matplotlib import rc
import matplotlib as mpl
import os


'''
User Inputs:

'date' when the measurment was started.
'readpath_temp' requires the file name to be specified. The file 
    should be in the correct date folder.
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
    for i in range(len(data2_freq)):
        if pm_data[i]>5e8: pm_data[i] = 10003.798e5 - pm_data[i]
    if trend==True: pm_data = signal.detrend(pm_data,type='linear')
    return pm_data

def transmission_maker(readpath_tran, writepath):
	'''
	Reading the files read from the NI Card. Read transmission, reflection and piezo data.
	''' 
	
	trans_data = np.genfromtxt(readpath_tran, usecols=(0,3,4,5), dtype=float,
	                           delimiter=',', skip_header = 23, skip_footer=0)
	refl_data = np.genfromtxt(readpath_tran, usecols=(0,6,7), dtype=float,
	                          delimiter=',', skip_header = 23, skip_footer=0)
	piezo_data = np.genfromtxt(readpath_tran, usecols=(0,1,2), dtype=float, 
		                       delimiter=',', skip_header = 23, skip_footer=0)
	
	# # src_data = np.genfromtxt(readpath_temp, usecols=(2), dtype=float, delimiter=',', skiprows=135000)
	# # time = np.genfromtxt(readpath_temp, usecols=(0), dtype=float, delimiter=',', skiprows=135000)
	# rfam_data = np.genfromtxt(readpath_temp, usecols=(0,3,4,5), dtype=float, delimiter=',', skip_header = 23, skip_footer=0)


	'''Save transmission data to variables.'''

	'''Multiply by 1000 to convert to mV'''
	time = trans_data[300000:,0]
	tc_trans = trans_data[300000:,1]*-1000
	src_trans = trans_data[300000:,2]*-1000
	lrc_trans = trans_data[300000:,3]*-1000
	tc_refl = refl_data[300000:,1]*-1000
	src_refl = refl_data[300000:,2]*-1000
	# lrc_refl = refl_data[:,3]*-1000
	m_piezo = piezo_data[300000:,1]*-1000
	g_piezo = piezo_data[300000:,2]*-1000


	'''Plotting the time series for transmission data'''

	plt.figure(figsize=(15,15))
	# plt.plot(data2_time,data2_freq, 'b', linewidth=2)
	# plt.plot(time, lrc_trans,'b.' , label='LRC transmission')
	# plt.plot(time[300000:], tc_trans[300000:],'b.' , label='TC transmission')
	# # plt.plot(data2_time,rinmaker(data2_freq,True)*500)
	plt.plot(time, m_piezo, 'b--',label='M piezo')
	plt.plot(time, g_piezo, 'r--', label='G piezo')
	# plt.plot(time[300000:], src_trans[300000:],'g.' , label='SRC transmission')
	# # plt.plot(time, rp,'r')
	plt.xticks(fontsize=20)
	plt.grid(b=True, which='major', color='black', linestyle='-', linewidth=1.2)
	plt.grid(b=True, which='minor', color='grey', linestyle='--', linewidth=0.7)
	plt.yticks(fontsize=20)
	plt.title('Time Series Data',fontsize=24)
	# plt.ylabel(r'Transmission Amplitude [$mV$]', fontsize=20)
	plt.ylabel(r'Beat Frequency [$Hz$]', fontsize=20)
	plt.xlabel(r'Time [$sec$]', fontsize=20)
	plt.legend(fontsize=20)
	# plt.savefig(os.path.join(writepath, 'transmission.png'))

	# plt.figure(figsize=(15,15))
	# plt.plot(data2_time,detr2, 'b', linewidth=2)
	# # plt.plot(data2_time[:-610000],data2_freq[:-610000])
	# # plt.plot(time, tc_trans,'b.' , label='TC transmission')
	# # plt.plot(time, src_trans,'g.' , label='SRC transmission')
	# # plt.plot(data2_time2,data2_freq2)
	# plt.plot(data2_time3,detr23)
	# # plt.plot(time, q, color ='orange',label='Q')
	# # plt.plot(time, i,label='I')
	# # plt.plot(time, rp,'r')
	# plt.grid(b=True, which='major', color='black', linestyle='-', linewidth=1.2)
	# plt.grid(b=True, which='minor', color='grey', linestyle='--', linewidth=0.7)
	# plt.xticks(fontsize=20)
	# plt.yticks(fontsize=20)
	# plt.title('Time Series Detrend Data',fontsize=24)
	# plt.ylabel(r'Frequency [$Hz$]', fontsize=20)
	# # plt.ylabel(r'Transmission Amplitude [$mV$]', fontsize=20)
	# plt.xlabel(r'Time [$sec$]', fontsize=20)
	# plt.legend(fontsize=20)
	# # # plt.savefig(os.path.join(writepath, 'transmission.png'))
	plt.show()


	'''PSD calculation for transmission/reflection and other signals'''

	# f, Pxx_den = signal.welch(rinmaker(detr2,True), 30.5,nperseg = len(detr2)/8, detrend = 'linear')
	# f_tc, Pxx_den_tc = signal.welch(rinmaker(tc_trans,True), 30.5,nperseg = len(tc_trans)/8, detrend = 'linear')
	# f_src, Pxx_den_src = signal.welch(rinmaker(src_trans,True), 30.5,nperseg = len(src_trans)/8, detrend = 'linear')
	# f_tc_2, Pxx_den_tc_2 = signal.welch(rinmaker(tc_refl,True), 30.5,nperseg = len(tc_refl)/8, detrend = 'linear')
	# f_src_2, Pxx_den_src_2 = signal.welch(rinmaker(src_refl,True), 30.5,nperseg = len(src_refl)/8, detrend = 'linear')
	f_tc, Pxx_den_tc = psdmaker(tc_trans,2,4,8, 30.5)
	f_src, Pxx_den_src = psdmaker(src_trans,2,4,8, 30.5)
	f_tc_2, Pxx_den_tc_2 = psdmaker(tc_refl,2,4,8, 30.5)
	f_src_2, Pxx_den_src_2 = psdmaker(src_refl,2,4,8, 30.5)
	# f_rin, Pxx_den_rin = signal.welch(rinmaker(rfam,True), 30.5,nperseg = len(q)/2, detrend = 'linear')
	# f_rfam, Pxx_den_rfam = signal.welch(rinmaker(dc_rfam,True), 30.5,nperseg = len(dc_rfam)/2, detrend = 'linear')
	# f_power, Pxx_den_power = signal.welch(rinmaker(power_laser,True), 30,nperseg = len(power_laser)/4, detrend = 'linear')
	# f_polar, Pxx_den_polar = signal.welch(rinmaker(polar_laser,True), 30,nperseg = len(power_laser)/4, detrend = 'linear')


	'''PSD plot for all noises in transmission/reflection etc. signals'''

	plt.figure(figsize=(20,12))
	# plt.loglog(f, np.sqrt(pxx)*4.25e-16,'r', label=label, linewidth='2.5')
	# plt.loglog(f_rin, np.sqrt(Pxx_den_rin),color='brown', label='rin', linewidth='2.5')
	# plt.loglog(f_rfam, np.sqrt(Pxx_den_rfam),color='orange', label='rfam dc', linewidth='2.5')
	plt.loglog(f_tc, np.sqrt(Pxx_den_tc),'b', label='TC trans noise',
			   alpha=0.7, linewidth='2.5')
	plt.loglog(f_tc_2, np.sqrt(Pxx_den_tc_2),'r', label='TC refl noise',
		       alpha=0.7, linewidth='2.5')
	# plt.loglog(f_power, np.sqrt(Pxx_den_power),color='black', label='Laser Power PD', linewidth='2.5')
	# plt.loglog(f_polar, np.sqrt(Pxx_den_polar),color='purple', label='Polarisation PD', linewidth='2.5')
	plt.loglog(f_src, np.sqrt(Pxx_den_src),'g', label='SRC trans noise',
		       alpha=0.7, linewidth='2.5')
	plt.loglog(f_src_2, np.sqrt(Pxx_den_src_2),'y', label='SRC refl noise',
		       alpha=0.7, linewidth='2.5')
	# plt.loglog(f[1:len(f)-1], lisa_req(f[1:len(f)-1]),'b', linewidth=2.5, linestyle='--', label='LISA requirement')
	#plt.ylim([1e1, 1e8])
	plt.xlim([1e-5, 2e1])
	plt.grid(b=True, which='major', linewidth = '1.5',color='black')
	plt.grid(which='minor')
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.legend(fontsize=28)
	# plt.title('LRC-SRC Noise Spectrum', fontsize=32)
	#plt.title('MOKU Noise Spectrum', fontsize=28)
	# plt.rc('text', usetex=True)
	#plt.rc('font', family='serif')
	plt.xlabel(r'Frequency [$Hz$]',fontsize=28)
	plt.savefig(os.path.join(writepath, 'transmission_noise.png'))
	plt.show()



	# return f_tc, f_src, Pxx_den_tc, Pxx_den_src
	return