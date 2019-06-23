#--- Package import
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import interpolate
import matplotlib.mlab as mlab
from matplotlib import rc
import matplotlib as mpl
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


def temp_all(readpath_temp):
	#--- Reading the temperature files saced from the multimeter data saved via LabView.

	# data_err = np.genfromtxt(os.path.join(writepath_temp,'err_mean.txt'),delimiter=',')
	data_temp = np.genfromtxt(readpath_temp, usecols=(1,2,3), delimiter=',', skip_header=22, skip_footer=0)
	# data_temp_4wire = np.genfromtxt(readpath_4wire, usecols=(1,2,3), delimiter=',', skip_header=23)


	'''Calculate temperature from resistence time series'''

	T_in = tempfinder(data_temp[:,1],'4wire')
	T_bridge = tempfinder(data_temp[:,2],'4wire')


	'''Plotting the time series for temperature data'''

	plt.figure(figsize=(15,15))
	# plt.plot(data_temp[:,0], signal.detrend(T_in,type='linear'), 'r', label='TC temp-4wire', linewidth=3)
	plt.plot(data_temp[:,0], T_in, 'b', label='TC-4 wire', linewidth=3)
	# plt.plot(data_temp2[:13600,1], 'b', label='Temp in', linewidth=3)
	# plt.plot(data_temp2[13940:,1], 'b', label='Temp in', linewidth=3)
	# plt.plot(data_temp[200:,0]-data_temp[0,0], signal.detrend(T_bridge[200:],type='linear')/10, 'b', label='Frame temp-4 wire', linewidth=3)
	# plt.plot(data2_time,detr2/-10000,'r')
	# plt.plot(data_temp[24000:,0]-data_temp[24000,0],signal.detrend(data_temp[24000:,2],type='linear'), 'b', label='frame temp', linewidth=3)
	# plt.plot(data_temp[24000:,0]-data_temp[24000,0],signal.detrend(data_temp[24000:,1],type='linear'),'g')
	plt.legend(fontsize=28)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.grid(b=True, which='major', linewidth = '1.5', linestyle='--',color='black')
	plt.grid(which='minor')
	plt.xlabel('Time (sec)', fontsize=28)
	plt.ylabel('Temperature (C)', fontsize=28)
	plt.title('Tank Temperature', fontsize=36)
	# plt.savefig(os.path.join(writepath_4wire,'tc_temp.png'))

	plt.figure(figsize=(15,15))
	# plt.plot(data_temp[:,0], p1, 'r', label='Temp 1', linewidth=3)
	plt.plot(data_temp[:,0], T_bridge, 'b', label='frame temp', linewidth=3)
	# plt.plot(data_temp[:,0], T_bridge2, 'b', label='Frame temp-4 wire', linewidth=3)
	# plt.plot(data_temp[:,0], signal.detrend(T_bridge,type='linear'), 'b', label='Frame temp-4 wire', linewidth=3)
	plt.legend(fontsize=28)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.grid(b=True, which='major', linewidth = '1.5', linestyle='--',color='black')
	plt.grid(which='minor')
	plt.xlabel('Time (sec)', fontsize=28)
	plt.ylabel('Temperature (C)', fontsize=28)
	plt.title('TC Temperature', fontsize=36)
	# plt.savefig(os.path.join(writepath_4wire,'frame_temp.png'))

	plt.show()


	fs = 0.2
	f_4w, pxx_4w = psdmaker(T_in, 1, 2, 4, 0.2)
	f_b, p_b = psdmaker(T_bridge, 1, 2, 4, 0.2)
	# f_b, p_b = signal.welch(T_in3, fs, nperseg = 13600, detrend = 'linear', return_onesided=True)

	cte = 2e-8
	cav_length = 0.25
	scaling = cav_length*cte
	string  = r'$\Delta l$ = '+ str(scaling) + r'$\Delta T$'

	fig, ax1 = plt.subplots(figsize=(15,10))
	# ax2 = ax1.twinx()

	# ax1.figsize=(20,15)
	# ax1.loglog(f_4wire[0:], np.sqrt((p_4wire[0:]+p_b[0:]+p_r[0:])/3),'r', linewidth=2.5, label='TC temp noise')
	# ax1.loglog(f_b[0:], np.sqrt(p_b[0:]),'g', linewidth=2.5, label='Frame temp noise')
	ax1.loglog(f_r[0:], np.sqrt(p_r[0:]),'b', linewidth=2.5, label='Frtemp noise')
	ax1.loglog(f_4w[0:], np.sqrt(p_4w[0:]),'g', linewidth=2.5, label='Frtemp noise')
	# ax1.loglog(f, np.sqrt(pxx)*,'r', label=label, linewidth='2.5')
	# ax1.loglog(f_b[0:], lisa_req(f_b[0:])/(cav_l*cte),'o', linewidth=2.5, linestyle='--', label='LISA requirement')
	ax1.set_xlabel(r'Frequency [$Hz$]',fontsize=20)
	ax1.set_ylabel(r'Temperature stability [$K / \sqrt{Hz}$]',fontsize=20)
	ax1.grid(b=True, which='major', linewidth = '1.5',color='black')
	ax1.grid(which='minor')
	ax1.set_xlim(1e-5,1e-1)
	ax1.tick_params(axis="both", labelsize=20)
	ymin,ymax = ax1.get_ylim()

	# plt.loglog(f, np.sqrt(pxx)*4.25e-16,color='orange', label=label, linewidth=2.5)
	# ax2.loglog(f_4wire[0:], lisa_req(f_4wire[0:]),color='blue', linewidth=2.5, linestyle='--', label='LISA requirement')
	# ax2.set_ylim(ymin*(cav_l*cte),ymax*(cav_l*cte))
	# ax2.set_xlim(1e-5,1e-1)
	# ax2.set_ylabel(r'Length Noise [m/ \sqrt{Hz}]',fontsize=20)
	# ax2.tick_params(axis="both", labelsize=20)
	# ax2.text(1.5e-2,1e-11,string,fontsize=20)
	#plt.yticks([1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8], fontsize=28)

	# h1, l1 = ax1.get_legend_handles_labels()
	# h2, l2 = ax2.get_legend_handles_labels()
	# ax1.legend(h1+h2, l1+l2, fontsize=20)
	plt.title('Temperature Noise', fontsize=28)
	# time = data_temp[800:,0]
	ax1.legend(fontsize=20)
	# plt.savefig(os.path.join(writepath_4wire,'temp_noise_'+fileId+'.png'))
	plt.show()

	# return f_4wire, f_r, p_4w, p_r
	return
