#--- Package import
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import interpolate
from matplotlib import rc
import os


import pyLISA_telescope as func_library


def temp_all(readpath_temp,writepath):
	#--- Reading the temperature files saced from the multimeter data saved via LabView.

	# data_err = np.genfromtxt(os.path.join(writepath_temp,'err_mean.txt'),delimiter=',')
	data_temp = np.genfromtxt(readpath_temp, usecols=(1,2,3), delimiter=',', 
		                      skip_header=22, skip_footer=0)
	# data_temp_4wire = np.genfromtxt(readpath_4wire, usecols=(1,2,3), delimiter=',', skip_header=23)


	'''Calculate temperature from resistence time series'''

	T_in = func_library.tempfinder(data_temp[:,1],'4wire')
	T_bridge = func_library.tempfinder(data_temp[:,2],'4wire')


	'''Plotting the time series for temperature data'''

	plt.figure(figsize=(15,15))
	# plt.plot((data_temp[:,0]-data_temp[0,0])/3600, signal.detrend(T_in,type='linear'), 'r', label='TC temp-4wire', linewidth=3)
	plt.plot((data_temp[:,0]-data_temp[0,0])/3600, T_in, 'b', 
		     label='TC-4 wire', linewidth=3)
	plt.legend(fontsize=28)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.grid(b=True, which='major', linewidth = '1.5', linestyle='--',color='black')
	plt.grid(which='minor')
	plt.xlabel('Time (Hours)', fontsize=28)
	plt.ylabel('Temperature (C)', fontsize=28)
	plt.title('Tank Temperature', fontsize=36)
	plt.savefig(os.path.join(writepath,'tc_temp.png'))

	plt.figure(figsize=(15,15))
	plt.plot((data_temp[:,0]-data_temp[0,0])/3600, T_bridge, 'b', 
		     label='frame temp', linewidth=3)
	# plt.plot((data_temp[:,0]-data_temp[0,0])/3600, signal.detrend(T_bridge,type='linear'), 'b', label='Frame temp-4 wire', linewidth=3)
	plt.legend(fontsize=28)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.grid(b=True, which='major', linewidth = '1.5', linestyle='--',color='black')
	plt.grid(which='minor')
	plt.xlabel('Time (Hours)', fontsize=28)
	plt.ylabel('Temperature (C)', fontsize=28)
	plt.title('TC Temperature', fontsize=36)
	plt.savefig(os.path.join(writepath,'frame_temp.png'))
	# plt.show()


	fs = 0.2
	f_4w, pxx_4w = func_library.psdmaker(T_in, 1, 2, 3, 0.2)
	f_b, p_b = func_library.psdmaker(T_bridge, 1, 2, 3, 0.2)

	cte = 2e-8
	cav_length = 0.25
	scaling = cav_length*cte
	string  = r'$\Delta l$ = '+ str(scaling) + r'$\Delta T$'

	fig, ax1 = plt.subplots(figsize=(15,10))
	ax2 = ax1.twinx()

	# ax1.figsize=(20,15)
	# ax1.loglog(f_b[0:], np.sqrt(p_b[0:]),'b', linewidth=2.5, label='Frtemp noise')
	ax1.loglog(f_4w[0:], np.sqrt(pxx_4w[0:]),'g', linewidth=2.5, label='Frtemp noise')
	ax1.loglog(f_4w[0:], np.sqrt(pxx_4w[0:]),'g', linewidth=2.5, linestyle='--')
	# ax1.loglog(f, np.sqrt(pxx)*,'r', label=label, linewidth='2.5')
	# ax1.loglog(f_4w[0:], lisa_req(f_4w[0:])/(cav_l*cte),'o', linewidth=2.5, 
	# 			linestyle='--', label='LISA requirement')
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
	ax1.legend(fontsize=20)
	plt.savefig(os.path.join(writepath,'temp_noise_'+fileId+'.png'))
	plt.show()

	# return f_4wire, f_r, p_4w, p_r
	return
