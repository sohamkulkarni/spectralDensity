#--- Package import

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import signal
# from scipy import interpolate
# import matplotlib.mlab as mlab
# from matplotlib import rc
# import matplotlib as mpl
import os


#--- Building all needed file paths


genpath = r'C:\Users\kulka\OneDrive - University of Florida\uf\research\PicometerStability\data'
# genpath = r'C:\Users\kulkarnisoham\OneDrive - University of Florida\uf\research\PicometerStability\data'
date = '2019_05_31'
fileId = '165920'
# duration = str(15)
# label = duration+' hour '+date
readpath = os.path.join(genpath, date, 'MokuPhasemeterStream_'+date+'_'+fileId+'.csv')
readpath_tran = os.path.join(genpath, date, '2019_05_29_tts_1.lvm')
readpath_temp = os.path.join(genpath, date, '2019_05_30_temp.lvm')
writepath = os.path.join(genpath,date,fileId)

if not os.path.isdir(writepath):
	print('Created new folder')
	os.mkdir(writepath)

import spectral_density

# timestamp, detrend_freq, ts_freq, f, pxx = spectral_density.spectral_den(readpath)
spectral_density.spectral_den(readpath) 

import transmission


transmission.transmission_maker(readpath_tran)

import temperature

temperature.temp_all(readpath_temp)

