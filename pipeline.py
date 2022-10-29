# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 20:44:22 2022

@author: Melissa
"""
import pandas as pd
import numpy as np
import json
from visualization_ecg import plot_ecg, plot_original_ecg, plot_ecg_fiducial_points, plot_ecg_fiducial_points2
from fiducial_point_detection import butterworth_bandpass_filter, signal_average, normalization, find_fiducial_points, find_R


"""
Se cargan los datos de ecg
"""
path='C:/Users/melis/Desktop/Bioseñales/ECG_veronica/ecg_70.txt'
ecg_70=pd.read_csv(path,sep=" ")

# Se transponen los datos  (68,240000) = (individuo, observaciones)
ecg_70=ecg_70.transpose()

# Se modifican los índices para que sean de 0 a 67
ecg_70.index = list(range(len(ecg_70)))

ecg1=ecg_70.iloc[1]


t=np.linspace(0,120,240000)
# Setting standard filter requirements.
order = 3
fs = 2000       
cutoff_low = 100
cutoff_high=1

ecg_filtered=butterworth_bandpass_filter(ecg1, cutoff_low, cutoff_high, fs, order)



gr_r = 0.8
gr2 = 0.1 * fs # max number of samples for RS distance  
gr3 = 0.32 * fs # max number of samples for ST distance = 0.32[s]*fs #640
gr4 = 0.05 * fs # max QR distance
gr5 = 0.2 * fs # max PQ distance
gr6 = 0.08 *fs # max PP1 distance (P1 - beginning of P wave) 
gr7 = 0.065 * fs  # max PP2 distance (P2 - end of P wave)
gr8 = 0.1 * fs # max TT1 distance (T1 - beginning of T wave)
gr9 = 0.1 * fs # max of TT2 distance (T2 - end of T wave)
gr10 = 0.04 * fs # max of SS2 distance (S2 - end of QRS complex)



ecg_normalized=normalization(ecg_filtered)
locs_R = find_R(ecg_normalized, height=0.8, distance=0.3*fs)


ecg_average=signal_average(ecg_normalized,locs_R,fs)

# detección de picos S 

fiducial = find_fiducial_points(ecg_average,fs,gr_r,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10)
#%%



path_fiducial='C:/Users/melis/Desktop/Bioseñales/ECG-Analysis/fiducial_points.json'
with open(path_fiducial) as f:
    fiducial_pooints = json.load(f)
segundos=2

# Lista de personas
Lista_id=[1]

plot_ecg_fiducial_points2(fiducial_pooints,Lista_id,segundos=segundos,fs=2000)
plt.show()
plot_ecg_fiducial_points(fiducial,segundos=segundos,fs=2000)




