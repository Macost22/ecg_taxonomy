# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:48:14 2022

@author: Melissa
"""
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, freqz
import matplotlib.pyplot as plt
import json
from scipy.signal import find_peaks



#path_fiducial='C:/Users/melis/Desktop/Bioseñales/ECG-Analysis/fiducial_points.json'
#with open(path_fiducial) as f:
#    fiducial = json.load(f)
    
"""
funciones para filtrar las señales de ecg con filtro Butterworth pasa bajas y altas
"""

def butterworth(cutoff, fs, order,btype):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype=btype, analog=False)
    return b, a

def butterworth_filter(data, cutoff, fs, order,btype):
    b, a = butterworth(cutoff, fs, order=order,btype=btype)
    y =filtfilt(b, a, data)
    return y

def butterworth_bandpass_filter(data,cutoff_low,cutoff_high,fs,order):
    b, a = butterworth(cutoff_low, fs, order,btype='low')
    y = butterworth_filter(data, cutoff_low, fs, order,btype='low')
    
    b1, a1 = butterworth(cutoff_high, fs, order,btype='high')
    data_filtered = butterworth_filter(y, cutoff_high, fs, order,btype='high')
    return data_filtered


"""
 averaging ecg signal
"""
def signal_average(signal, locs_R, fs):
    # longitud de cada segmento que se va a extraer (750 ms)
    segment=int(0.75*fs)
    # número de picos en la señal
    N=len(locs_R)
    # número de segmentos PQRST de la señal de ECG, que sera usados para promediar
    N1=10
    # el promedio inicia a partir del segundo beat (el primero puede no tratarse de un PQRST completo)
    N2=2
    
    PQRST = np.zeros((segment,N1))
    ECG= np.zeros((segment,N-N1-2))
    ECG1 = 0
    
    m=0
    # El k me permite identificar en que beat estoy ubicado para extraer el 
    #segmento de 1500 muestras
    k=1
    
    while k < (N-N1-1):
        for i in range(N1):
            start = int(locs_R[N2 + i - 1] - 0.25*fs)
            end = int(locs_R[N2 + i - 1] + 0.5*fs)
            PQRST[:,i] = signal[start:end] 
            ECG1 += PQRST[:,i] 
        
        ECG1= ECG1/N1
        ECG[:,m]=ECG1
        ECG1=0
        PQRST = np.zeros((segment,N1))
        N2 += 1
        m += 1  
        k += 1 
    
    ecg_average=ECG.flatten(order='F')
    return ecg_average

""" 
proceso de normalización 

"""
def normalization(signal):
    signal_max, signal_min = np.max(signal), np.min(signal)
    signal_normalized=(signal - signal_min)/(signal_max - signal_min)
    return signal_normalized

"""
Calculo de puntos fiduciales 
"""  
def find_R(signal,height, distance, mint): 
    locs_R, _ = find_peaks(signal, height=height, distance=distance)
    return locs_R + mint


def find_S(signal,locs_Rav,gr2):
    locs_S=[]
    for kk in range(len(locs_Rav)):
        start = locs_Rav[kk]
        end = int(start + gr2)
        ocs = np.argmin(signal[start:end])
        locs_S.append(ocs + locs_Rav[kk])
    return locs_S

def find_S2(signal,locs_S,gr10):
    locs_S2=[]
    for kk in range(len(locs_S)-1):
        start = locs_S[kk]
        end = int(start + gr10)
        ocs = np.argmax(signal[start:end])
        locs_S2.append(ocs + locs_S[kk])
    return locs_S2

def find_T(signal,locs_S,gr3):
    locs_T=[]
    for kk in range(len(locs_S)-1):
        start = locs_S[kk]
        end = int(start + gr3)
        ocs = np.argmax(signal[start:end])
        locs_T.append(ocs + locs_S[kk])
    return locs_T

def find_Q(signal, locs_Rav, gr4):
    locs_Q=[]
    for kk in range(len(locs_Rav)):
        start = locs_Rav[kk]
        end = int(start - gr4)
        ocs = np.argmin(signal[end:start])
        locs_Q.append(end + ocs)
    return locs_Q

def find_P(signal, locs_Q, gr5):
    locs_P=[]
    for kk in range(len(locs_Q)):
        start = locs_Q[kk]
        end = int(start - gr5)
        ocs = np.argmax(signal[end:start])
        locs_P.append(end + ocs)
    return locs_P
# Detección del inicio y final de las ondas P y T

def find_P1(signal, locs_P, gr6):
    locs_P1=[]
    for kk in range(len(locs_P)):
        start = locs_P[kk]
        end = int(start - gr6)
        ocs = np.argmin(signal[end:start])
        locs_P1.append(end + ocs)
    return locs_P1
    
def find_P2(signal, locs_P, gr7):
    locs_P2=[]
    for kk in range(len(locs_P)):
        start = locs_P[kk]
        end = int(start + gr7)
        ocs = np.argmin(signal[start:end])
        locs_P2.append(ocs + locs_P[kk])
    return locs_P2

def find_T1(signal, locs_T, gr8):
    locs_T1=[]
    for kk in range(len(locs_T)):
        start = locs_T[kk]
        end = int(start - gr8)
        ocs = np.argmin(signal[end:start])
        locs_T1.append(end + ocs)
    return locs_T1
    
def find_T2(signal, locs_T, gr9):
    locs_T2=[]
    for kk in range(len(locs_T)):
        start = locs_T[kk]
        end = int(start + gr9)
        ocs = np.argmin(signal[start:end])
        locs_T2.append(ocs + locs_T[kk])
    return locs_T2

    

def find_fiducial_points(signal,fs,gr_r,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10):
    fiducial_point=dict()
    
    n_valores = len(signal)
    tiempo = np.linspace(1,n_valores,n_valores)/fs
    mint = np.min(tiempo)*fs
    # Encontrar R 

    locs_Rav = find_R(signal, gr_r, 0.3*fs, 0)
    # Encontrar S
    locs_S = find_S(signal, locs_Rav, gr2)
    # Encontrar S2
    locs_S2 = find_S2(signal, locs_S, gr10)
    # Encontrar T
    locs_T = find_T(signal, locs_S, gr3)
    # Encontrar Q
    locs_Q = find_Q(signal, locs_Rav, gr4)
    # Encontrar P
    locs_P = find_P(signal, locs_Q, gr5)
    # Encontrar P1
    locs_P1 = find_P1(signal, locs_P, gr6)
    # Encontrar P2
    locs_P2 = find_P2(signal, locs_P, gr7)
    # Encontrar T1
    locs_T1 = find_T1(signal, locs_T, gr8)
    # Encontrar T2
    locs_T2 = find_T2(signal, locs_T, gr9) 
    
    fiducial_point={'locs_Rav':locs_Rav, 'locs_S':locs_S,'locs_S2':locs_S2, 'locs_T':locs_T, 'locs_Q':locs_Q, 
                    'locs_P':locs_P, 'locs_P1':locs_P1, 'locs_P2':locs_P2, 'locs_T1':locs_T1, 'locs_T2':locs_T2,'ecg_average':signal}

    return fiducial_point

    
#%%
path_fiducial='C:/Users/melis/Desktop/Bioseñales/ECG-Analysis/fiducial_points.json'
with open(path_fiducial) as f:
    fiducial_original = json.load(f)
    
#%%
path='C:/Users/melis/Desktop/Bioseñales/ECG_veronica/ecg_70.txt'
ecg_70=pd.read_csv(path,sep=" ")

# Se transponen los datos  (68,240000) = (individuo, observaciones)
ecg_70=ecg_70.transpose()
#%%
# Se modifican los índices para que sean de 0 a 67
ecg_70.index = list(range(len(ecg_70)))

ecg1=ecg_70.iloc[0]

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
locs_R = find_R(ecg_normalized, height=0.8, distance=0.3*fs,mint=0)


ecg_average=signal_average(ecg_normalized,locs_R,fs)
signal = ecg_average


# detección de picos S 

fiducial = find_fiducial_points(ecg_average,fs,gr_r,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10)

















#%%

from visualization_ecg import plot_ecg_fiducial_points, plot_ecg_fiducial_points2


plot_ecg_fiducial_points2(fiducial_original,[0],segundos=1,fs=2000)
plt.show()
plot_ecg_fiducial_points(fiducial,segundos=1,fs=2000)



