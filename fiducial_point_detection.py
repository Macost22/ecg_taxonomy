# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:48:14 2022

@author: Melissa
"""

import numpy as np
from scipy.signal import find_peaks
import pandas as pd
from scipy.signal import butter, filtfilt
from visualization_ecg import plot_ecg_fiducial_points, plot_ecg_fiducial_points2, plot_original_ecg
import json
import matplotlib.pyplot as plt
import neurokit2 as nk


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


def normalization(signal):
    signal_max, signal_min = np.max(signal), np.min(signal)
    signal_normalized=(signal - signal_min)/(signal_max - signal_min)
    return signal_normalized

def signal_average(signal, locs_R, fs):
    # longitud de cada segmento que se va a extraer (750 ms)
    segment=int(round(0.75*fs))
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
Calculo de puntos fiduciales 
"""  
def find_R(signal,height, distance, fs): 
    n_R = len(signal)*40/(fs*60)
    locs_R, _ = find_peaks(signal, height=height, distance=distance)
    while len(locs_R) < n_R:
        print("locs_R len: {} ".format(len(locs_R)))
        height -= 0.05
        print("nuevo gr_r: {}".format(height))
        locs_R, _ = find_peaks(signal, height=height, distance=distance)  
        print("nuevo locs_R len: {} \n".format(len(locs_R)))
     
        if height < 0.05:
            break       
    return locs_R 


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
    for kk in locs_P:
        start = kk
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

def find_fiducial_points_neurokit2(signal, gr_r, fs):

    n_valores = len(signal)
    stop = n_valores/fs
    tiempo = np.linspace(0,stop,n_valores)
    
    # Extract R-peaks locations
    _, locs_R = nk.ecg_peaks(signal, sampling_rate=fs)
    
    
    # Delineate
    if locs_R["ECG_R_Peaks"].size == 0:
        locs_R = find_R(signal, gr_r, 0.3*fs, fs)
        locs_R = {'ECG_R_Peaks':locs_R,'sampling_rate':fs}
        
    
    Signal, waves_peak = nk.ecg_delineate(signal, locs_R, sampling_rate=fs, method="dwt")

    waves_peak['ECG_R_Peaks'] = locs_R['ECG_R_Peaks']
    waves_peak['ECG'] = signal
    waves_peak['Tiempo'] = tiempo
    
    return waves_peak

def find_fiducial_points(signal,fs,gr_r,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10):
    fiducial_point=dict()
    
    n_valores = len(signal)
    stop = n_valores/fs
    #tiempo = np.linspace(1,n_valores,n_valores)/fs
    tiempo = np.linspace(0,stop,n_valores)

    # Encontrar R 
    _, locs_R = nk.ecg_peaks(signal, sampling_rate=fs)
    locs_R = locs_R["ECG_R_Peaks"]
    
    if locs_R.size == 0:
        locs_R = find_R(signal, gr_r, 0.3*fs, fs)
    # Encontrar S
    locs_S = find_S(signal, locs_R[1:], gr2)
    # Encontrar S2
    locs_S2 = find_S2(signal, locs_S, gr10)
    # Encontrar T
    locs_T = find_T(signal, locs_S, gr3)
    # Encontrar Q
    locs_Q = find_Q(signal, locs_R[1:], gr4)
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
    
    fiducial_point={'ECG_R_Peaks':locs_R[1:], 'ECG_S_Peaks':locs_S,'ECG_R_Offsets':locs_S2, 'ECG_T_Peaks':locs_T, 'ECG_Q_Peaks':locs_Q, 
                   'ECG_P_Peaks':locs_P, 'ECG_P_Onsets':locs_P1, 'ECG_P_Offsets':locs_P2, 'ECG_T_Onsets':locs_T1, 'ECG_T_Offsets':locs_T2,'ECG':signal,'Tiempo':tiempo}

    return fiducial_point

if __name__ == '__main__':
    #Se cargan los datos de ecg
    
    #path='C:/Users/melis/Desktop/Bioseñales/ECG_veronica/ecg_70.txt'
    path = 'C:/Users/melis/Desktop/Bioseñales/MIMIC/MIMIC_arritmia.txt'

    ecg = pd.read_csv(path,sep=" ",index_col=0)
    #ecg = pd.read_csv(path,sep=" ")
    #Se transponen los datos  (68,240000) = (individuo, observaciones)
    ecg = ecg.transpose()
    # Se modifican los índices para que sean de 0 a 67
    ecg.index = list(range(len(ecg)))
    #203 malo
    signal=ecg.iloc[7]
    fs=250
    Wn_low=60
    Wn_high=0.5

    
    gr_r = 0.8
    gr2 = 0.1 * fs # max number of samples for RS distance  
    gr3 = 0.32 * fs # max number of samples for ST distance = 0.32[s]*fs #640
    gr4 = 0.05 * fs # max QR distance
    gr5 = 0.2 * fs # max PQ distance
    gr6 = 0.08 *fs # max PP1 distance (P1 - beginning of P wave) 
    gr7 = 0.065 * fs  # max PP2 distance (P2 - end of P wave)
    gr8 = 0.08 * fs # max TT1 distance (T1 - beginning of T wave)
    gr9 = 0.08 * fs # max of TT2 distance (T2 - end of T wave)
    gr10 = 0.04 * fs # max of SS2 distance (S2 - end of QRS complex)
    
    t_start = 0
    t_end = 5
    plot_original_ecg(signal,t_start,t_end,fs)
    plt.show()
    
    # Filtrado de la señal
    signal_filtered = butterworth_bandpass_filter(signal, Wn_low, Wn_high, fs, 2)
    
    # Filtrado con neurokit
    signal_filtered_nk = nk.ecg_clean(signal, sampling_rate=fs, method="neurokit")

    
    # Normalización de la señal
    signal_normalized_nk = normalization(signal_filtered_nk)  
    signal_normalized = normalization(signal_filtered)  
   
    
    # Extracción de puntos fiduciales de la señal
    fiducial_neurokit =  find_fiducial_points_neurokit2(signal_normalized_nk, gr_r, fs)
    fiducial = find_fiducial_points(signal_normalized_nk,fs,gr_r,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10) 
    
    titulo1='Neurokit'
    titulo2='Algoritmo R'
    plot_ecg_fiducial_points(fiducial_neurokit,t_start,t_end,fs,titulo1)
    plt.show()
    plot_ecg_fiducial_points(fiducial,t_start,t_end,fs,titulo2)




    


