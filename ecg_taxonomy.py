# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 20:35:04 2022

@author: Melissa
"""
import os
import pandas as pd
import json
import numpy as np
from visualization_ecg import  plot_ecg_fiducial_points, plot_ecg_fiducial_points2, plot_original_ecg
from fiducial_point_detection import find_fiducial_points, find_R, butterworth_bandpass_filter, signal_average, normalization
import matplotlib.pyplot as plt


"""
Definición de funciones para realizar medidas en el ECG
"""
# Duración de la onda P
def duracion_P(paciente,fs):
    duracion_P = []
    
    P1=paciente['locs_P1']
    P2=paciente['locs_P2']
    
    P1 = [np.nan if x == "NA" else x for x in P1]
    P2 = [np.nan if x == "NA" else x for x in P2]

    for i in range(len(P1)):
        duracion_p = (P2[i]-P1[i])/fs
        duracion_P.append(duracion_p*1000)
    return duracion_P 

# Amplitud onda P
def amplitud_P(paciente,fs):
    ECG = paciente['ecg_average']
    amplitud_P = []
    P = paciente['locs_P']
    P2 = fiducial['locs_P2']
    P = [np.nan if x == "NA" else x for x in P]
    for i in range(len(P)):
        amplitud_p = ECG[P[i]]
        amplitud_p2 = ECG[P2[i]]
        amplitud_P.append(amplitud_p - amplitud_p2)
    return amplitud_P

# Duración complejo QRS
def duracion_QRS(paciente,fs):
    duracion_QRS = []
    Q = paciente['locs_Q']
    S = paciente['locs_S']
    
    Q = [np.nan if x == "NA" else x for x in Q]
    S = [np.nan if x == "NA" else x for x in S]
    for i in range(len(Q)):
        duracion_qrs = ((S[i]-Q[i])/fs)
        duracion_QRS.append(duracion_qrs*1000)
    return duracion_QRS

# Amplitud T
def amplitud_T(paciente,fs):
    amplitud_T = []
    ECG = paciente['ecg_average']
    T = paciente['locs_T']
    T = [np.nan if x == "NA" else x for x in T]
   
    for i in range(len(T)):
        try:
            amplitud_t = ECG[T[i]]
            amplitud_T.append(amplitud_t)
        except:
            continue
    return amplitud_T

# Bloqueo AV
# Cuando la duración del segmento PR < 200 ms
def duracion_PR(paciente,fs):
    duracion_PR=[]
    P1=paciente['locs_P1']
    R=paciente['locs_Rav']
    
    P1 = [np.nan if x == "NA" else x for x in P1]
    R = [np.nan if x == "NA" else x for x in R]
    for i in range(len(P1)):
        duracion_pr = (R[i]-P1[i])/fs
        duracion_PR.append(duracion_pr*1000)
    return duracion_PR

# Latido atrial prematuro
# Cuando amplitud de P1 y amplitud de P2 son diferentes

def amplitud_P1_P2(paciente,fs):
    amplitud_P1 = []
    amplitud_P2 = []
    ECG = paciente['ecg_average']
    P1 = paciente['locs_P1']
    P2 = paciente['locs_P2']
    
    P1 = [np.nan if x == "NA" else x for x in P1]
    P2 = [np.nan if x == "NA" else x for x in P2]
   
    for i in range(len(P1)):
        try:
            amplitud_p1 = ECG[P1[i]]
            amplitud_p2 = ECG[P2[i]]
            amplitud_P1.append(amplitud_p1)
            amplitud_P2.append(amplitud_p2)
        except:
            continue
    return amplitud_P1, amplitud_P2

# Calculo de la frecuencia cardíaca y duracion RR
def HR_mean(paciente,fs):  
    """
    Calculo de los intervalos RR, para determinar la frecuencia cardíaca promedio de cada ECG
    
    Parámetros:
    -----------
    R = paciente a analizar
    fs = int
        Frecuencia de muestreo
    Return
    -----------
    Frecuencia cardíaca media
    """
    R=paciente['locs_R']
    R = [np.nan if x == "NA" else x for x in R]
    RR=[]
    HR=[]
    for ind in range(len(R)-1):
        RR.append(R[ind+1]/fs - R[ind]/fs)
        HR.append(1/(R[ind+1]/fs - R[ind]/fs)*60)
    HR_mean=round(np.mean(HR))
    RR = list(map(lambda x: x * 1000, RR))
    RR =  np.round(RR,3)
    return HR_mean, RR
"""
TAXONOMY
"""

def taxonomy(paciente,fs):
    # La onda P debe durar menos de 120 ms
    duracionP = np.mean(duracion_P(paciente,fs))
    print('Duración onda P = {} ms'.format(round(duracionP,2)))
    
    # La amplitud de la onda P debe ester entre 0.15 y 0.2 mV
    amplitudP = np.mean(amplitud_P(paciente,fs))
    print('Amplitud onda P = {} mV'.format(round(amplitudP,2)))
    
    # La duración del complejo QRS debe estar entre 80 y 120 ms
    duracionQRS = np.mean(duracion_QRS(paciente,fs))
    print('Duración de QRS = {} ms'.format(round(duracionQRS,2)))
    
    # La amplitud de la onda T debe ser positiva
    amplitudT = np.mean(amplitud_T(paciente,fs))
    print('Amplitud onda T = {} mV'.format(round(amplitudT,2)))
    
    # El segmento PR debe durar menos de 200 ms 
    duracionPR = np.mean(duracion_PR(paciente, fs))
    print('Duración segmento PR = {} ms'.format(round(duracionPR,2)))
    
    amplitudP1, amplitudP2 = amplitud_P1_P2(paciente,fs)
    amplitudP1 = np.mean(amplitudP1)
    amplitudP2 = np.mean(amplitudP2)
    print('Amplitud P1 = {} y P2 = {} mV'.format(round(amplitudP1,2),round(amplitudP2,2)))
    
    # HRmean esta normal entre 60 y 100 ms, RR dura entre 600 y 1200 ms
    HRmean, RR = HR_mean(paciente,fs)
    print('Frecuencia cardíaca = {}'.format(round(HRmean,2)))
    print(HRmean)
    
    # El intervalo RR debe ser regular
    diffRR = np.diff(RR)
    
    if duracionPR > 200:
        print('Bloqueo AV \n')
    
    elif (amplitudP1 - amplitudP2) > 0.05:
        print('Latido atrial prematuro \n')
    
    elif duracionQRS > 120:
        print('Bloqueo de rama \n')
    
    elif HRmean < 60:
        print('Bradicardia \n')
    
    elif HRmean > 100:
        print('Taquicardia')
        if  duracionQRS < 120:
            print('Taquicardia supraventricular \n')
    else:
        print('sano')
             

def load_data_arrhythmia(file_path):
    ecg=pd.read_csv(file_path,sep=" ",index_col=0)
    # Se transponen los datos  (68,240000) = (individuo, observaciones)
    ecg=ecg.transpose()
    # Se modifican los índices para que sean de 0 a 67
    ecg.index = list(range(len(ecg)))
    #ecg1=ecg_70.iloc[1]
    return ecg

if __name__ == '__main__':
    
    # Estos tiene una fs= 250,  Wn_low = 60 y Wn_high = 0.5
    path_arritmia = 'C:/Users/melis/Desktop/Bioseñales/MIMIC/MIMIC_arritmia.txt'
    signals = load_data_arrhythmia(path_arritmia)
    signal = signals.iloc[0]
    
    fs=250
    Wn_low = 60
    Wn_high = 0.5
    # Parametros para detección de los puntos fiduciales respecto a picos R
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
    
    # Filtrado de la señal
    signal_filtered = butterworth_bandpass_filter(signal, Wn_low, Wn_high, fs, 3)
    
    # Normalización de la señal
    signal_normalized = normalization(signal_filtered)
    
    # Ubicación de los picos R en la señal
    locs_R = find_R(signal_normalized, height=0.8, distance=0.3*fs, fs=fs)
    
    # Promediado de la señal
    signal_av = signal_average(signal_normalized, locs_R, fs) 
    
    # Extracción de puntos fiduciales de la señal
    fiducial = find_fiducial_points(signal_av,fs,gr_r,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10)
    fiducial['locs_R'] = locs_R
    
    plot_original_ecg(signal_normalized,0,5,fs)
    plt.show()
    plot_ecg_fiducial_points(fiducial,0,5,fs)
    
    taxonomy(fiducial, fs)
    
#%%
#duracionP = duracion_P(fiducial,fs)
#duracionPmean = np.mean(duracion_P(fiducial,fs))


    



        

    
    
    