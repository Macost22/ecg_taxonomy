# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 20:35:04 2022

@author: Melissa
"""
import os
import json
from visualization_ecg import plot_ecg_fiducial_points2
import numpy as np


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
    P = [np.nan if x == "NA" else x for x in P]
    for i in range(len(P)):
        try:
            amplitud_p = ECG[P[i]]
            amplitud_P.append(amplitud_p)
        except:
            continue
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
    
    # La amplitud de la onda P debe ester entre 0.15 y 0.2 mV
    amplitudP = np.mean(amplitud_P(paciente,fs))
    
    # La duración del complejo QRS debe estar entre 80 y 120 ms
    duracionQRS = np.mean(duracion_QRS(paciente,fs))
    
    # La amplitud de la onda T debe ser positiva
    amplitudT = np.mean(amplitud_T(paciente,fs))
    
    # El segmento PR debe durar menos de 200 ms 
    duracionPR = np.mean(duracion_PR(paciente, fs))
    
    amplitudP1, amplitudP2 = amplitud_P1_P2(paciente,fs)
    amplitudP1 = np.mean(amplitudP1)
    amplitudP2 = np.mean(amplitudP2)
    
    # HRmean esta normal entre 60 y 100 ms, RR dura entre 600 y 1200 ms
    HRmean, RR = HR_mean(paciente,fs)
    
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
             


if __name__ == '__main__':
    
    # Encontrar directoria de trabajo 
    # Path fiducial_point.json
    path_fiducial_True=path_fiducial='C:/Users/melis/Desktop/Bioseñales/ECG-Analysis/fiducial_points_true.json'


    # Se carga el archivo fiducial_points.json
    with open(path_fiducial_True) as f:
        fiducial_true = json.load(f)
        
    fs=250
    for persona in range(len(fiducial_true)):
        paciente = fiducial_true[persona]
        
        print('Paciente {}'.format(persona))
        taxonomy(paciente, fs)
        
        #segundos=2
        # Lista de personas
        #Lista_id=[2]
        #plot_ecg_fiducial_points2(fiducial_true,Lista_id,segundos=segundos,fs=250)
    
    
    