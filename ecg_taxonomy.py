# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 20:35:04 2022

@author: Melissa
"""
import os
import json
from visualization_ecg import plot_ecg_fiducial_points2
import numpy as np
# Encontrar directoria de trabajo 

# Path fiducial_point.json
path_fiducial_True=path_fiducial='C:/Users/melis/Desktop/Bioseñales/ECG-Analysis/fiducial_points_true.json'


# Se carga el archivo fiducial_points.json
with open(path_fiducial_True) as f:
    fiducial_true = json.load(f)

segundos=2

# Lista de personas
Lista_id=[2]


plot_ecg_fiducial_points2(fiducial_true,Lista_id,segundos=segundos,fs=250)

#%%

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
    R=paciente['locs_Rav']
    R = [np.nan if x == "NA" else x for x in R]
    RR=[]
    HR=[]
    for ind in range(len(R)-1):
        RR.append(R[ind+1]/fs - R[ind]/fs)
        HR.append(1/(R[ind+1]/fs - R[ind]/fs)*60)
    HR_mean=round(np.mean(HR))
    return HR_mean, RR
#%%
fs=250
persona=fiducial_true[2]

duracionP = duracion_P(persona,fs)
amplitudP = amplitud_P(persona,fs)
duracionQRS = duracion_QRS(persona,fs)
amplitudT = amplitud_T(persona,fs)
duracionPR = duracion_PR(persona, fs)
amplitudP1P2 = amplitud_P1_P2(persona,fs)
HRmean, RR = HR_mean(persona,fs)
