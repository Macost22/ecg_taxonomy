# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 20:44:22 2022

@author: Melissa
"""
import pandas as pd
import json
import matplotlib.pyplot as plt
from visualization_ecg import  plot_ecg_fiducial_points, plot_ecg_fiducial_points2
from fiducial_point_detection import find_fiducial_points, find_R, butterworth_bandpass_filter, signal_average, normalization
from ecg_taxonomy import taxonomy

fiducial =[]


def main(signal, fs, Wn_low, Wn_high):
    
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
    
    return fiducial        
    
    
def load_data_personality_traits(file_path):
    ecg_70=pd.read_csv(file_path,sep=" ")
    # Se transponen los datos  (68,240000) = (individuo, observaciones)
    ecg_70=ecg_70.transpose()
    # Se modifican los índices para que sean de 0 a 67
    ecg_70.index = list(range(len(ecg_70)))
    #ecg1=ecg_70.iloc[1]
    return ecg_70

def load_data_arrhythmia(file_path):
    ecg=pd.read_csv(file_path,sep=" ",index_col=0)
    # Se transponen los datos  (68,240000) = (individuo, observaciones)
    ecg=ecg.transpose()
    # Se modifican los índices para que sean de 0 a 67
    ecg.index = list(range(len(ecg)))
    #ecg1=ecg_70.iloc[1]
    return ecg
    



if __name__ == '__main__':
    
    # Path de los ecg a analizar
    # Estos tiene una fs= 2000,  Wn_low = 100 y Wn_high = 1
    #path='C:/Users/melis/Desktop/Bioseñales/ECG_veronica/ecg_70.txt'
    
    # Estos tiene una fs= 250,  Wn_low = 60 y Wn_high = 0.5
    path_arritmia = 'C:/Users/melis/Desktop/Bioseñales/MIMIC/MIMIC_arritmia.txt'
    signals = load_data_arrhythmia(path_arritmia)
    
    fs=250
    
    # Se analiza cada ecg dentro del ciclo for, para así obtener los puntos fiduciales
    for signal in range(len(signals)):
            print('paciente: {}'.format(signal))
            ecg = signals.iloc[signal]
            ecg_fiducial = main(signal = ecg, fs = 250, Wn_low = 60, Wn_high = 0.5)
            fiducial.append(ecg_fiducial)
            
    # for persona in range(len(fiducial)):
    #     paciente = fiducial[persona]
        
    #     print('Paciente {}'.format(persona))
    #     taxonomy(paciente, fs)

    

                
                