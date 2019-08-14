#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:22:35 2019

@author: Sylvia_Madhow
"""
import numpy as np

def process_h5_filelist(file_list,fs = 192000,frame_ms = 100.0):
    frame_size = int(frame_ms*fs/1000.0)
    bioSound_list = []
    feature_vectors = []
    labels = []
    for file in all_h5:
        print(file)
        hf = h5py.File(file, 'r')
        label = np.array(hf.get('label'))
        signal = np.array(hf.get('signal'))
        num_points = len(signal)
        num_frames = int(np.floor(num_points/frame_size))
        for i in range(num_frames):
            print(i)
            labels.append(label)
            frame = signal[i*frame_size:(i+1)*frame_size]
            bioSound = ssnd.BioSound(frame,fs)      
            bioSound.ampenv()
            bioSound.spectrum()
            bioSound.spectroCalc()
            bioSound.fundest()
            bioSound.mpsCalc(window = 64)
            bioSound_list.append(bioSound)
            features = make_feature_vector(bioSound)
            feature_vectors.append(features)



def make_feature_vector(bioSound):
    # make a feature vector from the features in a biosound object by averaging
    # across the total time frame of the signal
    features = np.zeros(18)
# Temporal statistics
    if bioSound.meantime.size != 0:
        features[0] = bioSound.meantime
    if bioSound.stdtime.size != 0:
        features[1] = bioSound.stdtime
    if bioSound.skewtime.size != 0:
        features[2] = bioSound.skewtime
    if bioSound.kurtosistime.size != 0:
        features[3] = bioSound.kurtosistime
    if bioSound.entropytime.size != 0:
        features[4] = bioSound.entropytime
# Pitch statistics
    if bioSound.sal.size != 0:
        features[5] = bioSound.sal
    if bioSound.fund.size != 0:
        features[6] = bioSound.fund
    if bioSound.F1.size != 0:
        goodF1 = bioSound.F1[~np.isnan(bioSound.F1)]
        meanF1 = np.mean(goodF1)
        features[7] = meanF1
    if bioSound.F2.size != 0:
        goodF2 = bioSound.F2[~np.isnan(bioSound.F2)]
        meanF2 = np.mean(goodF2)
        features[8] = meanF2
    if bioSound.F3.size != 0:
        goodF3 = bioSound.F3[~np.isnan(bioSound.F3)]
        meanF3 = np.mean(goodF3)
        features[9] = meanF3
# Spectral statistics
    if bioSound.meanspect.size != 0:
        features[10] = bioSound.meanspect
    if bioSound.stdspect.size != 0:
        features[11] = bioSound.stdspect
    if bioSound.skewspect.size!= 0:
        features[12] = bioSound.skewspect
    if bioSound.kurtosisspect.size != 0:
        features[13] = bioSound.kurtosisspect
    if bioSound.entropyspect.size != 0:
        features[14] = bioSound.entropyspect
    if bioSound.q1.size != 0:
        features[15] = bioSound.q1
    if bioSound.q2.size != 0:
        features[16] = bioSound.q2
    if bioSound.q3.size != 0:
        features[17] = bioSound.q3
# Modulation statistics
    return features
        