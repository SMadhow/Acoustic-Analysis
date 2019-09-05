#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 13:14:11 2019

@author: Sylvia_Madhow
"""
import soundsig.sound as ssnd
import numpy as np
import matplotlib.pyplot as plt


def plot(biosound, DBNOISE=50, f_low=250, f_high=10000):
    # Plots a biosound in figures 1, 2, 3
    
        # Plotting Variables
        soundlen = np.size(biosound.sound)
        t = np.array(range(soundlen))
        t = t*(1000.0/biosound.samprate)

        # Plot the oscillogram + spectrogram
        plt.figure(1)
        plt.clf()
        # mngr = plt.get_current_fig_manager()
        # mngr.window.setGeometry(0, 260, 640, 545)
        
        
        # The oscillogram
        ax = plt.axes([0.45, 1.4, 0.855, 0.20])      
        ax.plot(t,biosound.sound, 'k')
        # plt.xlabel('Time (ms)')
        plt.xlim(0, t[-1])               
        # Plot the amplitude enveloppe  
        if biosound.tAmp.size != 0 :      
            ax.plot(biosound.tAmp*1000.0, biosound.amp, 'r', linewidth=2)
        ax.set_xticks([])
      
        # Plot the spectrogram
        ax = plt.axes([0.45, 0.75, 1.07, 0.6])
        cmap = plt.get_cmap('binary')
        
        if biosound.spectro.size != 0 :
            soundSpect = biosound.spectro
            if soundSpect.shape[0] == biosound.to.size:
                soundSpect = np.transpose(soundSpect)
            maxB = soundSpect.max()
            minB = maxB-DBNOISE
            soundSpect[soundSpect < minB] = minB
            minSpect = soundSpect.min()
            plt.imshow(soundSpect, extent = (biosound.to[0]*1000, biosound.to[-1]*1000, biosound.fo[0], biosound.fo[-1]), aspect='auto', interpolation='nearest', origin='lower', cmap=cmap, vmin=minSpect, vmax=maxB)
            plt.colorbar()
       
        plt.ylim(f_low, f_high)
        plt.xlim(0, t[-1])
        ax.set_xticks([])
        ax.set_yticks([])
        xlim = ax.get_xlim()

        
        # Power Spectrum
        
        ax = plt.axes([0.1,0.75,0.3,0.6])
        if biosound.psd.size != 0 :
            plt.plot(biosound.psd,biosound.fpsd, 'k-') 
            plt.ylabel('Frequency Hz')
            plt.xlabel('Power Linear')
            
        yl, yh, xl, xh = plt.axis()
        xl = 0.0
        xh = 10000.0
        plt.axis((yl, yh, xl, xh))
    
        if biosound.q1.size != 0:
            plt.plot([yl, yh], [biosound.q1, biosound.q1], 'k--')
            plt.plot([yl, yh], [biosound.q2, biosound.q2], 'k--')
            plt.plot([yl, yh], [biosound.q3, biosound.q3], 'k--') 
        
        if biosound.F1.size != 0:        
            F1Mean = biosound.F1[~np.isnan(biosound.F1)].mean()
            F2Mean = biosound.F2[~np.isnan(biosound.F2)].mean()
            F3Mean = biosound.F3[~np.isnan(biosound.F3)].mean()
            plt.plot([yl, yh], [F1Mean, F1Mean], 'r--', linewidth=2.0)
            plt.plot([yl, yh], [F2Mean, F2Mean], 'c--', linewidth=2.0)
            plt.plot([yl, yh], [F3Mean, F3Mean], 'b--', linewidth=2.0)
        ax.set_xlim(ax.get_xlim()[::-1])
        ax.set_ylim(xlim)
                     
    # Plot the fundamentals
        #plt.axes([0.1, 0.75, 0.855, 0.60])  
        if biosound.f0.size != 0 :
           # ax = plt.axes([0.1,0.1, 0.855, 0.6])
            ax = plt.axes([0.45,0.1, 0.855, 0.6])
            ax.plot(biosound.to*1000.0, biosound.f0, 'k', linewidth=3, label = 'fundamental')
            ax.plot(biosound.to*1000.0, biosound.f0_2, 'm', linewidth=3, label = 'fundamental 2')
            ax.plot(biosound.to*1000.0, biosound.F1, 'r--', linewidth=3, label = 'formant 1')
            ax.plot(biosound.to*1000.0, biosound.F2, 'c--', linewidth=3, label = 'formant 2')
            ax.plot(biosound.to*1000.0, biosound.F3, 'b--', linewidth=3, label = 'formant 3')
            plt.ylabel('Frequency (Hz)')
            plt.xlabel('Time (ms)')
            plt.legend()
            ax.set_xlim(xlim)
        plt.show()
            
           
    # Plot Power Spectrum
  #      plt.figure(3)
  #      plt.clf()
        # mngr = plt.get_current_fig_manager()
        # mngr.window.setGeometry(650, 260, 640, 545)

        # mngr = plt.get_current_fig_manager()
        # mngr.window.setGeometry(650, 260, 640, 545)

  
    # Table of results
        plt.figure(2)
        plt.clf()
        # mngr = plt.get_current_fig_manager()
        # mngr.window.setGeometry(320, 10, 640, 250)
        textstr = '%s  %s' % (biosound.emitter, biosound.type)
        plt.text(0.4, 1.0, textstr)
        if biosound.fund.size != 0:
            if biosound.fund2.size != 0:
                textstr = 'Mean Fund = %.2f Hz Mean Saliency = %.2f Mean Fund2 = %.2f PF2 = %.2f%%' % (biosound.fund, biosound.sal, biosound.fund2, biosound.voice2percent)
            else:
                textstr = 'Mean Fund = %.2f Hz Mean Saliency = %.2f No 2nd Voice Detected' % (biosound.fund, biosound.sal)
            plt.text(-0.1, 0.8, textstr)
            
        if biosound.fund.size != 0:
            textstr = '   Max Fund = %.2f Hz, Min Fund = %.2f Hz, CV = %.2f' % (biosound.maxfund, biosound.minfund, biosound.cvfund) 
            plt.text(-0.1, 0.7, textstr)
        textstr = 'Mean Spect = %.2f Hz, Std Spect= %.2f Hz' % (biosound.meanspect, biosound.stdspect)
        plt.text(-0.1, 0.6, textstr)
        textstr = '   Skew = %.2f, Kurtosis = %.2f Entropy=%.2f' % (biosound.skewspect, biosound.kurtosisspect, biosound.entropyspect)
        plt.text(-0.1, 0.5, textstr)
        textstr = '   Q1 F = %.2f Hz, Q2 F= %.2f Hz, Q3 F= %.2f Hz' % (biosound.q1, biosound.q2, biosound.q3 )
        plt.text(-0.1, 0.4, textstr)
        if biosound.F1.size != 0:
            textstr = '   For1 = %.2f Hz, For2 = %.2f Hz, For3= %.2f Hz' % (F1Mean, F2Mean, F3Mean )
            plt.text(-0.1, 0.3, textstr)
        textstr = 'Mean Time = %.2f s, Std Time= %.2f s' % (biosound.meantime, biosound.stdtime)
        plt.text(-0.1, 0.2, textstr)
        textstr = '   Skew = %.2f, Kurtosis = %.2f Entropy=%.2f' % (biosound.skewtime, biosound.kurtosistime, biosound.entropytime)
        plt.text(-0.1, 0.1, textstr)
        if biosound.rms.size != 0 and biosound.maxAmp.size != 0 :
            textstr = 'RMS = %.2f, Max Amp = %.2f' % (biosound.rms, biosound.maxAmp)
            plt.text(-0.1, 0.0, textstr)
        
        plt.axis('off')        
        plt.show()
        
    # Plot Modulation Power spectrum if it exists
    
        #ex = (spectral_freq.min(), spectral_freq.max(), temporal_freq.min(), temporal_freq.max())
        cmap = plt.get_cmap('binary')
        
        if biosound.mps.size != 0 :
        # Plot the modulation power spectrum
            plt.figure(3)
            ax = plt.axes([0.45,0.75,1.07,0.6])
            ex = (biosound.wt.min(), biosound.wt.max(), biosound.wf.min()*1e3, biosound.wf.max()*1e3)
            logMPS = 10.0*np.log10(biosound.mps)
            maxMPS = logMPS.max()
            minMPS = maxMPS-50
            logMPS[logMPS < minMPS] = minMPS
            plt.imshow(logMPS, interpolation='nearest', aspect='auto', origin='lower', cmap=cmap, extent=ex)
            plt.xlabel('Temporal Frequency (Hz)')
            plt.colorbar()
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.set_yticks([])
            plt.ylim((0,biosound.wf.max()*1e3))
            ylim = ax.get_ylim()
        
        # Plot the temporal modulation
            ax = plt.axes([0.45,1.4,0.855,0.20])
            temporal_dist = np.sum(10*np.log10(biosound.mps),axis=0)
            spectral_dist = np.sum(10*np.log10(biosound.mps),axis=1)
            spec_start = np.where(biosound.wf == 0)[0][0]
            spec_end = np.where(biosound.wf == biosound.wf.max())[0][0]
            plt.plot(biosound.wt,temporal_dist)
            ax.set_xlim(xlim)
            ax.set_xticks([])
        
        # Plot the spectral modulation
            ax = plt.axes([0.1,0.75,0.3,0.6])
            plt.plot(spectral_dist[spec_start:spec_end], 1e3*biosound.wf[spec_start:spec_end])
            #plt.plot(spectral_dist, 1e3*biosound.wf)
            plt.ylabel('Spectral Frequency (Cycles/KHz)')
            ax.set_ylim(ylim)
            ax.set_xlim(ax.get_xlim()[::-1])
        
        plt.pause(1)   # To flush the plots?
        plt.show()


def plot_mps(biosound):
    # Plot modulation power spectrum and flattened lines
    cmap = plt.get_cmap('binary')
    plt.figure(1)
    if biosound.mps.size != 0 :
        # Plot the modulation power spectrum
        ax = plt.axes([0.45,0.75,1.07,0.6])
        ex = (biosound.wt.min(), biosound.wt.max(), biosound.wf.min()*1e3, biosound.wf.max()*1e3)
        logMPS = 10.0*np.log10(biosound.mps)
        maxMPS = logMPS.max()
        minMPS = maxMPS-50
        logMPS[logMPS < minMPS] = minMPS
        plt.imshow(logMPS, interpolation='nearest', aspect='auto', origin='lower', cmap=cmap, extent=ex)
        plt.xlabel('Temporal Frequency (Hz)')
        plt.colorbar()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_yticks([])
        plt.ylim((0,biosound.wf.max()*1e3))
        ylim = ax.get_ylim()
        
        # Plot the temporal modulation
        ax = plt.axes([0.45,1.4,0.855,0.20])
        temporal_dist = np.sum(10*np.log10(biosound.mps),axis=0)
        spectral_dist = np.sum(10*np.log10(biosound.mps),axis=1)
        spec_start = np.where(biosound.wf == 0)[0][0]
        spec_end = np.where(biosound.wf == biosound.wf.max())[0][0]
        plt.plot(biosound.wt,temporal_dist)
        ax.set_xlim(xlim)
        ax.set_xticks([])
        
        # Plot the spectral modulation
        ax = plt.axes([0.1,0.75,0.3,0.6])
        plt.plot(spectral_dist[spec_start:spec_end], 1e3*biosound.wf[spec_start:spec_end])
        #plt.plot(spectral_dist, 1e3*biosound.wf)
        plt.ylabel('Spectral Frequency (Cycles/KHz)')
        ax.set_ylim(ylim)
        ax.set_xlim(ax.get_xlim()[::-1])
    plt.show()