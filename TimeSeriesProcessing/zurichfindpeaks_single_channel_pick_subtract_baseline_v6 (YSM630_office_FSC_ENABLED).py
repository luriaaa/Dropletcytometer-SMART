# -*- coding: utf-8 -*-
"""
ver 4.5
latest version of non-looped 

@author: Lu Ri

    changed FSC peak count gating to > median to > 80% percentile, increase the half window size to 100, pecentile to 85
    also, need to replace original FSC data by inverted, zero centered FSC data, otherwise the find_peaks doesnt know the negative peaks
    
    stopped using moving median to gate the threshold (doesnt make sense anymore, use median of entire frame instead)
    changed gating condition for droplet fluorescence filtering (from .any to .all)
    changed median thres from percentage to raw number (some sample the median could be wrong, use raw number better)
"""

import numpy as np
from scipy.signal import savgol_filter, find_peaks
from scipy.ndimage import rank_filter
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from math import sqrt
from numpy import argmax, argpartition, mean, median, zeros, std, percentile, where, concatenate
import time
import pandas as pd
#from pathlib import Path
import os
from flowio.create_fcs import create_fcs

class real_time_peak_detection():
    def __init__(self, array, sampleRate, duration, timeConst, startIdx,
                 thresfactor, fallfactor, minWidth):
        self.y = array
        self.length = len(array)
        self.threshold = thresfactor # 5 times std
        self.tc = timeConst
        self.alpha = -1/(sampleRate*np.log(timeConst))
        self.baseline = np.asarray(array[0]*np.ones(len(array))) # aka moving mean
        self.peakidx = zeros(len(array)) # 0 and 1 array
        self.peakdata = zeros(len(array)) # to hold peak data substracted baseline
        self.idx = startIdx # compute first moving mean from first x points
        self.minW = minWidth # minimal peak width as valid peak
        self.dur = sampleRate * duration # length limit for live thresholding 
        self.thresfall = fallfactor
        self.peakcount = 0
       
    def single_channel_findpeaks(self):
        # preprocess data, not applicable for realtime data
        # this is purely governed by the electronic noises, 
        # serves as a purpose to reduce electronic noise fluctuations, 
        # nothing to do with the droplet features (they are set by the time constant)
        data = savgol_filter(self.y, window_length=9, 
                                polyorder=2,deriv=0)
        
        # initial value matters
        detect = 0
        
        # sometimes the data wont center around 0, so need to compute first base
        # line and first mean
        if self.idx != 0:
            i = self.idx
            sigmaSqLast = std(data[:i])
            sigmaSqLast = sigmaSqLast * sigmaSqLast
            self.baseline[:i] = np.mean(self.baseline[:i])
        else:
            i = 0
            sigmaSqLast = 0
            
        while i < self.length:           
            # fetch mean from last value
            mean = self.baseline[i-1]
            next_val = data[i]
            alpha = self.alpha           
            
            # use exponetial weighted move mean
            diff = next_val - mean
            incr = alpha * diff       
            
            # update mean
            mean += incr 
            
            # update variance
            sigmaSq = (1-alpha)* (sigmaSqLast + diff*incr)
            sigmaSqLast = sigmaSq
            sigma = sqrt(sigmaSq)
            
            # update thresholds
            maxTmp = self.threshold * sigma
            thres = mean + maxTmp
            endthres = mean + self.thresfall * maxTmp
            
            # detection states
            if next_val > thres:
                detect = 1
                temp_idx = i + 1 # use a temp_idx, discard if not a valid peak
                while detect:
                    if temp_idx >= self.length:
                        i = temp_idx
                        break
                    elif data[temp_idx] > endthres: # modify if need fall factor
                        temp_idx += 1
                    else:
                        detect = 0
                        
                # after exit while loop, set those under peak point to baseline 
                self.baseline[i:temp_idx] = mean
                
                # only collect peak data if peak is wide enough
                width = temp_idx - i
                if width > self.minW:
                    self.peakidx[i:temp_idx] = 1
                    self.peakcount += 1
               
                i = temp_idx # move i ahead
            
            # if no peak, keep moving i ahead
            else:
                # dont forget to update baseline tracker
                self.baseline[i]=mean
                i += 1               

        # once loop over all data, update self.peakdata
        self.peakdata = (self.y-self.baseline)*self.peakidx
                
    # plot to see quality of peak identification
    def peak_data_plotter(self):    
        
        fig, axs = plt.subplots(4, sharex=True)
        x = np.arange(self.length)
        axs[0].plot(x, self.y)
        axs[1].plot(x, self.baseline)
        axs[2].plot(x, self.peakidx)
        axs[3].plot(x, self.peakdata)        
         
class peak_data_processer():
    

    def __init__(self, data, dropletmask, channels, 
                 max_width_coeff = 1.2, median_thres = 0.08, FSC_select = True,
                 FSC_chan = 5, FSC_pre_filt = False, FSC_thres = 300,
                 FSC_boundary_win = 59, FSC_boundary_short_shift = 7,
                 filehandle = None, trunc_start = 0, trunc_end = 0, pick_ch = 0,
                 droplet_cutoff = 0.5, col_names = None):
        
        # prepare fluorescence channel data
        self.filt_data = zeros(data.shape)
        if len(data.shape) > 1:
            for ch in channels:
                self.filt_data[:,ch] = savgol_filter(data[:,ch], window_length=9, 
                                polyorder=2,deriv=0)
        else:
            self.filt_data = savgol_filter(data, window_length=9, 
                                polyorder=2,deriv=0)        
        
        # prepare FSC data
        if FSC_pre_filt:
            self.FSC_data = savgol_filter(data[:,FSC_chan], window_length=5, 
                                polyorder=2,deriv=0)
        else:
            self.FSC_data = data[:,FSC_chan]
        
        # basic params
        self.thres_factor = median_thres
        self.FSC_thres = FSC_thres #prominence threshold in identifying the peaks
        self.channel = channels
        
        # Fetch droplet boundary parameters
        self.mask = dropletmask
        self.mwc = max_width_coeff
        self.win = FSC_boundary_win
        self.shift = FSC_boundary_short_shift
        self.FSC_selection = FSC_select
                
        # initialize 1-dimensional droplet attributes
        self.start_idxes = []
        self.droplet_widths = []
        self.FSC_start_idxes = []
        self.FSC_end_idxes = []
        self.FSC_peak_count = []
        self.abs_gating = []
        self.got_peak_gating = []
        self.FSC_gating = []
        self.have_cell = []
      
        # this is all-channels
        self.droplet_p2p = []
        self.droplet_AUC = []
        
        self.baseline_P2P = []
        self.baseline_AUC = []
        self.offset = []
        
        # for export, per 1 minute frame, passing this to save in buffer as metadata per min frame
        self.filehandle = filehandle
        self.trunc_start = trunc_start
        if trunc_end == 0:
            self.trunc_end = len(data)
        else:
            self.trunc_end = trunc_end
        self.pick_ch = pick_ch
        self.droplet_cutoff = droplet_cutoff
        self.total_extract_cell = 0
        self.positive_per = 0
        if col_names is None:
            self.col_names = [str(n) for n in channels]
        else:
            self.col_names = paired_name_list
                       
    def filter_droplets(self):
        # at this stage, only want start_idx and width for each selected droplet
        flag = 0
        width = 0
        
        for idx in range(len(self.mask)):
            if self.mask[idx] !=0 and flag ==0:
                flag = 1
                self.start_idxes.append(idx)
            
            if flag:
                if self.mask[idx] ==1:
                    width += 1
                else:
                    flag = 0
                    self.droplet_widths.append(width)
                    width = 0
                    
        self.start_idxes = np.asarray(self.start_idxes)
        self.droplet_widths = np.asarray(self.droplet_widths)
                
        median_width = np.median(self.droplet_widths)*self.mwc        
        droplet_to_delete = np.where(self.droplet_widths>median_width)
        self.start_idxes = np.delete(self.start_idxes, droplet_to_delete)
        self.droplet_widths = np.delete(self.droplet_widths, droplet_to_delete)
        print("Filtered droplet counts (size):", len(self.start_idxes))   
        
    def filter_droplet_by_empty_fluorescence(self):
        #threshold: use median of all; gating: delete <0.5*thres
        #to filter out abnormally low flurorescence empty droplets
        gating_thres = self.droplet_cutoff*np.median(self.droplet_p2p,axis=0) # [ch1,ch2....]
        
        
        boolmap = np.any(self.droplet_p2p,axis=1,where = self.droplet_p2p<gating_thres) #old version
       
        # boolmap=np.zeros(len(self.droplet_p2p),dtype=bool)
        # for i in range(len(boolmap)):
        #     boolmap[i]=np.all(self.droplet_p2p[i,:]<gating_thres)
       
        self.droplet_p2p=np.delete(self.droplet_p2p,boolmap,axis=0)
        self.droplet_AUC=np.delete(self.droplet_AUC,boolmap,axis=0)
        self.droplet_widths = np.delete(self.droplet_widths,boolmap)
        self.start_idxes = np.delete(self.start_idxes,boolmap)
        self.FSC_start_idxes=np.delete(self.FSC_start_idxes,boolmap)
        self.FSC_end_idxes=np.delete(self.FSC_end_idxes,boolmap)
        self.FSC_peak_count=np.delete(self.FSC_peak_count,boolmap)
        print("Filtered droplet counts (fl):", len(self.start_idxes))  
    
    def find_local_maxima(self, x, limit):
       x_dilate = rank_filter(x, -1, size=3)
       return (x_dilate == x) & (x>limit)
       
    def FSC_boundary_processor(self):        
        # shift FSC_data to median = 0
        FSC_median_shifted = self.FSC_data-median(self.FSC_data)
        # multiply with mask to extract "original" shape
       # FSC_extracted_orig = self.mask * FSC_median_shifted
        FSC_start_orig = np.array(self.start_idxes)
        FSC_end_orig = np.array(self.start_idxes+self.droplet_widths)
        total_count = len(FSC_start_orig)
        
        #extend 45 left and right
        win = self.win
        shift_1 = self.shift
        shift_2 = win-shift_1        
        FSC_extracted_corrected = -FSC_median_shifted
        thres = 2800
        
        # x = np.arange(len(self.FSC_data))
        # fig,ax = plt.subplots(2, sharex = True)
        # ax[0].plot(x,FSC_extracted_orig, color = 'k')
        # ax[0].scatter(self.start_idxes,np.ones(total_count)*1000, s = 10, c='r',marker='x')
        # ax[0].scatter((self.start_idxes+self.droplet_widths),np.ones(total_count)*1000, s = 10, c='b',marker='x')
        # # start = 1st >thres peak idx from left of FSC_start_orig
        # end = 1st > thres peak idx from right of FSC_end_orig
        # special treatment for left and right most one
        
        if (FSC_start_orig[0]-shift_1)<0:
            left = 0
        else:
            left = FSC_start_orig[0]-shift_1
            
        right = FSC_start_orig[0]+shift_2
        peaks = argmax(self.find_local_maxima(FSC_extracted_corrected[left:right],thres)>0)
        # if this returns an error, then increase the extender
        FSC_start_orig[0] = left + peaks
        
        for i in range(1, total_count):
            left = FSC_start_orig[i]-shift_1
            right = FSC_start_orig[i]+shift_2
            peaks = argmax(self.find_local_maxima(FSC_extracted_corrected[left:right],thres)>0)
            FSC_start_orig[i] = left + peaks
                        
        for i in range(total_count-1):
            left = FSC_end_orig[i]-shift_2
            right = FSC_end_orig[i]+shift_1
            peaks = argmax(self.find_local_maxima(FSC_extracted_corrected[left:right],thres)[::-1]>0)
            FSC_end_orig[i] = right + 1 - peaks 
                    
        if (FSC_end_orig[-1]+shift_1)>=len(FSC_extracted_corrected):
            right = len(FSC_extracted_corrected)
        else:
            right = FSC_end_orig[-1]+shift_1 
            
        left = FSC_end_orig[-1]-shift_2
        peaks= argmax(self.find_local_maxima(FSC_extracted_corrected[left:right],thres)[::-1]>0)
        FSC_end_orig[-1] = right + 1 - peaks 
        
        self.FSC_start_idxes = FSC_start_orig
        self.FSC_end_idxes = FSC_end_orig
        self.FSC_data = FSC_extracted_corrected
        
        # # visualize FSC plot       
        # ax[1].plot(x,FSC_extracted_corrected,color = 'k')
        # ax[1].scatter(FSC_start_orig,np.ones(total_count)*5000, s = 10, c='r',marker='x')
        # ax[1].scatter(FSC_end_orig,np.ones(total_count)*5000, s = 10, c='b',marker='x')           
    
    def FSC_cell_identifier(self):
        
        dropletcount=len(self.start_idxes)
        self.FSC_peak_count=zeros(dropletcount)

        for i in range(dropletcount):            
            peaks,_ = find_peaks(self.FSC_data[self.FSC_start_idxes[i]:self.FSC_end_idxes[i]],prominence = self.FSC_thres)
            self.FSC_peak_count[i]=len(peaks)

                    
    def get_p2p_and_AUC_value(self): 
        
        self.droplet_p2p = zeros([len(self.start_idxes),len(self.channel)])
        self.droplet_AUC = zeros([len(self.start_idxes),len(self.channel)])
        
        for i in range(len(self.start_idxes)):
            idx = self.start_idxes[i]
            width = self.droplet_widths[i]
            FSC_s_idx = self.FSC_start_idxes[i]
            FSC_e_idx = self.FSC_end_idxes[i]
            temp = np.asarray(self.filt_data[idx:(idx+width),:])
            temp_AUC = np.asarray(self.filt_data[FSC_s_idx:(FSC_e_idx+1),:])
            
            ch_idx = 0
            for ch in self.channel:             # here, self.channel is used be cause we can select only part of data   
                pts = 3
                temp_ch = temp[:,ch]
                # previously is minus min, but i think only take max value is enough
                self.droplet_p2p[i,ch_idx] = round(mean(temp_ch[argpartition(temp_ch, -pts)[-pts:]]))
                self.droplet_AUC[i,ch_idx] = round(mean(temp_AUC[:,ch]))
                ch_idx+=1 
                
    
    def rolling_median(self,data_1d,thres_factor,low_thres):
        median_of_all = median(data_1d)
        partialdata = data_1d[data_1d > 0.6*median_of_all]
        return max(((1+thres_factor)*median(partialdata)),low_thres)
    
    
    def empty_droplet_filter(self,plotting=True):
        ######### 1st layer ######## 
        # 1.08*moving median will be the adaptive threshold
        # if at least 1 channel above threshold (use .any()), then declair as True
        
        # verion 4 may have to use a hard frame median
        ######### 2nd layer ######## 
        # for the remaining "False" droplets, check if got ONLY 1 peak, in at least 1 channel
        # and that peak should be above certain prominence threshold (more than droplet baseline)
        ######### 3rd layer ########
        # check FSC data, see if there's something with no fluorescence signal but got FSC 
        
        # some preparation parameters
        #low_thres = 50 # 50 is approximately 1mV, if 0.08*median <50, then correct to this value
        total_count = len(self.droplet_p2p)
        channel_count = len(self.channel)
        thres_factor = np.array([0.08 for _ in range(channel_count)] )# default at 1.08*movemedian
        #movemed_win = 100
        FSC_half_win = 100

        self.abs_gating = zeros(total_count)
        self.got_peak_gating = zeros(total_count)
        self.FSC_gating = zeros(total_count)
        self.have_cell = zeros(total_count)
        

        # thres_factor can be ratio of median (for 1st sample) or raw value (thereafter)        
        if self.thres_factor is not None:
            thres_factor = np.array(self.thres_factor)
        if np.all(thres_factor<2):
            thres = np.around(np.median(self.droplet_p2p,axis=0)*(1+thres_factor))
        else:
            thres = thres_factor
            
        self.thres_factor=thres.tolist() #update thres_factor, including the 1st sample
        
        for i in range(total_count):
            # initiate boolean mask
            boolmask = zeros(channel_count)
            for ch in range(channel_count):
                if self.droplet_p2p[i,ch] >thres[ch]:
                    # assign got signal to any channel with significant signal
                    boolmask[ch]=1
                
            
            self.abs_gating[i] = boolmask.any()
            self.have_cell[i] = boolmask.any()
            
        log = sum(self.have_cell)
        print("thresholding method found ", log, " positive droplets")
           
        # peak method is to tidy up "loose-ends" of median thresholding methods, 
        # esp for channel 2 and channel 3, where the median level is high (insensitive to spikes)
        # for i in range(total_count):
        #     if self.have_cell[i]==0:
        #         boolmask = zeros(channel_count)
        #         partial_data = self.filt_data[self.start_idxes[i]:(self.start_idxes[i]+self.droplet_widths[i]),:]
        #         for ch in range(channel_count):
        #             peaks,_ = find_peaks(partial_data[:,ch],prominence=400)
        #             if len(peaks)>1:
        #                 boolmask[ch]=1
 
        #         self.got_peak_gating[i]=boolmask.any()
        #         self.have_cell[i]=boolmask.any()
                
        # print("peak method found ", sum(self.have_cell)-log, " positive droplets")
        # log=sum(self.have_cell)
   
        # shouldnt use median, because we are finding the "rare" cases
        # if using percentile, the window should be 50 or 100, easier for calculation
        perc = 85
        if self.FSC_selection:
            FSC_count_med=percentile(self.FSC_peak_count[:2*FSC_half_win],perc)
            
            # LEFT BOUNDARY 
            for i in range(FSC_half_win):
                if self.have_cell[i]==0 and self.FSC_peak_count[i]>FSC_count_med:
                    self.FSC_gating[i] = 1
                    self.have_cell[i]=1
            # middle part       
            for i in range(FSC_half_win, total_count-FSC_half_win):
                if self.have_cell[i]==0:
                    FSC_count_med = percentile(self.FSC_peak_count[(i-FSC_half_win):(i+FSC_half_win)],perc)
                    if self.FSC_peak_count[i]>FSC_count_med:
                        self.FSC_gating[i] = 1
                        self.have_cell[i]=1
            #right boundary            
            for i in range(total_count-FSC_half_win,total_count):
                # for the last points, dont need to update the percentile value anymore
                if self.have_cell[i]==0 and self.FSC_peak_count[i]>FSC_count_med:
                    self.FSC_gating[i] = 1
                    self.have_cell[i]=1
            
            print("FSC method found ", sum(self.have_cell)-log, " positive droplets")          
        
        
        # visualize the final result
        signals = np.zeros((len(self.filt_data),3)) # for gated droplet visualization
        for i in range(total_count):
            if self.abs_gating[i]==1:
                signals[self.start_idxes[i]:(self.start_idxes[i]+self.droplet_widths[i]),0]=1
            if self.got_peak_gating[i]==1:
                signals[self.start_idxes[i]:(self.start_idxes[i]+self.droplet_widths[i]),1]=1    
            if self.FSC_gating[i]==1:
                signals[self.start_idxes[i]:(self.start_idxes[i]+self.droplet_widths[i]),2]=1
        
        if plotting:
            x = np.arange(len(self.filt_data)) 
            channels = np.arange(channel_count+2)
            flg, ax = plt.subplots(channel_count+2, sharex = True) # +FSC and peakdetect_signal
            colors = ['c','g','y','r','brown']
            for ch in channels[:-2]:
                ax[ch].plot(x,self.fil_data[:,ch],color = colors[ch])
                ax[ch].set_ylim([0,25000])
            ax[channels[-2]].plot(x,self.FSC_data, color = 'k')
            ax[channels[-1]].plot(x,signals)
            
        self.total_extract_cell = sum(self.have_cell)
        self.positive_per = round(sum(self.have_cell)/total_count,3)
        print(100*self.positive_per, "% droplets has positive signal")
    
    
    def get_background(self):
        
        signals_p = np.array(self.droplet_p2p)
        signals_a = np.array(self.droplet_AUC)
        
        have_cell = np.array(self.have_cell)
        cell_val = signals_p[np.where(have_cell==1)[0],:]
        length = len(signals_p)
        baseline_p = np.zeros((cell_val.shape[0],cell_val.shape[1]))
        baseline_a = np.zeros((cell_val.shape[0],cell_val.shape[1]))
        half_win =7
        half_med_win = 3

        i=0
        for idx in where(have_cell==1)[0]:
            
            # left boundary
            if idx < half_win:
                base_idx_l = where(have_cell[0:idx]==0)
                base_idx_r = where(have_cell[(idx+1):(2*half_win)]==0)
                if base_idx_l[0].shape[0] == 0:
                    
                    # 1x5 channel medians
                    baseline_p[i,:] = median(signals_p[base_idx_r[0][:(2*half_med_win)],:],axis=0)
                    baseline_a[i,:] = median(signals_a[base_idx_r[0][:(2*half_med_win)],:],axis=0)
                    
                elif base_idx_l[0].shape[0]<half_med_win:
                    
                    baseline_p[i,:] = median(concatenate((signals_p[base_idx_l[0],:],
                                                   signals_p[(base_idx_r[0]+(idx+1))[:(2*half_win-len(base_idx_l[0]))],:])),axis=0)
                    baseline_a[i,:] = median(concatenate((signals_a[base_idx_l[0],:],
                                                   signals_a[(base_idx_r[0]+(idx+1))[:(2*half_win-len(base_idx_l[0]))],:])),axis=0)
                else:
                    baseline_p[i,:] = median(concatenate((signals_p[base_idx_l[0],:][-half_med_win:,:],
                                                   signals_p[(base_idx_r[0]+(idx+1))[:half_med_win],:])),axis=0)
                    baseline_a[i,:] = median(concatenate((signals_a[base_idx_l[0],:][-half_med_win:,:],
                                                   signals_a[(base_idx_r[0]+(idx+1))[:half_med_win],:])),axis=0)
                            
            elif idx<(length-half_win):
                base_idx_l = where(have_cell[(idx-half_win):idx]==0)
                base_idx_r = where(have_cell[(idx+1):(2*half_win)]==0)
                
                baseline_p[i,:] = median(concatenate((signals_p[((idx-half_win)+base_idx_l[0]),:][-half_med_win:,:],
                                               signals_p[(base_idx_r[0]+(idx+1))[:half_med_win],:])),axis=0)
                baseline_a[i,:] = median(concatenate((signals_a[((idx-half_win)+base_idx_l[0]),:][-half_med_win:,:],
                                               signals_a[(base_idx_r[0]+(idx+1))[:half_med_win],:])),axis=0)
                
            #right boundary
            else:
                base_idx_l = where(have_cell[(length-2*half_med_win):idx]==0)
                base_idx_r = where(have_cell[(idx+1):]==0)
                if base_idx_r[0].shape[0] == 0:
                    
                    # 1x5 channel medians
                    baseline_p[i,:] = median(signals_p[base_idx_l[0][-(2*half_med_win):],:],axis=0)
                    baseline_a[i,:] = median(signals_a[base_idx_l[0][-(2*half_med_win):],:],axis=0)
                    
                elif base_idx_r[0].shape[0]<half_med_win:
                    
                    baseline_p[i,:] = median(concatenate((signals_p[(base_idx_r[0]+idx+1),:],
                                                   signals_p[(base_idx_l[0]+(length-2*half_med_win))[:(2*half_win-len(base_idx_l[0]))],:])),axis=0)
                    baseline_a[i,:] = median(concatenate((signals_a[(base_idx_r[0]+idx+1),:],
                                                   signals_a[(base_idx_l[0]+(length-2*half_med_win))[:(2*half_win-len(base_idx_l[0]))],:])),axis=0)
                else:
                    baseline_p[i,:] = median(concatenate((signals_p[(base_idx_l[0]+length-2*half_med_win),:][-half_med_win:,:],
                                                   signals_p[(base_idx_r[0]+(idx+1))[:half_med_win],:])),axis=0)
                    baseline_a[i,:] = median(concatenate((signals_a[(base_idx_l[0]+length-2*half_med_win),:][-half_med_win:,:],
                                                   signals_a[(base_idx_r[0]+(idx+1))[:half_med_win],:])),axis=0)
            
            i = i+1
    
        #at last, get baseline to matrix
        self.baseline_P2P=baseline_p
        self.baseline_AUC=baseline_a
    
    def compile_data(self,cell_only_filter):
        
        ### do some modification, compile as list of tuple instead of np ###
        #a cell is a dict {p2p(pandas_dataframe),auc(),fsc(st,end,peakcount),baseline(p2p,auc),drop attr(start, witdth), condition (abs_detect,peak_detect_fsc_detect)}
        #return a list of dictionaries of "have_cell ==1" droplets
        
        # target name paired with fluorophore
        paired_name_list = self.col_names
        
        # clean up droplet values
        
        if cell_only_filter:
            idxes = self.have_cell==1
        else:
            idxes = np.ones(len(self.have_cell),dtype=bool)
            
        p2p = pd.DataFrame(self.droplet_p2p[idxes],columns = paired_name_list)
        auc = pd.DataFrame(self.droplet_AUC[idxes],columns = paired_name_list)
        fsc = pd.DataFrame(np.column_stack((self.FSC_start_idxes[idxes],self.FSC_end_idxes[idxes],self.FSC_peak_count[idxes])),columns = ["FSC_start","FSC_end","FSC_pkcnt"])
        bl_p = pd.DataFrame(self.baseline_P2P,columns = paired_name_list)
        bl_a = pd.DataFrame(self.baseline_AUC,columns = paired_name_list)
        drop_attr = pd.DataFrame(np.column_stack((self.start_idxes[idxes],self.droplet_widths[idxes])),columns = ["Drop_start","Drop_width"])
        cond = pd.DataFrame(np.column_stack((self.abs_gating[idxes],self.got_peak_gating[idxes],self.FSC_gating[idxes])),columns = ["ABS_gating","PEAK_gating","FSC_gating"])
        median_thres = pd.DataFrame(self.thres_factor).T
        median_thres = median_thres.rename(columns = lambda x: 'CH'+str(x))
        FSC_thres = pd.DataFrame([self.FSC_selection*self.FSC_thres],columns = ['FSC_thres'])
        cell_count = pd.DataFrame([self.total_extract_cell],columns = ['Cell_NO'])
        pos_per = pd.DataFrame({'POS_perc':[self.positive_per]})
        signal_name = pd.DataFrame({'Filename':[self.filehandle]})
        trunc_bounds = pd.DataFrame({'Data_start':[self.trunc_start],'Data_end':[self.trunc_end]})
        pick_channel = pd.DataFrame({'Pick_CH':[self.pick_ch]})
        cutoff = pd.DataFrame({'Pick_CH':[self.droplet_cutoff]})

        return {"P2P_data":p2p,
                "AUC_data":auc,
                "FSC_data":fsc,
                "P2P_baselines":bl_p,
                "AUC_baselines":bl_a,
                "Droplet_Attribute":drop_attr,
                "Gating_Condition":cond,
                "Median_threshold":median_thres,
                "FSC_threshold":FSC_thres,
                "Cell_NO": cell_count,
                "POS_percentage":pos_per,
                "Filename":signal_name,
                "Data_Boundary":trunc_bounds,
                "Pick_Channel":pick_channel,
                "Droplet_Cutoff":cutoff}


#### Step 0:Initializers, all the basic parameters for teh program
channel_list = [0,1,2,3,4]   
pick_ch = 2
# prepare names (column names)
channel_names = [("CH"+ str(ch) )for ch in channel_list]
fl_names = ["5FAM","MF570","PEDZ","PC55","PC7"]
target_names = ['NE','GzB','CD66b','CD3','CD31'] # fill accordingly wrt to experiment

paired_name_list = [dict(zip(channel_names, target_names))[key]+'_'+dict(zip(channel_names,fl_names))[key] for key in channel_names]

median_thres=[0.25,0.4,0.5,0.5,0.25]  #intial median_thres, fill in ratio, will be converted in to raw value
modifier = np.array([150,250,200,50,0])
FSC_thres = 1300
per_range = [0.01,0.05] # if too many cells, then limit the range of percentage
per = per_range[1]+0.01
FSC_sel =True

# memory holders
signals_buff = []  
signals = [0] *2
baselines = [0]*2

signalname = "WB_24H_PMA"
mother_dir = 'E:/CAMP Biochemical/20220210 C006/processed'
#mother_dir = 'P:/CAMP Staff/Lu Ri/20220420 CAR D14 4/processed'

### direct run methods ####
def load_and_process(channel_list, pick_ch, paired_name_list, modifier):
    global median_thres, FSC_thres, per_range, per, FSC_sel
    plt.close('all')
    root = tk.Tk()
    root.withdraw()
    file_name = filedialog.askopenfilename( filetypes =[('numpy', '*.npy')])
    print("Current File: ", os.path.basename(file_name))
    data = np.load(file_name)
    # to clean up "dirty starts"
    start = 0
    end = len(data)
    #end = 9028000
    while np.any(data[start:(start+60),pick_ch]>500):
        start=start+100
    data = data[start:end,:]
        
    # peakstart = mean + thresfactor * std, peakstop = mean+ fallfactor*std
    # time constant [s], if tc = 1 we will have a LP filter of 1Hz (-3db)
    detector = real_time_peak_detection(data[:,pick_ch], 
                                        sampleRate = 200000,
                                        duration = 60, timeConst = 0.01, 
                                        startIdx = 50,
                                        thresfactor = 4, fallfactor = 1, 
                                        minWidth = 80)
    # find peaks
    t = time.time()
    detector.single_channel_findpeaks()    
    print("Droplet counts:", detector.peakcount)
    print("Elapsed time:", int(time.time()-t))    
    # plot function, ensure you capture all droplets
    detector.peak_data_plotter() 
     

    processor = peak_data_processer(data = data, 
                                    dropletmask = detector.peakidx,
                                    channels = channel_list,
                                    max_width_coeff = 1.6, median_thres=median_thres, FSC_select = FSC_sel,
                                    FSC_chan = 5, FSC_pre_filt = False, FSC_thres = FSC_thres,
                                    FSC_boundary_win=69,FSC_boundary_short_shift=9,
                                    filehandle=os.path.basename(file_name).split('.')[0],
                                    trunc_start=start,trunc_end=end,pick_ch=pick_ch,
                                    droplet_cutoff = 0.6,
                                    col_names = paired_name_list)
    
    #### Step 3: find start_idx, width and delete large droplets ####
    processor.filter_droplets()        
    #### Step 3.5 process FSC data, still need opimize ####
    processor.FSC_boundary_processor()
    processor.FSC_cell_identifier()
        
    #### Step 4: find peak2peak of each droplet, check if need to adjust median_thres ####
    processor.get_p2p_and_AUC_value()  
    
    print(dict(zip(channel_names,np.median(processor.droplet_p2p,axis=0))))   
    processor.filter_droplet_by_empty_fluorescence()        
    #### Step 5: identify empty droplet ####
    processor.empty_droplet_filter(plotting = False)
    processor.get_background()
    
    median_thres = processor.thres_factor #update the thres_factor for next sample
    
    per = processor.positive_per
    #FSC_per = sum(processor.FSC_gating)/sum(processor.have_cell)

        
    per_range[0] = per-0.005
    per_range[1] = per + 0.005 # update per_range, every valid per (can exit while loop)
    per = per_range[1]+0.001
    return processor

def evaluate_result(median_thres, processor, final_plotting = False):
    if final_plotting:
        processor.empty_droplet_filter(plotting = True)
    flg, ax = plt.subplots(5, sharex = True)
    colors = ['c','g','y','r','brown']
    for i in range(5):
        ax[i].plot(np.arange(len(processor.droplet_p2p)),processor.droplet_p2p[:,i],color = colors[i])
        ax[i].plot(np.arange(len(processor.droplet_p2p)),
                   np.ones(len(processor.droplet_p2p))*(median_thres[i]),'k-')
        
        
        ax[i].plot(np.arange(len(processor.baseline_P2P)),processor.baseline_P2P[:,i],color = 'b')
        ax[i].set_ylim([0,median_thres[i]*2])
        
def export_data(signalname,mother_dir, signals, baselines, signals_buff):
    # compile metadata as well (threshold, FSC_threshold,baselines)
    # need to pass filename through processor because we only do this at last step
    metadata = pd.concat([pd.concat([result["Median_threshold"] for result in signals_buff], ignore_index=True),
                         pd.concat([result["FSC_threshold"] for result in signals_buff], ignore_index=True),
                         pd.concat([result["Cell_NO"] for result in signals_buff], ignore_index=True),
                         pd.concat([result["POS_percentage"] for result in signals_buff], ignore_index=True),
                         pd.concat([result["Filename"] for result in signals_buff], ignore_index=True),
                         pd.concat([result["Data_Boundary"] for result in signals_buff], ignore_index=True),
                         pd.concat([result["Pick_Channel"] for result in signals_buff], ignore_index=True),
                         pd.concat([result["Droplet_Cutoff"] for result in signals_buff], ignore_index=True)], axis=1)
    
    # at last, fetch either P2P or AUC value    
    signals[0]=pd.concat([result["P2P_data"] for result in signals_buff], ignore_index=True)# use "AUC_data" for AUC
    signals[1]=pd.concat([result["AUC_data"] for result in signals_buff], ignore_index=True)
    baselines[0] = pd.concat([result["P2P_baselines"] for result in signals_buff], ignore_index=True)
    baselines[1] = pd.concat([result["AUC_baselines"] for result in signals_buff], ignore_index=True)
    
    attributes = pd.concat([pd.concat([result["Droplet_Attribute"] for result in signals_buff], ignore_index=True),
                            pd.concat([result["FSC_data"] for result in signals_buff], ignore_index=True),
                            pd.concat([result["Gating_Condition"] for result in signals_buff], ignore_index=True) ],axis=1)
    
   
    #### Step 8: name your file and save everything
       
    os.chdir(mother_dir)
    final_directory = os.path.join(mother_dir, signalname)
    os.mkdir(final_directory)
    os.chdir(final_directory)

    
    attributename = signalname+"_Attr"
    metadataname = signalname+ "_meta"
    baselinename = signalname + "_Baseline"
    
    signals[0].to_pickle((signalname+"_P2P.pkl"))
    signals[1].to_pickle((signalname+"_AUC.pkl"))
    baselines[0].to_pickle((baselinename+"P2P.pkl"))
    baselines[1].to_pickle((baselinename+"AUC.pkl"))
    attributes.to_pickle((attributename+".pkl"))
    metadata.to_csv((metadataname+".csv")) # save metadata as csv because can easily open and check
    
    fcs_file_name = signalname + ".fcs"
    file_handle = open(fcs_file_name, 'wb')
    data_to_save = signals[0].to_numpy().flatten().tolist()
    create_fcs(data_to_save, paired_name_list, file_handle)
    file_handle.close()
    os.chdir(mother_dir)
    

if __name__ == "__main__":
    
    ### load and process data, return the processor class ###
    processor = load_and_process(channel_list, pick_ch, paired_name_list, modifier)
    
    ### if no problem, append to our signal buff ###
    signals_buff.append(processor.compile_data(cell_only_filter=1)) # mode 0: all droplets; mode 1: cell_droplet_only
    print("Last processed: ", processor.filehandle)
    
    
    ### if got problem, plot it and change the threshold accordingly ####
    evaluate_result(median_thres, processor)  # copy if you want to show    ,final_plotting = True
    
    modifier = np.array([40,500,1200,1200,80])
    #median_thres[0]+=100
    median_thres=(np.array(median_thres)+modifier).tolist()
    
    ### if you want to know how many cells we collected alr ###
    print("Current count: ", len(pd.concat([result["P2P_data"] for result in signals_buff], ignore_index=True)))
    
    ### finally, export data ###
    export_data(signalname,mother_dir, signals, baselines, signals_buff)

    

    #### Step 0: Trim data to ensure no truncated peaks or abnormal data ####
    #pick_ch = 0
    
    # fig, axs = plt.subplots(3)
    # x = np.arange(len(data)) 
    # head = 2000
    # x1 = np.arange(head)
    # axs[0].plot(x,data[:,pick_ch])
    # axs[1].plot(x1,data[:head,pick_ch])
    # axs[2].plot(x1,data[(len(data)-head):len(data),pick_ch])
        
    #### Step 7: plot scatter for all channel ####
    # channels_to_plot = 5 # 5 channel 
    # mode = 0 # 0 for p2p, 1 for AUC
    
        
    # fig, axs = plt.subplots(channels_to_plot,channels_to_plot,figsize=(15,15))
    # for i in range(channels_to_plot):
    #     for j in range(channels_to_plot):
    #         if i != j:
    #             axs[i,j].scatter(signals[mode].to_numpy()[:,i],signals[mode].to_numpy()[:,j],s=0.5,c='k') #0 for p2p, 1 for auc
    #             axs[i,j].set_xlim([10**2,10**6])
    #             axs[i,j].set_ylim([10**2,10**6])
    #             axs[i,j].set_yscale("log")
    #             axs[i,j].set_xscale("log")

