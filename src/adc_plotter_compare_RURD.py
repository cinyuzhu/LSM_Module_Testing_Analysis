#!/usr/bin/pytOAhon3

"""
Author: Sravan Munagavalasa
Date: December 2023
Description: A set of tools to quickly view the ADC data comming from the ACM
"""

import sys, os, csv
import argparse

import pandas as pd
import numpy as np

from scipy.signal import periodogram

import matplotlib.pyplot as plt

# Read the CSV file into a DataFrame
#fname = '/home/damicm/data/tmp/trace_adc1_bytes100000_2023_12_11_02_21_47.csv'

def main(fname, fname1, trace, psd, t0=0, delta_t=500):
    df = pd.read_csv(fname)
#    print("df ",df)
    df["t"] = df.index / 15 # 15MHz Clock    
    df1=[]
    adctrace1=[]
    if fname1 != "":
        df1 = pd.read_csv(fname1)
        df1["t"] = df1.index / 15 # 15MHz Clock
        adctrace1 = df1.val
    if trace:
        plot_trace(df, df1, t0, delta_t, title_=fname, title_1=fname1)
    
    if psd:
        plot_psd(df.val, adctrace1, npsd=5000, nsampint=150, fs=15.E6, title_=fname, title_1=fname1)

def plot_trace(df, df1, t0=0, delta_t=300, title_="ADC Trace", title_1=""):
    """ plot the adc trace from the start point to end point. Pedestal and signal as well"""
    title = title_
    color1 = "GREY: "
    color2 = "BLUE: "
    # transparency of second file
    transp = 0.5
    df = df[(df.t >= t0) & (df.t <= t0+delta_t)]


#    df_ped = df[df["mode"]==0]
#    df_sig = df[df["mode"]==1]
    df_ped = df[df["RD"]==1]
    df_sig = df[df["RU"]==1]

    plot = False
    if len(df1) >0:
        plot = True
        df1 = df1[(df1.t >= t0) & (df1.t <= t0+delta_t)]
        df_ped1 = df1[df1["mode"]==0]
        df_sig1 = df1[df1["mode"]==1]
        toffset = -(df_ped.iat[0,2]-df_ped1.iat[0,2]) 
        title = color1+": "+title_+"\n"+ color2+": "+title_1

    # Plotting the adc column within the specified range                                                                   
    plt.figure(figsize=(10, 5))
    plt.scatter(df.t,     df.val, color='dimgrey', label='trace_1')
    plt.scatter(df_ped.t, df_ped.val, color='peru', label='ped_1')
    plt.scatter(df_sig.t, df_sig.val, color='forestgreen',   label='sig_1')
    if plot:
        plt.scatter(df1.t-toffset,     df1.val, color='CornflowerBlue', label='trace_2', alpha=transp)
        plt.scatter(df_ped1.t-toffset, df_ped1.val, color='peru', label='ped_2',alpha=transp)
        plt.scatter(df_sig1.t-toffset, df_sig1.val, color='forestgreen',   label='sig_2', alpha=transp)
    # Adding labels and legend                                                                                             
    plt.xlabel('time [$\mu$s]')
    plt.ylabel('ADC Value')
    plt.title(title)
    plt.legend()


def plot_psd(adctrace, adctrace1, npsd=10000, nsampint=150, fs=15.E6, title_=None, title_1=None):
    # colors for the two files
    color1 = 'BLUE'
    color2 = 'RED'
    # transparency of second file
    transp = 0.6
    
    title = title_   
    # Calculate average PSD
    f0, Vxx0 =  make_psd(adctrace, fs)
    f, Vxx =  make_avg_psd(adctrace, fs, npsd)

    plot = False
    if len(adctrace1) > 0:
        plot = True
        f0_1, Vxx0_1 =  make_psd(adctrace1, fs)
        f_1, Vxx_1 =  make_avg_psd(adctrace1, fs, npsd)
        title = color1+": "+title_+"\n"+ color2+": "+title_1
        
    # Generate Figure for Plotting
    fig, axs = plt.subplots(2,2,figsize=(12,7.5))
    if title_:
        fig.suptitle(title, fontsize=16, fontweight='bold')

    # calculate CDS noise
    sig_CDS = calculate_CDS(adctrace, nsampint)
    sig_ADU = np.std(adctrace)
    if plot:
        sig_CDS1 = calculate_CDS(adctrace1, nsampint)
        sig_ADU1 = np.std(adctrace1)

    # Plot the ADC Trace
    ## plot the adc trace in units of  volts
    axs[0,0].plot(adctrace, '.', color=color1)
    if plot:
        axs[0,0].plot(adctrace1, '.', color=color2,  alpha=transp)
#    in the title the rms of the trace and the sigma_CDS 
        axs[0,0].set_title(r'BLUE: $\sigma$={:.2f} ADU;  '.format(sig_ADU)+r' $\sigma$_CDS={:.2e} ADU'.format(sig_CDS)+"\n"+r'RED $\sigma$={:.2f} ADU;  '.format(sig_ADU1)+r' $\sigma$_CDS={:.2e} ADU'.format(sig_CDS1),color="black")
    else:
        axs[0,0].set_title(r'$\sigma$={:.2f} ADU;  '.format(sig_ADU)+r' $\sigma$_CDS={:.2e} ADU'.format(sig_CDS), color="black")
    axs[0,0].set_xlabel('N Sample')
    axs[0,0].set_ylabel('Sample Value (ADU)')


    # Plot PSD in linear scale
    axs[1,1].plot(f[1:], Vxx[1:], color=color1)
    if plot:
        axs[1,1].plot(f_1[1:], Vxx_1[1:], color=color2, alpha=transp)
    axs[1,1].set_ylabel(r"ADU/$\sqrt{Hz}$ 1")
    axs[1,1].set_xlabel("Frequency (Hz)")

    # Plot PSD in log log scale
    axs[1,0].loglog(f[1:], Vxx[1:], color=color1)
    if plot:
        axs[1,0].loglog(f_1[1:], Vxx_1[1:], color=color2,  alpha=transp)
        #axs[1,0].set_title(r'ADC trace $\sigma$={:.2f} ADU;  '.format(sig_ADU1)+r' $\sigma$_CDS={:.2e} ADU'.format(sig_CDS1), color=color2)
    axs[1,0].set_ylabel(r"ADU/$\sqrt{Hz}$")
    axs[1,0].set_xlabel("Frequency (Hz)")
    
    # Plot PSD in log log scale                                                                                        
    axs[0,1].loglog(f0[1:], Vxx0[1:], color=color1)
    if plot:
        axs[0,1].loglog(f0_1[1:], Vxx0_1[1:], color=color2,  alpha=transp)
    axs[0,1].set_ylabel(r"ADU/$\sqrt{Hz}$")
    axs[0,1].set_xlabel("Frequency (Hz)")

    plt.tight_layout()

    #plt.show()


def calculate_CDS(adctrace, nsampint):
    nsamples = len(adctrace)
    nblocks = int(np.ceil(nsamples/(2*nsampint))-1)
    valCDS = np.zeros(nblocks, dtype=np.float64)
    for isampl in range (0,nblocks):
        valCDS[isampl] = np.mean(adctrace[isampl+nsampint:isampl+2*nsampint])-np.mean(adctrace[isampl:isampl+nsampint])
    sigmaCDS = np.std(valCDS)
    return sigmaCDS



def make_psd(adctrace, samplingFreq=15.e6):
    f, Pxx = periodogram(adctrace-np.mean(adctrace), fs=samplingFreq, scaling='density');
    Vxx=np.sqrt(Pxx)
    return f,Vxx

def make_avg_psd(adctrace, samplingFreq=15.e6, npsd=10000):
    Pxx_avg = np.zeros(int(npsd/2)+1, dtype=np.float64)
    nsamples = len(adctrace)
    nblocks = int(np.ceil(nsamples/npsd)-1)
    for isampl in range (0,nblocks):
        istart = isampl*npsd
        iend = istart+npsd
        y_unsigned = adctrace[istart:iend]
        f, Pxx = periodogram(y_unsigned-np.mean(y_unsigned), fs=samplingFreq, scaling='density');
        Pxx_avg = Pxx_avg + Pxx/nblocks

    Vxx=np.sqrt(Pxx_avg)
    return f,Vxx

def calculate_CDS(adctrace,nsampint):
    nsamples = len(adctrace)
    nblocks = int(np.ceil(nsamples/(2*nsampint))-1)
    valCDS = np.zeros(nblocks, dtype=np.float64)
    for isampl in range (0,nblocks):
        valCDS[isampl] = np.mean(adctrace[isampl+nsampint:isampl+2*nsampint])-np.mean(adctrace[isampl:isampl+nsampint])
    sigmaCDS = np.std(valCDS)
    return sigmaCDS


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process some parameters.")
    parser.add_argument('--fname', nargs="+", type=str, help='The file name as a string.')
    parser.add_argument('--trace', action='store_true', help='Enable trace (boolean).')
    parser.add_argument('--psd', action='store_true', help='Enable PSD (boolean).')
    parser.add_argument('--t0', type=float, default=0, help='Start time as a number.')
    parser.add_argument('--deltat', type=float, default=500, help='Time interval as a number.')

    args = parser.parse_args()
    
    f = args.fname[0]
    f1 =""
    if len(args.fname) == 2:
        f1 = args.fname[1]
    main(f, f1, args.trace, args.psd, args.t0, args.deltat)

    plt.show()
