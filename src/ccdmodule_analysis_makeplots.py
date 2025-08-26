# -*- coding: utf-8 -*-
"""
=======================================
Make plots
=======================================

*By P. Privitera 
-------------------

"""
##############################################################################                             
# Get Numpy, and relevant Scipy                                                                            
##############################################################################                             
import numpy as np
import numpy.ma as ma
import scipy.fftpack
from scipy import stats
from scipy.stats import norm, poisson
from scipy import signal
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import truncnorm
from math import exp, pi, sqrt, log

# Set up matplotlib and use a nicer set of plot parameters
import matplotlib.pyplot as plt

import ccdtest_lib_plots
from ccdtest_lib_plots import *

import ccdmodule_lib_functions_LBC as ccd
from ccdmodule_lib_functions_LBC import *
import ccdmodule_analysis_initvariables as inivar
from ccdmodule_analysis_initvariables import *
##############################################################################
# get iMinuit
# everything in iminuit is done through the Minuit object, so we import it
from iminuit import Minuit
from math import exp, pi, sqrt, sin, cos, atan, log
from scipy.special import factorial

##############################################################################
# FUNCTIONS
#############################################################################
def one_gauss(mu, sigma):
    funzione = np.zeros(len(binspectr))
    funzione = funzione + np.exp(-0.5*(binspectr-mu)**2/sigma**2)/(sqrt(2*pi)*sigma)
    return funzione
def one_gauss_p(x,par):
    mu = par[0]
    sigma = par[1]
    funzione = np.zeros(len(x))
    funzione = funzione + np.exp(-0.5*(x-mu)**2/sigma**2)/(sqrt(2*pi)*sigma)
    return funzione

def NLL_onegauss(mu, sigma):
    arr = one_gauss(mu, sigma)
    arr[arr<=0.]= 0.0001
    temp = nobs * np.log(arr)
    nll =  np.sum(arr) - np.sum(temp)
    return nll

def my_gauss_offset_f(mu, sigma, cal, off):
    funzione = np.zeros(len(binspectr))
    for i in range(0,2):
        #pamp = poisson.pmf(i, mu)
        # if you want a poisson fit comment pamp below
        if i==0:
            pamp = 1.-mu
        else:
            pamp = mu
        mui = float(i)*cal-off
        funzione = funzione + pamp*np.exp(-0.5*(binspectr-mui)**2/sigma**2)/(sqrt(2*pi)*sigma)
    return funzione

def my_gauss_offset_f_p(x,par):
    funzione = np.zeros(len(x))
    for i in range(0,2):
#        pamp = poisson.pmf(i, par[0])
        # if you want a poisson fit comment pamp below
        if i==0:
            pamp = 1.-par[0]
        else:
            pamp = par[0]
        mui = float(i)*par[2]-par[3]
        sigma = par[1]
        funzione = funzione + pamp*np.exp(-0.5*(x-mui)**2/sigma**2)/(sqrt(2*pi)*sigma)
    return funzione

def NLL(mu, sigma, cal, off):
    arr = my_gauss_offset_f(mu, sigma, cal, off)
    arr[arr<=0.]= 0.0001 #0.01
    temp = nobs * np.log(arr)
    nll =  np.sum(arr) - np.sum(temp)
    return nll

def fit_singlee_resolution(ccd_i):
    global nobs,binspectr
    ccdid = inivar.ccd_label[ccd_i]
    image_data = inivar.read_data(ccd_i)
    cluster_mask = inivar.read_mask(ccd_i)

#    spectrum = image_data[(cluster_mask==-1) & (image_data>-3.) & (image_data<10)]
    sigma_e = inivar.ccd_calib_sige[ccd_i]
    cutehigh = 1.+3*sigma_e
    spectrum = image_data[(cluster_mask==-1) & (image_data>-0.6) & (image_data<cutehigh)]
    nobs,bins = np.histogram(spectrum,1000,range=(-1.,cutehigh))
    binspectr = bins[:-1] + np.diff(bins)/2 
    m = Minuit(NLL, mu=0.05, sigma=sigma_e, cal=1., off=0.0)
#    m = Minuit(NLL, mu=0.0005, sigma=sigma_e, cal=1., off=0.0)

    m.print_level = 1
    m.errordef = 0.5
    m.errors = (0.001,sigma_e*0.1,0.01,0.01)
    m.limits = [(0., None), (0.05,0.3),(0.,None),(None,None)]
##    m.fixto('sigma',sigma_e)
##    m.fixto('cal',1.)
#    m.fixto('off',0.)
    
    m.migrad()
    print("               xxxx Single e resolution: CCD ",ccdid," xxxxxxxxxxxxx")
    for p in m.parameters:
        print("{} = {:.5f} +/- {:.5f}".format(p, m.values[p], m.errors[p]))
    xlam = m.values["mu"]
    dxlam = m.errors["mu"]
    inivar.dc[ccd_i] = xlam
    inivar.edc[ccd_i] = dxlam
# now make plots
    par = np.zeros(4)
    par[0] = xlam
    par[1] = m.values["sigma"]
    par[2] = m.values["cal"]
    par[3] = m.values["off"]

    if inivar.iPlotSkip:
        hist = Plot1DHisto(spectrum,100,-0.6,1.5,'bar','lin',my_gauss_offset_f_p,par)
        ###    hist = Plot1DHisto(spectrum,1000,-1,4.5,'bar','log',0,0)
        PlotLabels(hist,'ACM '+inivar.acmid+' - CCD '+ccdid+' Charge distribution '+str(ccd.nskips)+' skips','Charge (e-)','N pixels','')    
        #plt.yscale('log')

    inivar.calib_fit[ccd_i] = inivar.ccd_calib[ccd_i]*par[2]
    inivar.sigmae_fit[ccd_i] = par[1]/par[2]
    return 0

def fit_singleskip_resolution(ccd_i,spectrum):
    global nobs,binspectr
    ccdid = inivar.ccd_label[ccd_i]
    sigma_e = inivar.sigmaADU
    xmin = min(spectrum)
    xmax = max(spectrum)
    nobs,bins = np.histogram(spectrum,1000,range=(xmin,xmax))
    binspectr = bins[:-1] + np.diff(bins)/2
    # minuit fit
    m = Minuit(NLL_onegauss, mu=0.0005, sigma=sigma_e)

    m.print_level = 1
    m.errordef = 0.5
    m.errors = (0.001,sigma_e*0.1)
    m.limits = [(0., None), (0.05,10.)]
##    m.fixto('sigma',sigma_e)
    
    m.migrad()
    print("               xxxx Single Skip resolution: CCD ",ccdid," xxxxxxxxxxxxx")
    for p in m.parameters:
        print("{} = {:.5f} +/- {:.5f}".format(p, m.values[p], m.errors[p]))
# now make plots
    par = np.zeros(2)
    par[0] = m.values["mu"]
    par[1] = m.values["sigma"]

    if inivar.iPlotSingle:
        xmin = min(spectrum)
        xmax = max(spectrum)
        hist = Plot1DHisto(spectrum,100,xmin,xmax,'bar','lin',one_gauss_p,par)
        PlotLabels(hist,'ACM '+inivar.acmid+' - CCD '+ccdid+' Charge distribution Single skip','Charge (e-)','N pixels','')    
        #plt.yscale('log')
        plt.show()

    return par[0],par[1]
