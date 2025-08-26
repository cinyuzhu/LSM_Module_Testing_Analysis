# -*- coding: utf-8 -*-
"""
=======================================
Calibrate a CCD image, from ADU to e- 
=======================================

*By P. Privitera 
-------------------

"""
##############################################################################                             
# Get Numpy, and relevant Scipy                                                                            
##############################################################################                             
import numpy as np
import scipy.fftpack
from scipy import stats
from scipy.stats import norm, poisson
from math import exp, pi, sqrt

import ccdtest_lib_plots
from ccdtest_lib_plots import *

import ccdmodule_lib_functions_LBC as ccd
from ccdmodule_lib_functions_LBC import *
import ccdmodule_analysis_makeplots as ccdplt
from ccdmodule_analysis_makeplots import *

import ccdmodule_analysis_initvariables as inivar
##############################################################################
# FUNCTIONS
#############################################################################
def calibrate_simple(data, ADUlow, ADUhigh):
# determine the zero electron ADU
    mu, sigma = norm.fit(data[(data < ADUhigh) & (data > ADUlow)].flatten())
    return mu,sigma

##############################################################################                             

def calibrate_ccdimage(ccd_i):
##############################################################################    
    # approximate resolution in ADU                                                                                                      
    sigmaADU = inivar.sigmaADU 
    ccd_id = inivar.ccd_label[ccd_i]
    calib = inivar.ccd_calib[ccd_i]
    if ccd_i==0:
        image_data0 = ccd.image_data_A
    elif ccd_i==1:
        image_data0 = ccd.image_data_B
    elif ccd_i==2:
        image_data0 = ccd.image_data_C
    elif ccd_i==3:
        image_data0 = ccd.image_data_D

#    sigma_ADU = sigmaADU
# calibrate the image: 
# pedestal = ADU value of zero electron peak
# calib_data = calibration constant ADU per electron
# sigmaADU_data = avg sigma of single electrons peaks in ADU 
# fix to sigmaADU
    sigmaADU_data = sigmaADU
    calib_data = calib
    sigma_e = sigmaADU_data/calib
# now fill the image                                                             
# image with pedestal subtracted and charge positive; calibrated in e-
    image_data = np.zeros((ccd.nrows, ccd.ncolCCD), dtype=np.float64)
# row by row calibration
    avgped = 0
    nsigma = 2
    for y in range(0,ccd.nrows):
        ADUmedian = np.median(image_data0[y,inivar.ncol_rw_s:inivar.ncol_rw_e].flatten())
        ADUlow = ADUmedian-nsigma*sigma_e*calib_data
        ADUhigh = ADUmedian+nsigma*sigma_e*calib_data
        pedestal_rw, sigmaADU_rw = calibrate_simple(image_data0[y,inivar.ncol_rw_s:inivar.ncol_rw_e], ADUlow, ADUhigh)
        image_data[y,:] = (image_data0[y,:]-pedestal_rw)/calib_data
        avgped = avgped + pedestal_rw
    avgped = avgped/ccd.nrows

    spectrum = image_data[ (image_data<2.5*sigma_e) & (image_data>-2.5*sigma_e) ]*calib_data # in ADU
#    imagecut = image_data[:,inivar.ncol_rw_s:]
#    spectrum = imagecut[ (imagecut<4.*sigma_e) & (imagecut>-4.*sigma_e) ]*calib_data # in ADU
#    spectrum = image_data[ (image_data<4.*sigma_e) & (image_data>-4.*sigma_e) ]*calib_data # in ADU
    mu_fit,sigma_fit = ccdplt.fit_singleskip_resolution(ccd_i,spectrum)
    inivar.sigma_ee[ccd_i] = sigma_fit/calib_data #np.std(image_data[ (image_data<2.*sigma_e) & (image_data>-2.*sigma_e) ] )
    sigmaADU_data = inivar.sigma_ee[ccd_i]*calib_data
    print("               xxxxxxxxxxx Calibration: CCD ",ccd_id," xxxxxxxxxxxxx")
    print("Pedestal 0 electrons: ",avgped," ADU" )
    print("CALIBRATION CONSTANT: ",calib_data," ADU/e- ; SIGMA (ADU):",sigmaADU_data," ; RESOLUTION: ",inivar.sigma_ee[ccd_i]," e-" )
    # fill the corresponding image in elecrons                                                                                              
    if ccd_i==0:
        inivar.image_in_e_A = image_data
    elif ccd_i==1:
        inivar.image_in_e_B = image_data
    elif ccd_i==2:
        inivar.image_in_e_C = image_data
    elif ccd_i==3:
        inivar.image_in_e_D = image_data
    return 0
 
#############################################################################
# END OF FUNCTIONS     
