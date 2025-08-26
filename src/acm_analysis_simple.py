# -*- coding: utf-8 -*-
"""
=======================================
Read a skipper-averaged image from a FITS file, calibrate from ADU to e-, find clusters, measure dark current
=======================================

*By P. Privitera 
-------------------

"""

import time
import os
# input values from command line                                                                                                         
import sys
# acm id 
acm1 = sys.argv[1]
# ccd module: mod1 or mod2
iccd = sys.argv[2]
# file to analyze
imgfile = sys.argv[3]

##############################################################################                             
### import ccdmodule_analysis packages
import ccdmodule_analysis_calibration as ccdcal
from ccdmodule_analysis_calibration import *
import ccdmodule_analysis_makeplots as ccdplt
from ccdmodule_analysis_makeplots import *
import ccdmodule_analysis_initvariables as inivar
import ccdtest_lib_plots
from ccdtest_lib_plots import *
import ccdmodule_lib_functions_LBC as ccd
from ccdmodule_lib_functions_LBC import *
import adc_plotter_compare_RURD as adcplot

img = 1 # 1 = multiskip img ; 0 = single skip 
iPlotSingle = True #True
iPlotSkip = True #True

def make_mask():
    chargesat = 100
    sat_mask = np.zeros((ccd.nrows, ccd.ncolCCD), dtype=np.int32)
    ccd_list = inivar.ccd_list
    for ccd_i in ccd_list:
        ccdid = inivar.ccd_label[ccd_i]
        cluster_mask = inivar.read_mask(ccd_i)
        image_data = inivar.read_data(ccd_i)
        # mask saturated pixels
        sat_mask[image_data>chargesat]=-100
        # mask rows according to init input                                                                                       
        cluster_mask[0:inivar.nrow_cut_s,:] = -90
        cluster_mask[inivar.nrow_cut_e:,:] = -90
        cluster_mask[:,0:inivar.ncol_cut_s] = -90
        cluster_mask[:,inivar.ncol_cut_e:] = -90
        inivar.load_newmask(ccd_i,cluster_mask)
    for ccd_i in ccd_list:
        cluster_mask = inivar.read_mask(ccd_i)
        cluster_mask = cluster_mask + sat_mask
        inivar.load_newmask(ccd_i,cluster_mask)
    
def process_acm(imagefiles,iccd,img):
    # READ IMAGES
    inivar.imagefilename = imagefiles
    print("****** Read image ",imagefiles)
    readimage = ccd.ReadCCDImage_ACM(imagefiles, "data")
    # initialize variables
    inivar.iccd = iccd
    inivar.init()
    ccd_list = inivar.ccd_list
    if img == 0:
        inivar.sigmaADU = 5.
        inivar.ncol_rw_s = 300 # to avoid transient at the beginning of image
        inivar.ncol_cut_s = 617 # to avoid transient at the beginning of image
    else:
        inivar.sigmaADU = 0.2

    #calibrate
    for ccd_i in ccd_list:
        ccdid = inivar.ccd_label[ccd_i]
        print("********************************************************************************")
        print("***************************** CCD ",ccdid," ******************************************")
        print("********************************************************************************")
        if img==0:
            inivar.iPlotSingle = iPlotSingle            
        # standard calibration
        ccdcal.calibrate_ccdimage(ccd_i)

    # make mask
    make_mask()
    
    for ccd_i in ccd_list:
        ccdid = inivar.ccd_label[ccd_i]
        print("********************************************************************************")
        print("***************************** CCD ",ccdid," ******************************************")
        print("********************************************************************************")
        if img==0:
            inivar.iPlotSingle = iPlotSingle            
        if img==0:
            sigma_ee_single[ccd_i] = inivar.sigma_ee[ccd_i]
            sigma_ADU_single[ccd_i] = inivar.sigma_ee[ccd_i]*inivar.ccd_calib[ccd_i]
        # fit pixel charge distribution to get single e resolution
        if img>0:
            inivar.iPlotSingle = False
            inivar.iPlotSkip = iPlotSkip
            ccdplt.fit_singlee_resolution(ccd_i)
            calib_fit_500[ccd_i]=inivar.calib_fit[ccd_i]
            sigmae_fit_500[ccd_i]=inivar.sigmae_fit[ccd_i]
        plt.show()
    return 0

def save_results(acm,iccd):
    print("********************************************************************************",file=fileout)
    print("***************************** FINAL RESULTS ******************************************",file=fileout)
    print("*************************** ACM ",acm,"  ",iccd,"  *****************************************************",file=fileout)
    print("********************************************************************************",file=fileout)
    ccd_list = inivar.ccd_list
    for ccd_i in ccd_list:
        ccdid = inivar.ccd_label[ccd_i]
        print("CCD ",ccdid,": Single skip resolution: %.3g" % sigma_ADU_single[ccd_i]," ADU ;  %.3g" % sigma_ee_single[ccd_i]," e-",file=fileout) 
        print("CCD ",ccdid,":  Multi skips - Calibration: %.4g" % calib_fit_500[ccd_i]," ADU/e- ; Resolution: %.3g" % sigmae_fit_500[ccd_i]," e-",file=fileout)
        print("********************************************************************************",file=fileout)
    return 0

###############################################################################
########  PROGRAM STARTS HERE
###############################################################################
# for the results summary
global calib_fit,sigmae_fit,calib_fit_500,sigmae_fit_500,calib_fit_1000,sigmae_fit_1000,sigma_ee_single,sigma_ADU_single
calib_fit = [0.,0.,0.,0.]
sigmae_fit = [0.,0.,0.,0.]
calib_fit_500 = [0.,0.,0.,0.]
sigmae_fit_500 = [0.,0.,0.,0.]
sigma_ee_single = [0.,0.,0.,0.]
sigma_ADU_single = [0.,0.,0.,0.]
#
original_stdout = sys.stdout # Save a reference to the original standard output
#

#
# process Images for ACM Module
global fileout
fileout = sys.stdout
inivar.acmid = acm1
print("*********** STARTING PROCESSING ACM ",acm1," ",iccd,"**********************************")

imagefiles = imgfile
process_acm(imagefiles,iccd,img)

save_results(acm1,iccd)
    
print("*********** COMPLETED PROCESSING ACM ",acm1," ",iccd,"**********************************")

