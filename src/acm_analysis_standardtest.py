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

# directory where files are stored
#directory_path = "../Decoder/data/acm_TC2REF/"
#directory_path = "/Users/priviter/Work/DAMIC/LBC/UsefulCode/data/2025-06-22-acmstdtest1/"
directory_path = "/Users/priviter/Work/DAMIC/LBC/UsefulCode/data/test/test1/"
# directory where analysis results are stored
results_path = ""
# which analysis to make and plots to show                                                                                               
iPSD = False #True
iTrace = False#True
iProcImg = [True,True,True]   # single kip, 500 skips, 1000 skips  
iPlotSingle = False#True
iPlotSkip = False#True

def read_files(directory_path):
    files_only = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]
    for i in range(0,len(files_only)):
#        print(files_only[i])
        if "avg_img_conv_100x6300x1_bin1x1_"+acm1 in files_only[i]:
            imagefile_acm1[0] = files_only[i]
        if "avg_img_conv_10x630x500_bin1x1_"+acm1 in files_only[i]:
            imagefile_acm1[1] = files_only[i]
        if "avg_img_conv_10x630x1000_bin1x1_"+acm1 in files_only[i]:
            imagefile_acm1[2] = files_only[i]
        if "trace_ch0_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            tracefile_acm1[0] = files_only[i]
        if "trace_ch1_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            tracefile_acm1[1] = files_only[i]
        if "trace_ch2_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            tracefile_acm1[2] = files_only[i]
        if "trace_ch3_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            tracefile_acm1[3] = files_only[i]
        if "psd_ch0_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            psdfile_acm1[0] = files_only[i]
        if "psd_ch1_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            psdfile_acm1[1] = files_only[i]
        if "psd_ch2_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            psdfile_acm1[2] = files_only[i]
        if "psd_ch3_"+acm1 in files_only[i] and 'csv' in files_only[i]:
            psdfile_acm1[3] = files_only[i]
    ### find CCD module corresponding to this acm
    #substring = "mod"
    #index = imagefile_acm1[0].find(substring)
    #if index != -1:
    #    iccd = imagefile_acm1[0][index:index+4]        
    #else:
    #    print("CCD Module not found in file name!")
    #    exit()

    return 0

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
    else:
        inivar.sigmaADU = 0.2
        inivar.ncol_rw_s = 100 # to avoid transient at the beginning of image
    for ccd_i in ccd_list:
        ccdid = inivar.ccd_label[ccd_i]
        print("********************************************************************************")
        print("***************************** CCD ",ccdid," ******************************************")
        print("********************************************************************************")
        if img==0:
            inivar.iPlotSingle = iPlotSingle            
        # standard calibration
        ccdcal.calibrate_ccdimage(ccd_i)
        if img==0:
            sigma_ee_single[ccd_i] = inivar.sigma_ee[ccd_i]
            sigma_ADU_single[ccd_i] = inivar.sigma_ee[ccd_i]*inivar.ccd_calib[ccd_i]
        # fit pixel charge distribution to get dark current
        if img>0:
            inivar.iPlotSingle = False
            inivar.iPlotSkip = iPlotSkip
            ccdplt.fit_singlee_resolution(ccd_i)
        if img==1:
            calib_fit_500[ccd_i]=inivar.calib_fit[ccd_i]
            sigmae_fit_500[ccd_i]=inivar.sigmae_fit[ccd_i]
        elif img==2:
            calib_fit_1000[ccd_i]=inivar.calib_fit[ccd_i]
            sigmae_fit_1000[ccd_i]=inivar.sigmae_fit[ccd_i]
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
        print("CCD ",ccdid,":  500 skips - Calibration: %.4g" % calib_fit_500[ccd_i]," ADU/e- ; Resolution: %.3g" % sigmae_fit_500[ccd_i]," e-",file=fileout)
        print("CCD ",ccdid,": 1000 skips - Calibration: %.4g" % calib_fit_1000[ccd_i]," ADU/e- ; Resolution: %.3g" % sigmae_fit_1000[ccd_i]," e-",file=fileout) 
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
calib_fit_1000 = [0.,0.,0.,0.]
sigmae_fit_1000 = [0.,0.,0.,0.]
sigma_ee_single = [0.,0.,0.,0.]
sigma_ADU_single = [0.,0.,0.,0.]
global imagefile_acm1,tracefile_acm1,psdfile_acm1
imagefile_acm1 = ["","",""]
tracefile_acm1 = ["","","",""]
psdfile_acm1 = ["","","",""]

#
original_stdout = sys.stdout # Save a reference to the original standard output
#

# Read the files in the directory
read_files(directory_path)
print("ACM ",acm1," Image files: ",imagefile_acm1)
print("ACM ",acm1," Trace files: ",tracefile_acm1)
print("ACM ",acm1," PSD files: ",tracefile_acm1)
#

#
# process Images for ACM Module
global fileout
inivar.acmid = acm1
print("*********** STARTING PROCESSING ACM ",acm1," ",iccd,"**********************************")
for image in range(0,3):
    imagefiles = imagefile_acm1[image]
    if imagefiles !="" and iProcImg[image]: process_acm(directory_path+imagefiles,iccd,image)
if imagefile_acm1==["","",""]:
    print("NO FILES for ACM ",acm1," !!!")
else:
    # save results in file
    fileresults = results_path+"Summary_"+acm1+"_"+iccd+".txt"
    print("Results for ACM ",acm1," in file ",fileresults)
    fileout = open(fileresults, "w")
    save_results(acm1,iccd)
    print("CHECK tracks in single-skip image: ",directory_path+imagefile_acm1[0],file=fileout)
    print("************************************************************************************",file=fileout)
    print("**************************** ADC TRACES ********************************************",file=fileout)
    print("Normal?  YES/NO",file=fileout)
    print("COMMENTS: ",file=fileout)
    print(" ",file=fileout)
    print(" ",file=fileout)
    print("************************************************************************************",file=fileout)
    print("**************************** PSD  ********************************************",file=fileout)
    print("Normal?  YES/NO",file=fileout)
    print("COMMENTS: ",file=fileout)
    print(" ",file=fileout)
    print(" ",file=fileout)
    print("************************************************************************************",file=fileout)
    print("**************************** TRACKS  ********************************************",file=fileout)
    print("Normal?  YES/NO",file=fileout)
    print("COMMENTS: ",file=fileout)
    print(" ",file=fileout)
    print(" ",file=fileout)   
    fileout.close()
    
print("*********** COMPLETED PROCESSING ACM ",acm1," ",iccd,"**********************************")

# print the results on stdout
print(" ")
print(" ")
with open(fileresults, 'r') as file:
    # Read the content of the file
    file_content = file.read()    
    # Print the content
    print(file_content)

# Now look at the traces
if iTrace:
    # Display traces 
    ref_trace=["","","",""] # reference files to compare with
    trace = True
    psd = False
    t0 = 0
    delta_t = 5000
    for ich in range(0,4):
        if tracefile_acm1[ich] != "":
            adcplot.main(directory_path+tracefile_acm1[ich], ref_trace[ich], True, False, t0, delta_t)
    plt.show()

# Now look at the PSD
if iPSD:
    ref_psd=["","","",""] # reference files to compare with
    # Display psd 
    for ich in range(0,4):
        if psdfile_acm1[ich] != "":
            adcplot.main(directory_path+psdfile_acm1[ich], ref_psd[ich], False, True, t0, delta_t)
    plt.show()
