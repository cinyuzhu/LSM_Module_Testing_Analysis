import numpy as np
import ccdmodule_lib_functions_LBC as ccd
from ccdmodule_lib_functions_LBC import *
import ccdmodule_analysis_calibration as ccdcal
from ccdmodule_analysis_calibration import *
import ccdmodule_analysis_initvariables as inivar

def init():    
##############################################################################                                                                   
# Parameters of the LBC Module CCDs                                                                                                    
# 6144 col x 1536 rows                                                                                                            
# 8 extended serial pixels
# so total readout cols for each CCD = 6144 + 8 + 8 = 6160
    global nccds
    nccds = 4
    nrows = ccd.nrows
    ncolCCD = ccd.ncolCCD

    global nimages
    #
    global acmid
    global ccd_i, ccd_label
    ccd_label = ['A','B','C','D']
    global iccd,ccd_calib,ccd_calib_sige,ccd_list

    global nrow_cut_s,nrow_cut_e,ncol_cut_s,ncol_cut_e

    if iccd == "mod1":
#        ccd_calib = [1.246, 1.211, 1.171, 1.129] # LBC pa08
#        ccd_calib = [1.293,1.330,1.383,1.365] # TC2 125
#        ccd_calib = [1.268, 1.220, 1.367, 1.300] # LPNHE 127
        ccd_calib = [1.349,1.322,1.276,1.279]  # TC LSM - CCD PA06
        ccd_calib_sige = np.array([0.16, 0.16, 0.16, 0.16])
    elif iccd == "mod2":
#        ccd_calib = [1.198, 1.235, 1.00, 1.108] # LBC pa07       
        ccd_calib = [1.248,1.32,1.335,1.299] # TC2 124 
#        ccd_calib = [1.216, 1.161, 1.357, 1.288] # LPNHE 127
        ccd_calib_sige = np.array([0.16, 0.16, 0.16, 0.16])
    else:
        print("ONLY CCD MODULES: mod1 or mod2 !!!!")
        exit()        
    ccd_list = [0,1,2,3] # 0=A, 1=B, 2=C, 3=D
    nrows_s = 0 
    nrows_e = nrows 
    # range of rows to be masked because of known noise or problems                                                            
    # for all CCDs                                                                                                             
    nrow_cut_s = 1
    nrow_cut_e = nrows
    # range of columns to be masked because of known noise or problems                                                            
    ncol_cut_s = 10
    ncol_cut_e = ncolCCD
    
    # which analysis to make
    global iPSD, iTrace, iProcImg
    iPSD = True
    iTrace = True
    iProcImg = [True,True,True]   # single kip, 500 skips, 1000 skips
    
######## Calibration ##########
    global sigmaADU,sigma_ee,ncol_rw_s,ncol_rw_e ,nrow_rw_s,nrow_rw_e 
    # approximate resolution in ADU
    sigmaADU = 0.2
    sigma_ee = [0,0,0,0]
    # range used to calculate the pedestal row by row                                                                           
    ncol_rw_s =  10
    ncol_rw_e =  ncolCCD
    nrow_rw_s =  0
    nrow_rw_e =  nrows
    
####### Resolution fit ########
    global calib_fit, sigmae_fit
    calib_fit = [0.,0.,0.,0.]
    sigmae_fit= [0.,0.,0.,0.]

###### Dark current measurement #############
    global dc,edc
    dc = [0.,0.,0.,0.]
    edc = [0.,0.,0.,0.]

########### PLOTTING ###############################
    global iPlotSingle,iPlotSkip
    iPlotSingle = False 
    iPlotSkip = False 
    
    global imagefilename

    global image_in_e_A,image_in_e_B,image_in_e_C,image_in_e_D
    image_in_e_A = np.zeros((nrows, ncolCCD), dtype=np.float64)
    image_in_e_B = np.zeros((nrows, ncolCCD), dtype=np.float64)
    image_in_e_C = np.zeros((nrows, ncolCCD), dtype=np.float64)
    image_in_e_D = np.zeros((nrows, ncolCCD), dtype=np.float64)

    ##### Mask ########
    global cluster_mask_A, cluster_mask_B, cluster_mask_C, cluster_mask_D
    cluster_mask_A = np.full((nrows, ncolCCD), -1, dtype=np.int32)
    cluster_mask_B = np.full((nrows, ncolCCD), -1, dtype=np.int32)
    cluster_mask_C = np.full((nrows, ncolCCD), -1, dtype=np.int32)
    cluster_mask_D = np.full((nrows, ncolCCD), -1, dtype=np.int32)
    return 0

def read_mask(ccd_i):
    if ccd_i==0:
        cluster_mask = inivar.cluster_mask_A
    elif ccd_i==1:
        cluster_mask = inivar.cluster_mask_B
    elif ccd_i==2:
        cluster_mask = inivar.cluster_mask_C
    elif ccd_i==3:
        cluster_mask = inivar.cluster_mask_D
    return cluster_mask

def read_data(ccd_i):
    if ccd_i==0:
        image_data = image_in_e_A
    elif ccd_i==1:
        image_data = image_in_e_B
    elif ccd_i==2:
        image_data = image_in_e_C
    elif ccd_i==3:
        image_data = image_in_e_D
    return image_data

def load_newmask(ccd_i,cluster_mask):
    if ccd_i==0:
        inivar.cluster_mask_A = cluster_mask
    elif ccd_i==1:
        inivar.cluster_mask_B = cluster_mask
    elif ccd_i==2:
        inivar.cluster_mask_C = cluster_mask
    elif ccd_i==3:
        inivar.cluster_mask_D = cluster_mask
    return 0
