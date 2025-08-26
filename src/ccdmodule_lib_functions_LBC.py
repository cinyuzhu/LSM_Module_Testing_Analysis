##############################################################################                             
###  THESE ARE COMMON FUNCTION USED FOR THE CCD TEST PYTHON PACKAGE
###  by P. Privitera
##############################################################################                               
##############################################################################                                                    
# used to calculate time, before and after an operation: time.time()                                                              
import time
from datetime import datetime
##############################################################################                                                    
# Get Numpy, and relevant Scipy                                                                                                   
##############################################################################                                                    
import numpy as np
import scipy.fftpack
from scipy import stats
from scipy.stats import norm, poisson
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
##############################################################################                                                    
# Set up matplotlib and use a nicer set of plot parameters                                                                        
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
##############################################################################                                                    
# get the FITS package from astropy                                                                                               
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
##############################################################################   
import ccdmodule_analysis_initvariables as inivar
##############################################################################   

def ReadCCDImage_ACM(imagename, imgtyp):
    global image_data_A,image_data_B,image_data_C,image_data_D
    global image_var_A,image_var_B,image_var_C,image_var_D
    global ccd_mask_A,ccd_mask_B,ccd_mask_C,ccd_mask_D
    # Open the data image
    image_file = fits.open(imagename)

    # For a CCD module with Dsub-50 connector:  LBC data
    if imgtyp =="data":
        image_data_A = image_file[1].data
        image_data_B = image_file[2].data
        image_data_D = image_file[3].data
        image_data_C = image_file[4].data
        # For a CCD module with Samtec connector: to be checked
        #image_data_D = image_file[1].data
        #image_data_C = image_file[2].data
        #image_data_A = image_file[3].data
        #image_data_B = image_file[4].data
    elif imgtyp == "var":
        # For a CCD module with Dsub-50 connector                                                                                    
        image_var_A = image_file[1].data
        image_var_B = image_file[2].data
        image_var_D = image_file[3].data
        image_var_C = image_file[4].data
    elif imgtyp == "mask":
        # For a CCD module with Dsub-50 connector                                                                                
        ccd_mask_A = image_file[1].data
        ccd_mask_B = image_file[2].data
        ccd_mask_D = image_file[3].data
        ccd_mask_C = image_file[4].data
    else:
        print("UNKNOWN IMAGE TYPE!!!")
        exit()        
        
    global hdr, nallcolumns, nrows, ncolCCD, nskips, npbin, nsbin, ampl, exposure_time, read_time, readstart

    hdr = image_file[0].header

    print(repr(hdr))
    print()

    acm_nrow = int(hdr['NROW'])
    acm_ncol = int(hdr['NCOL'])
    acm_nsamp = int(hdr['NDCM'])

    # 
    nallcolumns = acm_ncol * acm_nsamp
    nrows = acm_nrow
    ncolCCD = acm_ncol
    nskips = acm_nsamp
    npbin = int(hdr['NPBIN'])
    nsbin = int(hdr['NSBIN'])

    exposure_time = float(hdr['EXPOSURE'])
    
    dateini = hdr['DATEINI']
    dateend = hdr['DATEEND']
    readstart = datetime.fromisoformat(dateini)
    readend = datetime.fromisoformat(dateend)
    readtime = readend-readstart
    read_time = readtime.total_seconds()-exposure_time

    # total time in days                                                                                                                        
    global tottime, ncolumns
    tot_time = (exposure_time+read_time)/86400.
    ncolumns = acm_ncol

    return 0

def MakeAvgSkipImg(iskipstart,iskipend,image_data):
    skipper_avg = np.zeros((nrows, ncolumns), dtype=np.float64)
#    t0 = time.time()
    xp = -1
    for x in range(0, nallcolumns, nskips):
        xp = xp+1
        xeffstart = x + iskipstart
        xeffend = x + iskipend
        skipper_avg[:,xp] = np.mean(image_data[:,xeffstart:xeffend],axis=1)
#    t1 = time.time()
#    print("Time to make avg skip image ",t1-t0)
    return skipper_avg

def MakeRmsSkipImg(iskipstart,iskipend,image_data):
    skipper_rms = np.zeros((nrows, ncolumns), dtype=np.float64)
#    t0 = time.time()                                                                                        
    xp = -1
    for x in range(0, nallcolumns, nskips):
        xp = xp+1
        xeffstart = x + iskipstart
        xeffend = x + iskipend
        skipper_rms[:,xp] = np.std(image_data[:,xeffstart:xeffend],axis=1)
#    t1 = time.time()                                                                                      
#    print("Time to make rms skip image ",t1-t0)     
    return skipper_rms

############################################################################## 
# END OF FUNCTIONS                                                                                           
############################################################################## 
