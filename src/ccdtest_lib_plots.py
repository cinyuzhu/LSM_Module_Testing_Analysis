# Set up matplotlib and use a nicer set of plot parameters                                                                  
import matplotlib
import matplotlib.pyplot as plt
matplotlib.get_backend()
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import numpy as np
import scipy.fftpack
from scipy import stats
from math import exp, pi, sqrt, sin, cos, atan
###############################################################################3

def Plot1DHisto(data,nbins,xlow,xhigh,htype,Log,pdf,par):
# htype="bar" or " " 
# 1D histogram
    global ax,fig
    fig = plt.figure(figsize=[10,5])
    ax = plt.subplot(111)
# create the binned histogram
    databinned, bins = np.histogram(data,nbins,range=[xlow,xhigh])
    centers = bins[:-1] + np.diff(bins)/2 
    if htype=="bar":
        ax.hist(bins[:-1], bins, weights=databinned, histtype='bar')
    else:
        errors = np.sqrt(databinned)
        ax.errorbar(centers, databinned, yerr=errors, fmt='bo', ms=4)

    if Log=='log':
        plt.yscale('log')

    if pdf!=0:
# superimpose function
        ndatatot = np.sum(databinned)
        x = np.linspace(xlow,xhigh,nbins*20)
        binwidth = bins[1]-bins[0]
        pdf_int = np.sum(pdf(x,par))*(x[1]-x[0])
        pdf_norm = pdf(x,par)/pdf_int * ndatatot * binwidth 
        ax.plot(x,pdf_norm,'r-')

    return ax

def Plot2DHisto(datax,datay,errorx,errory,Logx,Logy,pdf,par):
# 2D plot                                                                                                                    
    fig = plt.figure(figsize=[10,5])
    ax = plt.subplot(111)
    ax.errorbar(datax, datay, xerr= errorx, yerr=errory, fmt='bo', ms=4)

    if Logx=='log':
        plt.xscale('log')
    if Logy=='log':
        plt.yscale('log')

    if pdf!=0:
# superimpose function   
        xlow = np.amin(datax)
        xhigh = np.amax(datax)
        x = np.linspace(xlow,xhigh,2000)
        ax.plot(x,pdf(x,par),'r-')

    return ax

def PlotProfileHisto(datax,datay,errorx,errory,nbinprof,Logx,Logy,pdf,par):
# Profile plot                                     
    fig = plt.figure(figsize=[10,5])
    ax = plt.subplot(111)

    xavgbin, xstdbin, yavgbin, ystdbin, npixbin =  profile_w(datax,datay,errorx,errory,nbinprof) 
#    xavgbin, xstdbin, yavgbin, ystdbin, npixbin =  profile(datax,datay,nbinprof)
    ax.errorbar(xavgbin, yavgbin,
            xerr=xstdbin,                                                                                                                                                                             
            yerr=ystdbin,
                 fmt='bo', ms=4)

    if Logx=='log':
        plt.xscale('log')
    if Logy=='log':
        plt.yscale('log')

    if pdf!=0:
# superimpose function                                                                                                                                                                                                                                       
        xlow = np.amin(datax)
        xhigh = np.amax(datax)
        x = np.linspace(xlow,xhigh,2000)
        ax.plot(x,pdf(x,par),'r-')

    return ax

def profile(x,y,nbin):
    xavgbin = np.zeros(nbin, dtype=np.float64)
    xstdbin = np.zeros(nbin, dtype=np.float64)
    yavgbin = np.zeros(nbin, dtype=np.float64)
    ystdbin = np.zeros(nbin, dtype=np.float64)
    xbincenter = np.zeros(nbin, dtype=np.float64)
    xmin = np.amin(x)
    xmax = np.amax(x)
    binwidth = (xmax-xmin)/float(nbin)
    for ibin in range(0,nbin):
        xlow = xmin+ibin*binwidth
        xhigh = xlow + binwidth
        xbincenter[ibin] = xlow+binwidth/2.
        cutbin = (x>=xlow) & (x<xhigh)
        yavgbin[ibin] = np.mean(y[cutbin])
        ystdbin[ibin] = np.std(y[cutbin])/np.sqrt(len(y[cutbin]))
        xavgbin[ibin] = np.mean(x[cutbin])
        xstdbin[ibin] = np.std(x[cutbin])/np.sqrt(len(x[cutbin]))
    return xavgbin, xstdbin, yavgbin,ystdbin,xbincenter

def profile_w(x,y,sigx,sigy,nbin):
    xavgbin = np.zeros(nbin, dtype=np.float64)
    xstdbin = np.zeros(nbin, dtype=np.float64)
    yavgbin = np.zeros(nbin, dtype=np.float64)
    ystdbin = np.zeros(nbin, dtype=np.float64)
    xbincenter = np.zeros(nbin, dtype=np.float64)
    xmin = np.amin(x)
    xmax = np.amax(x)
    binwidth = (xmax-xmin)/float(nbin)
    for ibin in range(0,nbin):
        xlow = xmin+ibin*binwidth
        xhigh = xlow + binwidth
        xbincenter[ibin] = xlow+binwidth/2.
        cutbin = (x>=xlow) & (x<xhigh)
        weig = 1./(sigy[cutbin]*sigy[cutbin])
        yavgbin[ibin] = np.average(y[cutbin],weights=weig)
        ystdbin[ibin] = 1./sqrt(np.sum(weig)) #(np.std(y[cutbin])/np.sqrt(len(y[cutbin]))
        xavgbin[ibin] = np.mean(x[cutbin])
        xstdbin[ibin] = np.std(x[cutbin])/np.sqrt(len(x[cutbin]))
    return xavgbin, xstdbin, yavgbin,ystdbin,xbincenter



def PlotLabels(ax,HTitle,Xlabel,Ylabel,filename):
    ax.set_title(HTitle, fontsize= 16, fontweight="bold")
    ax.set_xlabel(Xlabel)
    ax.set_ylabel(Ylabel)
    if filename != '':
        fig.savefig(filename)
    return 0

def Plot1DHistoS(data,nbins,xlow,xhigh,htype,Log,pdf,par):
# htype="bar" or " " 
# 1D histogram
#    fig = plt.figure(figsize=[10,5])
#    ax = plt.subplot(111)
# create the binned histogram
    databinned, bins = np.histogram(data,nbins,range=[xlow,xhigh])
    centers = bins[:-1] + np.diff(bins)/2 
    if htype=="bar":
        ax.hist(bins[:-1], bins, weights=databinned, histtype='bar')
    else:
        errors = np.sqrt(databinned)
        ax.errorbar(centers, databinned, yerr=errors, fmt='bo', ms=4)

    if Log=='log':
        plt.yscale('log')

    if pdf!=0:
# superimpose function
        ndatatot = np.sum(databinned)
        x = np.linspace(xlow,xhigh,nbins*20)
        binwidth = bins[1]-bins[0]
        pdf_int = np.sum(pdf(x,par))*(x[1]-x[0])
        pdf_norm = pdf(x,par)/pdf_int * ndatatot * binwidth 
        ax.plot(x,pdf_norm,'r-')

    return ax
