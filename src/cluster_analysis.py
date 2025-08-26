import ROOT
from ROOT import TChain
from array import array
import numpy as npy
import pandas as pd
import sys
import os
import natsort
from pathlib import Path
from multiprocessing import Pool
import re
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit


def process_root_file(file_path):
    ROOT.gInterpreter.GenerateDictionary("vector<vector<int> >", "vector")
    ROOT.gInterpreter.GenerateDictionary("vector<vector float> >", "vector")
    cluster_data_list = []
    try:
        Energy = ROOT.vector('float')()
        RUNID = array('i', [0])
        nfile = ROOT.vector('int')()
        pixels_x = ROOT.vector(ROOT.vector('int'))()
        pixels_y = ROOT.vector(ROOT.vector('int'))()
        pixels_z = ROOT.vector(ROOT.vector('float'))()
        pixels_E = ROOT.vector(ROOT.vector('float'))()
        has_seed = ROOT.vector('int')()
        #excluded = ROOT.vector('int')()
        maxX = ROOT.vector('float')()
        maxY = ROOT.vector('float')()
        minX = ROOT.vector('float')()
        minY = ROOT.vector('float')()
        PosX = ROOT.vector('float')()
        PosY = ROOT.vector('float')()
        meanX = ROOT.vector('float')()
        meanY = ROOT.vector('float')()
        wSTD_X = ROOT.vector('float')()
        wSTD_Y = ROOT.vector('float')()
        fSTD_X = ROOT.vector('float')()
        fSTD_Y = ROOT.vector('float')()
        readout_start = ROOT.vector('float')()
        readout_end = ROOT.vector('float')()

        t4dat_s = TChain("clustersRec")
        t4dat_s.Add(file_path)

        t4dat_s.SetBranchAddress("Energy", Energy)
        t4dat_s.SetBranchAddress("RUNID", RUNID)
        t4dat_s.SetBranchAddress("nfile", nfile)
        t4dat_s.SetBranchAddress("pixels_x", pixels_x)
        t4dat_s.SetBranchAddress("pixels_y", pixels_y)
        t4dat_s.SetBranchAddress("pixels_z", pixels_z)
        t4dat_s.SetBranchAddress("pixels_E", pixels_E)
        t4dat_s.SetBranchAddress("has_seed", has_seed)
        #t4dat_s.SetBranchAddress("excluded", excluded)
        t4dat_s.SetBranchAddress("maxX", maxX)
        t4dat_s.SetBranchAddress("maxY", maxY)
        t4dat_s.SetBranchAddress("minX", minX)
        t4dat_s.SetBranchAddress("minY", minY)
        t4dat_s.SetBranchAddress("PosX", PosX)
        t4dat_s.SetBranchAddress("PosY", PosY)
        t4dat_s.SetBranchAddress("meanX", meanX)
        t4dat_s.SetBranchAddress("meanY", meanY)
        t4dat_s.SetBranchAddress("wSTD_X", wSTD_X)
        t4dat_s.SetBranchAddress("wSTD_Y", wSTD_Y)
        t4dat_s.SetBranchAddress("fSTD_X", fSTD_X)
        t4dat_s.SetBranchAddress("fSTD_Y", fSTD_Y)
        t4dat_s.SetBranchAddress("readout_end", readout_end)
        t4dat_s.SetBranchAddress("readout_start", readout_start)
        clus = ROOT.TH2F("clus", "hXY", 500, -495, 5, 500, -495, 5)

        for i in range(t4dat_s.GetEntries()):
            t4dat_s.GetEntry(i)
            ne = Energy.size()
            runid = RUNID[0]
            for l in range(ne):
                if has_seed[l] == 1:
                    ene = Energy[l]
                    numpix = pixels_x[l].size()
                    x_pix = [pixels_x[l][j] for j in range(numpix)]
                    y_pix = [pixels_y[l][j] for j in range(numpix)]
                    z_pix = [pixels_z[l][j] for j in range(numpix)]
                    pix_e = [pixels_E[l][j] for j in range(numpix)]

                    cluster_data = {
                        'RUNID': runid,
                        'nfile': nfile[l],
                        'has_seed': has_seed[l],
                        #'excluded': excluded[l],
                        'PosX': PosX[l],
                        'PosY': PosY[l],
                        'meanX': meanX[l],
                        'meanY': meanY[l],
                        'xmin': int(minX[l]),
                        'xmax': int(maxX[l]),
                        'ymin': int(minY[l]),
                        'ymax': int(maxY[l]),
                        'npix': numpix,
                        'fpix': numpix / ((maxX[l] - minX[l] + 1) * (maxY[l] - minY[l] + 1)),
                        'energy': ene,
                        'pixels_x': npy.asarray(x_pix),
                        'pixels_y': npy.asarray(y_pix),
                        'pixels_z': npy.asarray(z_pix),
                        'pix_ene': npy.asarray(pix_e),
                        'wsigma_x': wSTD_X[l],
                        'wsigma_y': wSTD_Y[l],
                        'fsigma_x': fSTD_X[l],
                        'fsigma_y': fSTD_Y[l],
                        'readout_start': readout_start[l],
                        'readout_end': readout_end[l]
                    }
                    if ((fSTD_X[l] or fSTD_Y[l]) == -1):
                        cluster_data['fsigma_x'] = 0
                        cluster_data['fsigma_y'] = 0
                
                    cluster_data_list.append(cluster_data)
    except Exception as e:
        print(f"Error while processing {file_path}: {e}")

    return cluster_data_list

def filter_clusters(df,energy_min=0.03,xmin_cut=8,overscan=6152,ymin_cut=1):
    df_has_seed = df[(df['has_seed'] == 1)] #& (df['excluded'] == 0)]  
    clusters_filtered = df_has_seed[
        (df_has_seed['energy'] >= energy_min) &
        (df_has_seed['xmin'] >= xmin_cut) &
        (df_has_seed['xmax'] <= overscan)
    ].copy()

    return clusters_filtered 

def serial_register_events(clusters_df):
    serial_events_x = clusters_df[
        (clusters_df['xmax'] - clusters_df['xmin'] > 2) &
        (clusters_df['fsigma_x'] > 0.5) &
        (clusters_df['fsigma_y'] < 0.2)
    ]
    
    return serial_events_x

def remove_clusters(clusters_df, clusters_to_remove):
    indexes_to_remove = clusters_to_remove.index
    clusters_df = clusters_df.drop(indexes_to_remove)
    return clusters_df

def draw_clusters(clusters_df,figname,im_x=1600,im_y=6400,):
    try:
        hEXY = ROOT.TH2F("hEXY", "hXY", im_y, 0, im_y, im_x, 0, im_x)
        for index, row in clusters_df.iterrows():
            np = len(row['pixels_x'])
            for j in range(np):
                hEXY.Fill(row['pixels_x'][j], row['pixels_y'][j], row['pix_ene'][j])


        c1 = ROOT.TCanvas("c1", "c1",1200,700)
        hEXY.Draw("colz")
        c1.Update()
        c1.Draw()
        input('Press Key to Proceed')
        c1.SaveAs(figname)
        c1.Close()
    except Exception as e:
        print(f"Error while drawing clusters: {e}")

def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    g1 = A1 * npy.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
    g2 = A2 * npy.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
    return g1 + g2

def plot_ene_spectrum(df_list, mode, figname):
    mu1 = 0
    mu2 = 0
    if mode == 'low':
        ene_cut_low = 4.5
        ene_cut_high = 7
        bin_width = 0.1
        plt.xlim(ene_cut_low, ene_cut_high)
        plt.xticks(npy.arange(ene_cut_low,ene_cut_high+bin_width,bin_width))
        plt.title(f'Low Energy ({ene_cut_low}-{ene_cut_high} keV) Spectrum')
        plt.ylabel('Counts')
        plt.xlabel('Energy [keV]')
    elif mode == 'all':
        ene_cut_low = 1
        ene_cut_high = 20
        bin_width = 0.5
        plt.xlim(ene_cut_low, ene_cut_high)
        plt.xticks(npy.arange(ene_cut_low,ene_cut_high+bin_width,bin_width))
        plt.title(f'Low Energy ({ene_cut_low}-{ene_cut_high} keV) Spectrum')
        plt.ylabel('Counts')
        plt.xlabel('Energy [keV]')
    for i, df in enumerate(df_list):
        if df.shape[0] > 1:
            counts, bins = npy.histogram(df[(df['energy'] >= ene_cut_low) & (df['energy'] <= ene_cut_high)]['energy'], bins=npy.arange(ene_cut_low, ene_cut_high + bin_width, bin_width))
            bin_centers = (bins[:-1] + bins[1:]) / 2
            plt.bar(bin_centers, counts, width=bin_width, alpha=0.6, label="Histogram")
            if mode == 'low':
                initial_guesses = [200, 5.65, 0.1, 20, 6.2, 0.15]
                bounds = ([0, 5, 0, 0, 6.0, 0], [npy.inf, 6.6, 0.2, npy.inf, 6.8, 0.2])
                popt, pcov = curve_fit(double_gaussian, bin_centers, counts, p0=initial_guesses, bounds=bounds)
                A1, mu1, sigma1, A2, mu2, sigma2 = popt
                perr = npy.sqrt(npy.diag(pcov))
                mu1_err = perr[1]
                mu2_err = perr[4]
                
                x_fit = npy.linspace(ene_cut_low, ene_cut_high, 1000)
                g1 = A1 * npy.exp(-0.5 * ((x_fit - mu1) / sigma1) ** 2)
                g2 = A2 * npy.exp(-0.5 * ((x_fit - mu2) / sigma2) ** 2)
                print(f"Gaussian 1: mu1={mu1:.2f} +/- {mu1_err:.2f}")
                print(f"Gaussian 2: mu2={mu2:.2f} +/- {mu2_err:.2f}")
                plt.plot(x_fit, g1, 'g--', label=f"Gaussian 1: mu1={mu1:.2f} +/- {mu1_err:.2f}")
                plt.plot(x_fit, g2, 'b--', label=f"Gaussian 2: mu2={mu2:.2f} +/- {mu2_err:.2f}")
                
                
                
    plt.legend()
    plt.savefig(figname, dpi=300)
    plt.show()
    return mu1, mu2

def plot_energy_sigma(dataframes, labels, colors, figname='energysigmaxy.png'):
    
    for df, label, col in zip(dataframes, labels, colors):
        sigma_xy = ((df['fsigma_x']**2 + df['fsigma_y']**2) / 2)**0.5
        plt.scatter(df['energy'], sigma_xy, label=label, color = col)
    
    plt.xticks(npy.arange(4.5,7.1,0.1))
    plt.xlim(4.5,7)
    plt.ylim(0,2)
    plt.axhspan(0,0.5,xmin=(5.4-4.5)/(7-4.5), xmax=(6.8-4.5)/(7-4.5),color='red',alpha=0.3,label='Front Events')
    plt.axhspan(0.8,1.3,color='green',alpha=0.3,label='Back Events')
    plt.title(r'Energy ${\sigma_{xy}}$ Distribution')
    plt.ylabel(r'${\sigma_{xy}}$ [pixels]')
    plt.xlabel('Energy [keV]')
    plt.grid()
    plt.legend()
    plt.savefig(figname,dpi=300)
    plt.show()


def plot_x_energy_sigma(dataframes, labels, colors, figname='energypositionxsigmaxy.png'):
    
    for df, label, col in zip(dataframes, labels, colors):
        df = df[(df['energy'] < 7) & (df['energy'] > 4.5)]
        sigma_xy = ((df['fsigma_x']**2 + df['fsigma_y']**2) / 2)**0.5
        plt.scatter(df['meanX'], sigma_xy, label=label, color = col)
    
    plt.xticks(npy.arange(0,6600,200))
    plt.xlim(0,6400)
    plt.ylim(0,2)
    plt.axhspan(0,0.5,xmin=(3800)/(6400), xmax=(6144)/(6400),color='red',alpha=0.3,label='Front Events')
    plt.axhspan(0.8,1.3,xmin=(2400)/(6400), xmax=(4800)/(6400),color='green',alpha=0.3,label='Back Events')
    plt.title(r'Energy ${\sigma_{xy}}$ Distribution')
    plt.ylabel(r'${\sigma_{xy}}$ [pixels]')
    plt.xlabel('MeanX [Column Number]')
    plt.grid()
    plt.legend()
    plt.savefig(figname,dpi=300)
    plt.show()
   
def plot_pixels_ene(df_list, labels, col, figname):
    try:
        for df, label,col in zip(df_list, labels, colors):
            pix_vals = []
            for index, row in df.iterrows():
                np = len(row['pix_ene'])
                for j in range(np):
                    pix_vals.append(row['pix_ene'][j])
                     
            pix_vals = npy.asarray(pix_vals)
            max_pix_val = round(max(pix_vals),2)
            plt.xlabel('Energy of pixels [keV]')
            plt.ylabel('Counts')
            plt.title('Pix Vals of CCDs')
            plt.hist(pix_vals,bins=npy.arange(1,50.2,0.2),color=col,alpha=0.3,label='{} [Max pixel value: {} keV]'.format(label,max_pix_val))
        plt.yscale('log')
        plt.legend()
        plt.savefig(figname,dpi=300)
        plt.show()
    except Exception as e:
        print(f"Error while drawing clusters: {e}")
    
def make_composite(dataframes, label, canvas_size = 20, sigma_xy_cut = 'front', lower_condition = 4.5, higher_condition = 7):
    try:
        hEXY = ROOT.TH2F("hEXY", "hXY", canvas_size*2, -canvas_size, canvas_size, canvas_size*2, -canvas_size, canvas_size)
        for df, label in zip(dataframes, labels):
            clusters_df = df[(df['energy']>lower_condition) & (df['energy']<higher_condition)]
            clusters_df_copy = clusters_df.copy()
            clusters_df_copy['sigma_xy'] = ((clusters_df_copy['fsigma_x']**2 + clusters_df_copy['fsigma_y']**2) / 2)**0.5
            if (sigma_xy_cut == 'front'):
                filtered_df = clusters_df_copy[(clusters_df_copy['sigma_xy']<0.5) & (clusters_df_copy['meanX']>3800)]
            elif (sigma_xy_cut == 'back'):
                filtered_df = clusters_df_copy[(clusters_df_copy['sigma_xy']>0.8) & (clusters_df_copy['sigma_xy']<1.3) & (clusters_df_copy['meanX']<4800) & (clusters_df_copy['meanX']>2400)]
            hEXY.Reset()
            num_clus = filtered_df.shape[0]
            for index,row in filtered_df.iterrows():
                x = row['pixels_x'] - row['meanX'] 
                y = row['pixels_y'] - row['meanY'] 
                np = len(row['pixels_x'])
                for j in range(np):
                    hEXY.Fill(x[j], y[j], row['pix_ene'][j])    
            c1 = ROOT.TCanvas("c1", "c1",1200,700)
            center = ROOT.TPad("center", "center", 0, 0, 0.6, 0.6)
            center.Draw()
            right = ROOT.TPad("right", "right", 0.55, 0, 1, 0.6);
            right.Draw()
            top = ROOT.TPad("top", "top", 0, 0.55, 0.6, 1);
            top.Draw()
            
            center.cd()
            hEXY.Draw("colz")
            top.cd()
            projh2x = hEXY.ProjectionX()
            projh2x.Draw("bar")
            for bin in range(1, projh2x.GetNbinsX() + 1):
                if projh2x.GetBinContent(bin) > 0:
                    first_bin = bin
                    break
            for bin in range(1, projh2x.GetNbinsX() + 1):
                if projh2x.GetBinContent(bin) > 0:
                    last_bin = bin
            xmin = projh2x.GetXaxis().GetBinCenter(first_bin) - 0.5
            xmax = projh2x.GetXaxis().GetBinCenter(last_bin) + 0.5
            fitxr = ROOT.TF1("fitxr", "gaus")
            projh2x.Fit("fitxr", "LQIE0", "", -0.5, xmax)
            r_sigma = fitxr.GetParameter(2)
            r_sigma_err = fitxr.GetParError(2)
            fitxl = ROOT.TF1("fitxl", "gaus")
            projh2x.Fit("fitxl", "LQIE", "", xmin, 0.5)
            l_sigma = fitxl.GetParameter(2)
            l_sigma_err = fitxl.GetParError(2)
            fitxr.SetLineColor(ROOT.kGreen)
            fitxr.SetLineWidth(2)
            fitxr.Draw("same")
            right.cd()
            projh2y = hEXY.ProjectionY()
            projh2y.Draw("hbar")
            for bin in range(1, projh2y.GetNbinsX() + 1):
                if projh2y.GetBinContent(bin) > 0:
                    first_bin = bin
                    break
            for bin in range(1, projh2y.GetNbinsX() + 1):
                if projh2y.GetBinContent(bin) > 0:
                    last_bin = bin
            ymin = projh2y.GetXaxis().GetBinCenter(first_bin) - 0.5
            ymax = projh2y.GetXaxis().GetBinCenter(last_bin) + 0.5
            fityt = ROOT.TF1("fityt", "gaus")
            projh2y.Fit("fityt", "LQIE0", "", -0.5, ymax)
            t_sigma = fityt.GetParameter(2)
            t_sigma_err = fityt.GetParError(2)
            fityb = ROOT.TF1("fityb", "gaus")
            projh2y.Fit("fityb", "LQIE0", "", ymin, 0.5)
            b_sigma = fityb.GetParameter(2)
            b_sigma_err = fityb.GetParError(2)
            fityt.SetLineColor(ROOT.kGreen)
            fityt.SetLineWidth(2)
            fityt.Draw("same")
            
            x_sigma = (l_sigma + r_sigma)/2
            x_sigma_err = (npy.square(l_sigma_err)+npy.square(r_sigma_err))**0.5/2
            y_sigma = (b_sigma + t_sigma)/2
            y_sigma_err = (npy.square(b_sigma_err)+npy.square(t_sigma_err))**0.5/2
            sigma = (npy.square(x_sigma)/2 + npy.square(y_sigma)/2)**0.5
            sigma_err = (1/(2*sigma))*(((x_sigma*x_sigma_err)**2+(y_sigma*y_sigma_err)**2)**0.5)
            print("Sigma: {:.2f} +/- {:.2f}".format(sigma,sigma_err))

            
            c1.cd()
            t = ROOT.TLatex()
            t.SetTextFont(42)
            t.SetTextSize(0.02)
            t.DrawLatex(0.6, 0.88, "Composite Cluster of {}keV < {} Clusters < {}keV and its X and Y projections".format(lower_condition, label, higher_condition))
            t.DrawLatex(0.6, 0.84, "Left Sigma: {:.2f} +/- {:.2f}, Right Sigma: {:.2f} +/- {:.2f}".format(l_sigma, l_sigma_err, r_sigma, r_sigma_err))
            t.DrawLatex(0.6, 0.80, "Assymetry X: {:.2f} +/- {:.2f}".format(r_sigma-l_sigma, npy.sqrt(npy.power(l_sigma_err,2)/l_sigma + npy.power(r_sigma_err,2)/r_sigma)))
            t.DrawLatex(0.6, 0.72, "Top Sigma: {:.2f} +/- {:.2f}, Bottom Sigma: {:.2f} +/- {:.2f}".format(t_sigma, t_sigma_err, b_sigma, b_sigma_err))
            t.DrawLatex(0.6, 0.68, "Assymetry Y: {:.2f} +/- {:.2f}".format(b_sigma-t_sigma, npy.sqrt(npy.power(t_sigma_err,2)/t_sigma + npy.power(b_sigma_err,2)/b_sigma)))
            t.DrawLatex(0.6, 0.62, "Sigma: {:.2f} +/- {:.2f}".format(sigma,sigma_err))
            
            c1.Update()
            c1.Draw()
            input('Press Key to Proceed')
            c1.SaveAs('Composite_cluster__{}_{}.png'.format(sigma_xy_cut,label))
            c1.Close()
    except Exception as e:
        print(f"Error while drawing clusters: {e}")

def plot_cti_front_events(df, labels, col, mu, figname):
    try:
        sigma_xy = ((df['fsigma_x']**2 + df['fsigma_y']**2) / 2)**0.5
        df = df[(df['energy'] < 7) & (df['energy'] > 4.5) & (df['meanX'] > 3800) & (sigma_xy < 0.5)]
        
        mean_x = df['meanX'].values
        pix_vals = npy.concatenate(df['pix_ene'].values)
        max_pix_vals = df['pix_ene'].apply(max).values
        
        fig, axs = plt.subplots(1, 2, figsize=(8, 6))
        scatter = axs[1].scatter(mean_x, max_pix_vals, color=col, alpha=0.8, label=labels)
        axs[1].set_xlabel('Mean X')
        axs[1].set_ylabel('Max PixVals [keV]')
        axs[1].set_title('Max PixVals vs Mean X')
        axs[1].hlines(mu, 3800, 6400, colors='red', linestyles='--', label='Threshold')
        axs[1].set_xticks(npy.arange(3800, 6600, 200))
        axs[1].set_ylim(4.5, 7)
        axs[1].set_xlim(3800, 6400)
        axs[0].hist(max_pix_vals, bins=npy.arange(4.5, 7, 0.1), orientation='horizontal', color=col, alpha=0.8, edgecolor='black')
        axs[0].hlines(mu, xmin=0, xmax=max(npy.histogram(max_pix_vals, bins=npy.arange(4.5, 7, 0.1))[0]), colors='red', linestyles='--', label='Threshold')
        axs[0].set_xlabel('Counts')
        axs[0].set_ylabel('Max PixVals [keV]')
        axs[0].set_title('Projection of Max PixVals')
        axs[0].set_ylim(4.5, 7)
        plt.title('CTI of Fe-55 front events')
        plt.savefig(figname,dpi=300)
        plt.show()
    except Exception as e:
        print(f"Error while drawing clusters: {e}")
    
def print_cluster_statistics(text, clusters_df, columns_to_print, print_all_info=False):
    if print_all_info == True and not clusters_df.empty:
        clusters_df_copy = clusters_df.copy()
        clusters_df_copy['readout_start'] = clusters_df_copy['readout_start'].apply(lambda x: datetime.fromtimestamp(x).strftime('%Y-%m-%d %H:%M:%S'))
        print("Clusters in {} :\n".format(text), clusters_df_copy[columns_to_print])
    print("Number of clusters in {}:".format(text), clusters_df.shape[0])
    
def inspect_clusters(clusters_df):    
    for index, row in clusters_df.iterrows():
        x_range = max(row['pixels_x']) - min(row['pixels_x']) + 1
        y_range = max(row['pixels_y']) - min(row['pixels_y']) + 1
        hEXY = npy.zeros((x_range, y_range))
        x = row['pixels_x'] - min(row['pixels_x'])
        y = row['pixels_y'] - min(row['pixels_y'])
        hEXY[x, y] += row['pix_ene']

        fig, ax = plt.subplots()
        X, Y = npy.meshgrid(npy.arange(x_range+1), npy.arange(y_range+1))
        im = ax.pcolormesh(X, Y, hEXY.T, cmap='jet')
        x_sum = npy.sum(hEXY, axis=0)  
        y_sum = npy.sum(hEXY, axis=1)
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="1%", pad="1%")
        cbar = fig.colorbar(im, cax=cax, ticks=[0, npy.max(hEXY)], orientation="vertical")
        cbar.ax.tick_params(labelsize=5) 
        axtop = divider.append_axes("bottom", size=1.7, pad=1.0, sharex=ax)
        axtop.step(npy.arange(len(y_sum))+0.5,y_sum,where='mid')   
        axtop.set_ylabel("Energy [keV]")
        axtop.set_xlabel("Column Number")
        axtop.set_title('Column-wise sum',fontweight="bold")
        axright = divider.append_axes("right", size=2.0, pad=1.0, sharey=ax)
        axright.step(x_sum,npy.arange(len(x_sum)),where='pre')
        axright.set_xlabel("Energy [keV]")
        axright.set_ylabel("Row Number")
        axright.set_title('Row-wise sum',fontweight="bold")
        plt.tight_layout()
        ax.set_aspect('auto', adjustable='box', anchor='C')
        ax.set_title('Cluster (Energy: {} keV)'.format(round(row['energy'],2)), fontweight="bold")
        ax.set_xlabel("Column Number")
        ax.set_ylabel("Row Number")
        print_cluster_statistics("Cluster inspection", pd.DataFrame([row]), columns_to_print, print_all_info=True)
        plt.show()


if __name__ == "__main__":
    try:
        args = sys.argv
    except IndexError:
        raise SystemExit(f"Usage: {sys.argv[0]} <Directory or File> <ccd_name_ext_name>")

    input_path = args[1]
    ccd_name_ext_name = args[2]
    if os.path.isdir(input_path):
        directory = input_path
        root_files = []
        for root, _, filenames in os.walk(directory):
            filenames = natsort.natsorted(filenames)
            file_paths = [os.path.join(root, filename) for filename in filenames if filename.endswith('.root')]
            root_files.extend(file_paths)
    else:
        root_files = [input_path]

    # Use multiprocessing to process ROOT files concurrently
    for root_file in root_files:
        print(f"Processing file: {root_file}")
        cluster_data_list = process_root_file(root_file)
        if not cluster_data_list:
            print(f"Skipping {root_file} due to empty data.")
            continue
        df = pd.DataFrame.from_records(cluster_data_list)
        df = df.drop(columns=['RUNID'])
        
        ccd_df = df
        
        # Define columns to print
        columns_to_print = ['energy', 'nfile', 'xmin', 'xmax', 'ymin', 'ymax', 'npix', 'fpix', 'fsigma_x', 'fsigma_y', 'readout_start']
        
        clusters = filter_clusters(ccd_df)
        
        serial_events_x = serial_register_events(clusters)
        clusters = remove_clusters(clusters,serial_events_x)
        
        # Print cluster statistics for different datasets and conditions
        print_cluster_statistics('Clusters CCD [4.5-7keV]', clusters[(clusters['energy'] > 4.5) & (clusters['energy'] < 7)], columns_to_print, print_all_info=False)
        print_cluster_statistics('Clusters CCD ', clusters, columns_to_print, print_all_info=False)
        
        # Prepare data for energy vs. sigma plots
        colors = ['blue']
        dataframes = [clusters]
        labels = [f'{ccd_name_ext_name}']  
        
        # Draw 2D histograms for all clusters
        draw_clusters(clusters, f'2dhist_{ccd_name_ext_name}.png')

        # Plot energy spectra
        mu1, mu2 = plot_ene_spectrum([clusters], 'low', f'energyspectrum_{ccd_name_ext_name}.png')

        # Plot energy vs. sigma_xy
        plot_energy_sigma(dataframes, labels, colors, f'energy_sigma_{ccd_name_ext_name}.png')
        
        # Plot position vs. sigma_xy
        plot_x_energy_sigma(dataframes, labels, colors, f'xpos_energy_sigma_{ccd_name_ext_name}.png')
        
        #Plot CTI of front events
        plot_cti_front_events(clusters, labels, colors, mu1, f'ctifrontevents_{ccd_name_ext_name}.png')
        
        #Plot a composite cluster made from all the clusters
        make_composite(dataframes,labels,sigma_xy_cut='front')
        make_composite(dataframes,labels,sigma_xy_cut='back')

        #Plot pixels_ene distribution
        plot_pixels_ene(dataframes, labels, colors, f'pixvals_{ccd_name_ext_name}.png')
    
    
