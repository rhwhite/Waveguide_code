# Module to FFT and extract wavenumbers 0, 1 and 2 using np.fft
# Written by rhwhite rachel.white@cantab.net
import numpy as np
import xarray as xr
import math

def getFFT(var,norm):
    # Normalize with 1
    var_fft = np.fft.fft(var,axis=-1,norm=norm)
    return(var_fft)


def splitFFT012(FFTin,N):
    return np.absolute(FFTin[...,0])/N,np.absolute(FFTin[...,1])/N,np.absolute(FFTin[...,2])/N

def splitFFT(FFTin,N,maxwn):
    FFTout = []
    for i in range(0,maxwn):
        FFTout.append(np.absolute(FFTin[...,i])/N)
    return FFTout


def get_mean_std(varin):
    varmean = np.mean(varin,axis=0)
    varstd = np.std(varin,axis=0)
    varerror = np.std(varin,axis=0)/math.sqrt(len(varin[:,0]))
    return(varmean,varstd,varerror)

def get_mean_std_year(varin):
    varmean = np.mean(varin,axis=(0,1))
    varstd = np.std(varin,axis=(0,1))

    ## Need to calculate annual means first, and then create standard deviation
    # and calculate standard error on this
    seas_mean = np.mean(varin,axis=1)
    varerror = np.std(seas_mean,axis=(0))/math.sqrt(seas_mean.shape[0])
    return(varmean,varstd,varerror)


def getFFTstats(invar,maxwn):
    try:
        length = len(invar.lon)
    except AttributeError:
        length = len(invar.longitude)
    temp_fft = getFFT(invar,None)

    temp = splitFFT(temp_fft,length,maxwn)

    means_out = []
    std_out = []
    for i in range(0,maxwn):
        means_out.append(get_mean_std(temp[i])[0])
        std_out.append(get_mean_std(temp[i])[1])

    return np.array(means_out), np.array(std_out)

def stat_trans_stats(Daily,lat,lon):

    for iyear in range(0,nyears):

        DJF_daily = Daily[iyear]
        DJF_stat = DJF_daily.mean(axis=0)

        DJF_trans = DJF_daily - DJF_stat

        # make xr again
        DJF_stat = xr.DataArray(DJF_stat, dims=['lat','lon'],
                           coords=[Z_Dec.lat,Z_Dec.lon])
        DJF_trans = xr.DataArray(DJF_trans, dims=['time','lat','lon'],
                           coords=[np.arange(0,DJF_trans.shape[0]),Z_Dec.lat,Z_Dec.lon])

        if iyear == 0:
            DJF_xr_stat = DJF_stat
            DJF_xr_trans = DJF_trans
        else:
            DJF_xr_stat = xr.concat([DJF_xr_stat,DJF_stat],dim='year')
            DJF_xr_trans = xr.concat([DJF_xr_trans,DJF_trans],dim='year')

    # get FFT stats and return
    DJF_stat_stats = getstats(DJF_xr_stat)
    DJF_trans_stats = getstats_year(DJF_xr_trans)

    return DJF_stat_stats, DJF_trans_stats



def getstats(invar,maxwn):
    try:
        length = len(invar.lon)
    except AttributeError:
        length = len(invar.longitude)
    temp_fft = getFFT(invar,None)

    waves = splitFFT(temp_fft,length,maxwn)

    wavestats = {}
    for i in range(0,maxwn):
        wavestats[i] = get_mean_std(waves[i])

    return wavestats 
#wave0_stats, wave1_stats, wave2_stats

def getstats_year(invar,maxwn):
    try:
        length = len(invar.lon)
    except AttributeError:
        length = len(invar.longitude)
    temp_fft = getFFT(invar,None)
    wave0,wave1,wave2 = splitFFT(temp_fft,length,maxwn)
    wave0_stats = get_mean_std_year(wave0)
    wave1_stats = get_mean_std_year(wave1)
    wave2_stats = get_mean_std_year(wave2)

    return wave0_stats, wave1_stats, wave2_stats

