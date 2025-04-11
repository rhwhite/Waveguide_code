# import
import os, errno
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
#import Ngl
import math
import timeit

from copy import copy

from scipy import stats,signal
from scipy.interpolate import interp1d,CubicSpline
import scipy.fftpack as fftpack
from scipy.signal import argrelextrema
from sklearn.linear_model import LinearRegression
import cartopy.crs as ccrs
import cftime


# My packages
import rhwhitepackages3
from rhwhitepackages3.readwrite import *
from rhwhitepackages3.stats import regressmaps
#from rhwhitepackages3.CESMmanip import *
from rhwhitepackages3.wavenumber import *
#from rhwhitepackages3.CESMconst import *
from rhwhitepackages3.processing import *
from rhwhitepackages3.clim_tools import *
from rhwhitepackages3.physconst import *

def calc_wavepacket_3D_tur(datain,len_min,len_max):
    # Subtract zonal mean
    test_data = datain - datain.mean(dim='longitude')
    lats = test_data.latitude
    lons = test_data.longitude
    times = test_data.time

    wavepacket = xr.DataArray(np.zeros([len(times),len(lats),len(lons)]),
                              coords={'time':times,'latitude':lats,'longitude':lons},
                              dims = ('time','latitude','longitude'))

    templats = lats.sel(latitude=slice(85,20))
    ilatstart = np.where(lats == templats[0])[0]
    ilatend = np.where(lats == templats[-1])[0]

    for latsel in np.arange(ilatstart,ilatend+1):
        # Find wavenumbers for km values at given latitude
        # include 0.001 factor to convert rearth from m to km
        circum = (2.0 * np.pi * (0.001 * rearth * np.cos(np.deg2rad(lats.isel(latitude=latsel)))))
        s2 = (circum/len_min).values
        s1 = (circum/len_max).values

        indata = test_data.isel(latitude=latsel)
        inlat = indata.latitude.values

        # Tukey transform the data
        nlons = len(indata.longitude)
        tukey_transform, std_transform,t = fourier_Tukey_hilow_3D(indata,nlons,s1,s2,ndegs=360)
        # Now hilbert transfrom (For now, ignoring the semi-geostrophic filter)
        hilbert_transform = fourier_Hilbert_3D(tukey_transform,nlons)

        wavepacket[:,latsel,:] = 2.0 *np.abs(hilbert_transform)[0]

    return(wavepacket)

def calc_wavepacket_3D_std(datain,len_min,len_max):
    # Subtract zonal mean
    test_data = datain - datain.mean(dim='longitude')
    lats = test_data.latitude
    lons = test_data.longitude
    times = test_data.time

    wavepacket = xr.DataArray(np.zeros([len(times),len(lats),len(lons)]),
                              coords={'time':times,'latitude':lats,'longitude':lons},
                              dims = ('time','latitude','longitude'))

    templats = lats.sel(latitude=slice(85,20))
    ilatstart = np.where(lats == templats[0])[0]
    ilatend = np.where(lats == templats[-1])[0]

    for latsel in np.arange(ilatstart,ilatend+1):
        # Find wavenumbers for km values at given latitude
        # include 0.001 factor to convert rearth from m to km
        circum = (2.0 * np.pi * (0.001 * rearth * np.cos(np.deg2rad(lats.isel(latitude=latsel)))))
        s2 = (circum/len_min).values
        s1 = (circum/len_max).values

        indata = test_data.isel(latitude=latsel)
        inlat = indata.latitude.values

        # Tukey transform the data
        nlons = len(indata.longitude)
        tukey_transform, std_transform,t = fourier_Tukey_hilow_3D(indata,nlons,s1,s2,ndegs=360)
        # Now hilbert transfrom (For now, ignoring the semi-geostrophic filter)
        hilbert_transform = fourier_Hilbert_3D(std_transform,nlons)

        wavepacket[:,latsel,:] = 2.0 *np.abs(hilbert_transform)[0]

    return(wavepacket)

def calc_wavepacket_3D_fixed_std(datain,wn_min,wn_max,hemi='NH'):
    # Subtract zonal mean
    test_data = datain - datain.mean(dim='longitude')
    lats = test_data.latitude
    lons = test_data.longitude
    times = test_data.time

    wavepacket = xr.DataArray(np.zeros([len(times),len(lats),len(lons)]),
                              coords={'time':times,'latitude':lats,'longitude':lons},
                              dims = ('time','latitude','longitude'))

    
    if hemi == 'NH':
        templats = lats.sel(latitude=slice(85,20))
    elif hemi == 'SH':
        templats = lats.sel(latitude=slice(-20,-85))
 
    ilatstart = np.where(lats == templats[0])[0]
    ilatend = np.where(lats == templats[-1])[0]

    for latsel in np.arange(ilatstart,ilatend+1):
        # Find wavenumbers for km values at given latitude
        # include 0.001 factor to convert rearth from m to km
        circum = (2.0 * np.pi * (0.001 * rearth * np.cos(np.deg2rad(lats.isel(latitude=latsel)))))
        s2 = wn_max
        s1 = wn_min

        indata = test_data.isel(latitude=latsel)
        inlat = indata.latitude.values

        # Tukey transform the data
        nlons = len(indata.longitude)
        tukey_transform, std_transform,t = fourier_Tukey_hilow_3D(indata.values,nlons,s1,s2,ndegs=360)
        # Now hilbert transfrom (For now, ignoring the semi-geostrophic filter)
        hilbert_transform = fourier_Hilbert_3D(std_transform,nlons)

        wavepacket[:,latsel,:] = 2.0 *np.abs(hilbert_transform)[0]

    return(wavepacket)


def calc_wavepacket(datain,len_min,len_max):
    # Subtract zonal mean
    test_data = datain - datain.mean(dim='longitude')
    lats = test_data.latitude
    lons = test_data.longitude
    times = test_data.time

    wavepacket = xr.DataArray(np.zeros([len(lats),len(lons)]),
                              coords={'latitude':lats,'longitude':lons},
                              dims = ('latitude','longitude'))

    templats = lats.sel(latitude=slice(85,20))
    ilatstart = np.where(lats == templats[0])[0]
    ilatend = np.where(lats == templats[-1])[0]

    for latsel in np.arange(ilatstart,ilatend+1):
        # Find wavenumbers for km values at given latitude
        # include 0.001 factor to convert rearth from m to km
        circum = (2.0 * np.pi * (0.001 * rearth * np.cos(np.deg2rad(lats.isel(latitude=latsel)))))
        s2 = (circum/len_min).values
        s1 = (circum/len_max).values

        indata = test_data.isel(latitude=latsel)
        inlat = indata.latitude.values

        # Tukey transform the data
        nlons = len(indata.longitude)
        tukey_transform, std_transform,t = fourier_Tukey_hilow(indata,nlons,s1,s2,ndegs=360)
        # Now hilbert transfrom (For now, ignoring the semi-geostrophic filter)
        hilbert_transform = fourier_Hilbert(tukey_transform,nlons)

        wavepacket[latsel,:] = 2.0 *np.abs(hilbert_transform)[0]

    return(wavepacket)
