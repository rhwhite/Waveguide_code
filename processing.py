# Module for performing various processing of data
# Written by rhwhite rachel.white@cantab.net
import numpy as np
import xarray as xr
from scipy import stats,signal
import scipy.fftpack as fftpack
from scipy.interpolate import interp1d,CubicSpline

# Lanczos function weight calculation. confirmed it matches that from NCL when normalization added.
def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma

    # Needed to add a normalization factor to make the weights add up to 1.
    norm = np.sum(w)
    w = w/norm
    return w[1:-1]

def Lanczos_filter(datain,nwghts,fca):
    # Calculate Lanczos weights
    lp_wghts = low_pass_weights(nwghts,fca)

    timefilt = np.ndarray(datain.shape)
    try:
        lons = datain.lon
    except AttributeError:
        lons = datain.longitude

    try:
        lats = datain.lat
    except AttributeError:
        lats = datain.latitude

    nlons = len(lons)
    nlats = len(lats)
    ntimes = len(datain.time)

    temp = xr.DataArray(datain.values,coords={'time':datain.time,
                                             'longitude':lons.values,'latitude':lats.values},
                                      dims = ('time','latitude','longitude'))

    # convolve data with weights
    for ilat in np.arange(0,nlats):
        for ilon in np.arange(0,nlons):
            timefilt[:,ilat,ilon] = np.convolve(datain[:,ilat,ilon], lp_wghts, 'same')

    timefilt[0:int((nwghts-1)/2),:,:] = np.nan
    timefilt[int(ntimes - (nwghts-1)/2):ntimes,:,:] = np.nan

    # save to xarray
    xr.DataArray(np.random.randn(2, 3), coords={'x': ['a', 'b']}, dims=('x', 'y'))

    xrtimefilt = xr.DataArray(timefilt,coords={'time':datain.time,
                                             'longitude':lons.values,'latitude':lats.values},
                                      dims = ('time','latitude','longitude'))
    return(xrtimefilt)


def fourier_Tukey(indata,nlons,peak_freq,ndegs=360):
    X_fft = fftpack.fft(indata)

    f_s = nlons
    freqs = fftpack.fftfreq(len(indata)) * f_s
    t = np.linspace(0, ndegs, f_s, endpoint=False)

    filt_fft = X_fft.copy()
    filt_fft[np.abs(freqs) > peak_freq] = 0
    filtered_sig = fftpack.ifft(filt_fft)

    # create Tukey window to smooth the wavenumbers removed (so no exact cutoff at k=5, 
    #which will change at different latitudes)
    # Window is 2 wavenumbers more than the peak, but multiplied by 2 because the Tukey window is symmetric
    M = (peak_freq + 2)*2
    alpha = 0.3 # co-sine weighting covers 30% of the window
    tukeyWin = signal.tukey(M, alpha=0.3, sym=True)[int(M/2):M]

    turfilt_fft = X_fft.copy()
    n = len(turfilt_fft)
    turfilt_fft[0:int(M/2)] = turfilt_fft[0:int(M/2)]*tukeyWin
    turfilt_fft[int(M/2):n-int(M/2)] = 0
    turfilt_fft[n-int(M/2):n] = turfilt_fft[n-int(M/2):n]*tukeyWin[::-1]
    tur_filtered_sig = fftpack.ifft(turfilt_fft)

    return(tur_filtered_sig,filtered_sig,t)

def fourier_Tukey_hilow_3D(indata,nlons,s1,s2,ndegs=360):
    try:
        ntimes = len(indata.times)
    except AttributeError:
        ntimes = 1

    X_fft = fftpack.fft(indata)

    f_s = nlons
    freqs = fftpack.fftfreq(nlons) * f_s
    t = np.linspace(0, ndegs, f_s, endpoint=False)

    filt_fft = X_fft.copy()
    filt_fft[:,np.abs(freqs) > s2] = 0
    filt_fft[:,np.abs(freqs) < s1] = 0
    filtered_sig = fftpack.ifft(filt_fft)

    # create Tukey window to smooth the wavenumbers removed (so no exact cutoff at any particular wavenumber, 
    # which will change at different latitudes)
    # Window is 2 wavenumbers more than the peak, but multiplied by 2 because the Tukey window is symmetric
    limit1 = np.amax([0,s1 - 1.5])
    limit2 = s2 + 1.5

    diff = limit2-limit1
    M = 100

    tukeymax = int(np.amax([limit1,limit2])) + 10

    tukeydata = signal.tukey(M, alpha=0.3, sym=True)

    # Add zeros at front and back of window for interpolation
    tukeydata = np.insert(tukeydata,0,np.zeros(40))
    tukeydata = np.append(tukeydata,np.zeros(40))

    # wavenumber values for Tukey window with added zeros
    xs = np.arange(limit1-40.0*diff/M,limit2+ 39.999*diff/M,diff/M)

    # model tukey window with cubic splines
    wn_cs = np.linspace(0,tukeymax-1,tukeymax)
    CS = CubicSpline(xs,tukeydata)

    tukeyWin = CS(wn_cs)
    turfilt_fft = X_fft.copy()
    n = len(turfilt_fft[0,:])
    turfilt_fft[:,0:tukeymax] = turfilt_fft[:,0:tukeymax]*tukeyWin
    turfilt_fft[:,tukeymax:n-tukeymax] = 0
    turfilt_fft[:,n-tukeymax:n] = turfilt_fft[:,n-tukeymax:n]*tukeyWin[::-1]
    tur_filtered_sig = fftpack.ifft(turfilt_fft)

    return(tur_filtered_sig,filtered_sig,t)


def fourier_Tukey_hilow(indata,nlons,s1,s2,ndegs=360):
    X_fft = fftpack.fft(indata)

    f_s = nlons
    freqs = fftpack.fftfreq(len(indata)) * f_s
    t = np.linspace(0, ndegs, f_s, endpoint=False)

    filt_fft = X_fft.copy()
    filt_fft[np.abs(freqs) > s2] = 0
    filt_fft[np.abs(freqs) < s1] = 0
    filtered_sig = fftpack.ifft(filt_fft)

    # create Tukey window to smooth the wavenumbers removed (so no exact cutoff at any particular wavenumber, 
    # which will change at different latitudes)
    # Window is 2 wavenumbers more than the peak, but multiplied by 2 because the Tukey window is symmetric
    limit1 = np.amax([0,s1 - 1.5])
    limit2 = s2 + 1.5

    diff = limit2-limit1
    M = 100

    tukeymax = int(np.amax([limit1,limit2])) + 10
    tukeydata = signal.tukey(M, alpha=0.3, sym=True)

    # Add zeros at front and back of window for interpolation
    tukeydata = np.insert(tukeydata,0,np.zeros(40))
    tukeydata = np.append(tukeydata,np.zeros(40))

    # wavenumber values for Tukey window with added zeros
    xs = np.arange(limit1-40.0*diff/M,limit2+ 39.999*diff/M,diff/M)

    # model tukey window with cubic splines
    wn_cs = np.linspace(0,tukeymax-1,tukeymax)
    CS = CubicSpline(xs,tukeydata)

    tukeyWin = CS(wn_cs)
    turfilt_fft = X_fft.copy()
    n = len(turfilt_fft)
    turfilt_fft[0:tukeymax] = turfilt_fft[0:tukeymax]*tukeyWin
    turfilt_fft[tukeymax:n-tukeymax] = 0
    turfilt_fft[n-tukeymax:n] = turfilt_fft[n-tukeymax:n]*tukeyWin[::-1]
    tur_filtered_sig = fftpack.ifft(turfilt_fft)

    return(tur_filtered_sig,filtered_sig,t)


def fourier_Hilbert_3D(indata,nlons,peak_freq=0,ndegs=360):
    X_fft = fftpack.fft(indata)
    f_s = nlons
    freqs = fftpack.fftfreq(nlons) * f_s
    t = np.linspace(0, ndegs, f_s, endpoint=False)

    filt_fft = X_fft.copy()
    if peak_freq != 0:
        filt_fft[:,np.abs(freqs) > peak_freq] = 0
    
    # Hilbert transform: set negatives = 0
    filt_fft[:,freqs<0] = 0
    filtered_sig = fftpack.ifft(filt_fft)

    return(filtered_sig,t)

def fourier_Hilbert(indata,nlons,peak_freq=0,ndegs=360):
    X_fft = fftpack.fft(indata)

    f_s = nlons
    freqs = fftpack.fftfreq(len(indata)) * f_s
    t = np.linspace(0, ndegs, f_s, endpoint=False)

    filt_fft = X_fft.copy()
    if peak_freq != 0:
        filt_fft[np.abs(freqs) > peak_freq] = 0
    
    # Hilbert transform: set negatives = 0
    filt_fft[freqs<0] = 0
    filtered_sig = fftpack.ifft(filt_fft)

    return(filtered_sig,t)

