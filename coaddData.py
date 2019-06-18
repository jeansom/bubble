from astropy.io import fits
import numpy as np

def coaddData(exp_pre, keys, wave_bins):
    exp_hdul = [ fits.open(pre+'x1d.fits') for pre in exp_pre ]
    exp_data = [ exp_hdul[i][1].data for i in range(len(exp_pre)) ]
    # DQ cuts
    DQ_data = [{} for i in range(len(exp_pre))]
    for i in range(len(exp_pre)):
        keepinds = np.logical_or(exp_data[i]['DQ'] == 0, exp_data[i]['DQ'] == 1)
        for key in keys:
            if key == 'EXPTIME': DQ_data[i][key] = np.ones(len(keepinds[keepinds].flatten()))*exp_data[i][key]
            elif key == 'NETCOUNTS': DQ_data[i][key] = np.array(exp_data[i]['NET'][keepinds])*exp_data[i]['EXPTIME']
            else: DQ_data[i][key] = np.array(exp_data[i][key][keepinds])

    # Rebin data into same wavelength bins and average (with exposure weighting)
    comb_data = {}
    binned_data = [{} for i in range(len(exp_pre))]

    exp_all = np.array([ exp_data[i]['EXPTIME'] for i in range(len(exp_pre)) ]).flatten()

    for i in range(len(exp_pre)):
        for key in keys:
            binned_data[i][key] = np.histogram(DQ_data[i]['WAVELENGTH'], weights=DQ_data[i][key], bins=wave_bins)[0]

    for key in keys:
        data_key = [ binned_data[i][key] for i in range(len(exp_pre)) ]
        comb_data[key] = np.average(data_key, axis=0, weights=exp_all)
    return exp_hdul, exp_data, DQ_data, binned_data, comb_data