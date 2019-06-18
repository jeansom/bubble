from astropy.io import fits
from astropy.table import Table
import astropy.stats as stats
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import scipy.integrate as integrate

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def getIndices(line, wave_cen):
    cont_ind = False*np.ones(len(wave_cen))
    for c in line['continuum']:
        cont_ind = np.logical_or(cont_ind, (c[0] < wave_cen)*(wave_cen < c[1]))
    line_ind = (wave_cen>line['line'][0])*(wave_cen<line['line'][1])
    entire_ind = (wave_cen > min(wave_cen[np.logical_or(cont_ind, line_ind)]))*(wave_cen < max(wave_cen[np.logical_or(cont_ind, line_ind)]))
    return cont_ind, line_ind, entire_ind

def getLineFit(line, wave_cen, comb_data, indices):
    coeff, var_matrix, p = [ [] for i in range(3) ]
    for i, ind in enumerate(indices):
        cont_ind, line_ind, _ = getIndices(line, wave_cen[ind])
        
        countsflux_conv = comb_data[ind]['NETCOUNTS']/np.maximum(1e-50, comb_data[ind]['FLUX']) # Conversion factor, flux to counts
        countsflux_conv[np.isnan(countsflux_conv)] = 0

        p.append(np.polyfit(wave_cen[ind][cont_ind], comb_data[ind]['FLUX'][cont_ind]*1e14, deg=2))
        c, v = curve_fit(gauss, wave_cen[ind][line_ind], comb_data[ind]['FLUX'][line_ind]*1e14-np.polyval(p[-1], wave_cen[ind][line_ind]), p0=line['p0'])
        coeff.append(c); var_matrix.append(v)
    return coeff, var_matrix, p

def getLineCounts(coeff, p, line, wave_cen, comb_data, indices, par=True):
    flux, counts, counts_errlow, counts_errup = [ np.empty(len(indices)) for i in range(4) ]
    for i, ind in enumerate(indices):
        countsflux_conv = comb_data[ind]['NETCOUNTS']/np.maximum(1e-50, comb_data[ind]['FLUX']) # Conversion factor, flux to counts
        
        _, _, entire_ind = getIndices(line, wave_cen[ind])

        if par: flux[i] = (np.sqrt(2*np.pi)*coeff[i][0]*coeff[i][2]*1e-14)
        else: flux[i] = integrate.simps(comb_data[ind]['FLUX'][entire_ind]-1e-14*np.polyval(p[i], wave_cen[ind][entire_ind]), wave_cen[ind][entire_ind])

        if par: counts[i] = integrate.simps((comb_data[ind]['FLUX'][entire_ind]-1e-14*np.polyval(p[i], wave_cen[ind][entire_ind]))*countsflux_conv[entire_ind], wave_cen[ind][entire_ind])
        else: counts[i] = integrate.simps(gauss(wave_cen[ind], *coeff[i])[entire_ind]*countsflux_conv[entire_ind]/1e14, wave_cen[ind][entire_ind])

        counts_errlow[i], counts_errup[i] = np.abs(stats.poisson_conf_interval(counts[i], interval='frequentist-confidence') - counts[i])
    return flux, counts, counts_errlow, counts_errup