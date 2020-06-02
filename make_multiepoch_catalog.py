import numpy as np
import sys
sys.path.insert(0,'/home/jill/python/')
import AstroTools as at
import CRRPphot as crrp
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from astropy.time import Time

channels = np.array(['I1', 'I2'])
n_filters = 2
target = 'NGC6402'
n_epochs = 12

os.chdir('{}'.format(target))


# Setup array structures using input star list
input_list = crrp.dao.read_mag('median.off')

ids = input_list['id']
sort_ids = ids.argsort()
n_stars = len(input_list['id'])

print(n_stars)

multiepoch_mags = np.zeros((n_stars, n_epochs, n_filters))
multiepoch_errs = np.zeros((n_stars, n_epochs, n_filters))
multiepoch_chi = np.zeros((n_stars, n_epochs, n_filters))
multiepoch_sharp = np.zeros((n_stars, n_epochs, n_filters))
multiepoch_times = np.zeros((n_stars, n_epochs, n_filters))

n_frames = np.zeros((n_stars, n_filters))
mean_mags = np.zeros((n_stars, n_filters))
mean_mag_errs = np.zeros((n_stars, n_filters))
mean_mag_std = np.zeros((n_stars, n_filters))
chi = np.zeros(n_stars)
sharp = np.zeros(n_stars)

multiepoch_mags[:,:,:] = np.nan
multiepoch_errs[:,:,:] = np.nan
multiepoch_chi[:,:,:] = np.nan
multiepoch_sharp[:,:,:] = np.nan
multiepoch_times[:,:,:] = np.nan
mean_mags[:,:] = np.nan
mean_mag_errs[:,:] = np.nan
mean_mag_std[:,:] = np.nan
chi[::] = np.nan
sharp[:] = np.nan

for cc in range(n_filters):
    ch = channels[cc]

    diag = PdfPages('Diagnostics/{}_{}_mag_err.pdf'.format(target, ch))

    for ee in range(n_epochs):
        epoch = ee + 1
        epoch_list = '{}_{}_e{}_mosaic_dn.alf'.format(ch, target, epoch)
        data = crrp.dao.read_alf(epoch_list)

        # Make diagnostic plot
        fig = plt.figure()
        plt.scatter(data['mag'], data['err'], s=1, alpha=0.3, color='k')
        plt.xlabel('{} mag'.format(ch))
        plt.ylabel('{} err'.format(ch))
        diag.savefig(fig, rasterized=True)
        plt.close()

        # identify which stars from the master list were identified in this epoch
        stars = sort_ids[np.searchsorted(ids, data['id'], sorter=sort_ids)]

        # open fits image to get observation time
        img_file = '{}_{}_e{}_mosaic_dn.fits'.format(ch, target, epoch)
        hdul = fits.open(img_file)
        hdr = hdul[0].header
        t = Time(hdr['DATE_OBS'], format='isot', scale='utc')
        mjd = t.mjd

        multiepoch_mags[stars, ee, cc] = data['mag']
        multiepoch_errs[stars, ee, cc] = data['err']
        multiepoch_chi[stars, ee, cc] = data['chi']
        multiepoch_sharp[stars, ee, cc] = data['sharp']
        multiepoch_times[stars, ee, cc] = mjd

        n_frames[stars, cc] += 1

    diag.close()


# Use masked arrays to deal with non-detections
masked_mags = np.ma.MaskedArray(multiepoch_mags, mask=np.isnan(multiepoch_mags))
masked_errs = np.ma.MaskedArray(multiepoch_errs, mask=np.isnan(multiepoch_mags))
masked_chi = np.ma.MaskedArray(multiepoch_chi, mask=np.isnan(multiepoch_mags))
masked_sharp = np.ma.MaskedArray(multiepoch_sharp, mask=np.isnan(multiepoch_mags))

# Compute the weighted mean, standard error, and standard deviation
# Note the standard deviation is not around the weighted mean
mean_mags[:,:] = np.ma.average(masked_mags, axis=1, weights=1./masked_errs**2)
mean_mag_errs[:,:] = 1./np.sqrt(np.ma.sum(1./masked_errs**2, axis=1))
mean_mag_std[:,:] = np.ma.std(masked_mags, axis=1)
mean_mags[mean_mags == 0] = np.nan
mean_mag_errs[mean_mag_errs == 1] = np.nan
mean_mag_std[mean_mag_std == 0] = np.nan
# compute the average chi and sharp value (in all filters) for each star
chi_per_filt = np.ma.average(masked_chi, axis=1)
chi[:] = np.nanmean(chi_per_filt, axis=1)
sharp_per_filt = np.ma.average(masked_sharp, axis=1)
sharp[:] = np.nanmean(sharp_per_filt, axis=1)


# Make diagnostic plot
fig, ax = plt.subplots(3,2, sharex='col', sharey='row', figsize=(8,7))
ax[2,0].set_xlabel('[3.6] mag')
ax[2,1].set_xlabel('[4.5] mag')
ax[0,0].set_ylabel('$\sigma$')
#ax[1,0].set_ylabel('std dev')
ax[1,0].set_ylabel('$\chi$')
ax[2,0].set_ylabel('sharp')
ax[0,0].set_ylim(0,0.1)
#ax[1,0].set_ylim(0,0.5)
ax[1,0].set_ylim(0,5)
ax[2,0].set_ylim(-1,1)
ax[0,0].set_xlim(np.nanmin(mean_mags[:,0]), np.nanmax(mean_mags[:,0]))
ax[0,1].set_xlim(np.nanmin(mean_mags[:,1]), np.nanmax(mean_mags[:,1]))
at.AstroPlots.plot_2D_density(mean_mags[:,0], mean_mag_errs[:,0], plt_axes=ax[0,0],
    xlim=[12,22], ylim=[0,0.1], cmap=plt.cm.bone)
at.AstroPlots.plot_2D_density(mean_mags[:,1], mean_mag_errs[:,1], plt_axes=ax[0,1],
    xlim=[12,22], ylim=[0,0.1], cmap=plt.cm.bone)
at.AstroPlots.plot_2D_density(mean_mags[:,0], chi_per_filt[:,0], plt_axes=ax[1,0],
    xlim=[12,22], ylim=[0,5.0], cmap=plt.cm.bone)
at.AstroPlots.plot_2D_density(mean_mags[:,1], chi_per_filt[:,1], plt_axes=ax[1,1],
    xlim=[12,22], ylim=[0,5.0], cmap=plt.cm.bone)
at.AstroPlots.plot_2D_density(mean_mags[:,0], sharp_per_filt[:,0], plt_axes=ax[2,0],
    xlim=[12,22], ylim=[-1,1], cmap=plt.cm.bone)
at.AstroPlots.plot_2D_density(mean_mags[:,1], sharp_per_filt[:,1], plt_axes=ax[2,1],
    xlim=[12,22], ylim=[-1,1], cmap=plt.cm.bone)
plt.savefig('Diagnostics/{}_mean_catalog.pdf'.format(target), format='pdf',
    rasterized=True)


# Make multiepoch fits file
hdu = fits.BinTableHDU.from_columns(
[fits.Column(name='DAO_ID', format='10A', array=ids),
fits.Column(name='X', format='D', array=input_list['x']),
fits.Column(name='Y', format='D', array=input_list['y']),
fits.Column(name='FILTERS', format='2A', array=np.array(['I1', 'I2'])),
fits.Column(name='WAVES', format='2D', array=np.array([3.6, 4.5])),
fits.Column(name='OBS_TIME', format='{}D'.format(n_epochs*n_filters), \
    dim='({},{})'.format(n_filters, n_epochs), array=multiepoch_times),
fits.Column(name='MAGS', format='{}D'.format(n_epochs*n_filters), \
    dim='({},{})'.format(n_filters, n_epochs), array=multiepoch_mags),
fits.Column(name='ERRS', format='{}D'.format(n_epochs*n_filters), \
    dim='({},{})'.format(n_filters, n_epochs), array=multiepoch_errs),
fits.Column(name='CHI', format='{}D'.format(n_epochs*n_filters), \
    dim='({},{})'.format(n_filters, n_epochs), array=multiepoch_chi),
fits.Column(name='SHARP', format='{}D'.format(n_epochs*n_filters), \
    dim='({},{})'.format(n_filters, n_epochs), array=multiepoch_sharp),
fits.Column(name='N_FRAMES', format='2D', array=n_frames),
])
hdu.writeto('multiepoch.fits', overwrite=True)

# Make mean catalog fits file
hdu = fits.BinTableHDU.from_columns(
[fits.Column(name='DAO_ID', format='10A', array=ids),
fits.Column(name='X', format='D', array=input_list['x']),
fits.Column(name='Y', format='D', array=input_list['y']),
fits.Column(name='FILTERS', format='2A', array=np.array(['I1', 'I2'])),
fits.Column(name='WAVES', format='2D', array=np.array([3.6, 4.5])),
fits.Column(name='MEAN_MAG', format='2D', array=mean_mags),
fits.Column(name='MEAN_MAG_ERR', format='2D', array=mean_mag_errs),
fits.Column(name='MEAN_MAG_STD', format='2D', array=mean_mag_std),
fits.Column(name='N_FRAMES', format='2D', array=n_frames),
fits.Column(name='CHI', format='D', array=chi),
fits.Column(name='SHARP', format='D', array=sharp)
])
hdu.writeto('catalog.fits', overwrite=True)

# Make txt version of mean catalog
dt = np.dtype([('id', int), ('x', float), ('y', float), ('m1', float),
    ('e1', float), ('m2', float), ('e2', float),
    ('n', int), ('n2', int), ('chi', float), ('sharp', float)])

data_save = np.c_[ids, input_list['x'], input_list['y'], mean_mags[:,0], \
    mean_mag_errs[:,0], mean_mags[:,1], mean_mag_errs[:,1], \
    n_frames[:,0], n_frames[:,1], chi, sharp]

#data_save = np.array(zip(ids, input_list['x'], input_list['y'], mean_mags[:,0], \
#    mean_mag_errs[:,0], mean_mag_std[:,0], mean_mags[:,1], mean_mag_errs[:,1], \
#    mean_mag_std[:,1], n_frames[:,0], n_frames[:,1], chi, sharp), dtype=dt)
mean_mags[:,0][np.isnan(mean_mags[:,0])] = 99.999
mean_mags[:,1][np.isnan(mean_mags[:,1])] = 99.999
mean_mag_errs[:,0][np.isnan(mean_mag_errs[:,0])] = 9.999
mean_mag_errs[:,1][np.isnan(mean_mag_errs[:,1])] = 9.999
data_save2 = np.c_[ids, input_list['x'], input_list['y'], mean_mags[:,0], \
    mean_mag_errs[:,0], mean_mags[:,1], mean_mag_errs[:,1], \
    n_frames[:,0], n_frames[:,1], chi, sharp]
np.savetxt('catalog.txt', data_save2, fmt='%9i '+2*'%8.3f '+4*'%6.3f '+'%2i %2i %6.3f %6.3f')
