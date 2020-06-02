import numpy as np
import sys
import os
import shutil
sys.path.insert(0, '/home/jill/python/')
import daophot_tools as dao
from astropy.io import fits

target = sys.argv[1]

channels = ['I1', 'I2']
n_filters = 2
n_epochs = 12
n_tests = 40

os.chdir(target+'/ArtificialStar')

f = open('artificial_stars.txt', 'w')

for i in range(n_tests):

    dt = np.dtype([('id', int), ('x', float), ('y', float), ('mag_true', float)])
    input_list = np.loadtxt('median_fake_{:02}.add'.format(i+1),
        dtype=dt, skiprows=3)

    ids = input_list['id']
    n_stars = len(input_list['id'])

    multiepoch_mags = np.zeros((n_stars, n_epochs, n_filters))
    multiepoch_errs = np.zeros((n_stars, n_epochs, n_filters))

    n_found = np.zeros((n_stars, n_filters))
    n_expected = np.zeros((n_stars, n_filters))
    mean_mags = np.zeros((n_stars, n_filters))
    mean_mag_errs = np.zeros((n_stars, n_filters))

    multiepoch_mags[:,:,:] = np.nan
    multiepoch_errs[:,:,:] = np.nan
    mean_mags[:,:] = np.nan
    mean_mag_errs[:,:] = np.nan

    for cc in range(n_filters):
        ch = channels[cc]

        for ee in range(n_epochs):
            epoch = ee + 1
            epoch_stem = '{}_e{}_fake_{:02}'.format(ch, epoch, i+1)

            # check to see if artificial star is in fov of this image
            dt = np.dtype([('x', float), ('y', float)])
            img_xy = np.loadtxt(epoch_stem+'.add', usecols=(1,2), skiprows=3,
                dtype=dt)
            img = fits.open(epoch_stem+'.fits')
            img_data = img[0].data
            x = img_xy['x'].astype(int)
            y = img_xy['y'].astype(int)
            n = len(x)
            for s in range(n):
                # Check if stars are within the image limits
                if (y[s] >= 0) & (y[s] < img_data.shape[0]) & (x[s] >= 0) & \
                    (x[s] < img_data.shape[1]):
                    # Check if star is too close to the edge
                    #if np.any(data[y[i]-10:y[i]+10, x[i]-10:x[i]+10] < -500):
                    #    use[ee] = 0
                    #else:
                    if img_data[y[s], x[s]] > -500:
                        n_expected[s,cc] += 1


            data = dao.read_dao.read_alf(epoch_stem+'.alf')

            for s in range(n_stars):
                star = data['id'] == ids[s]
                if len(data['id'][star]) == 1:
                    multiepoch_mags[s, ee, cc] = data['mag'][star]
                    multiepoch_errs[s, ee, cc] = data['err'][star]
                    n_found[s, cc] += 1



    # Use masked arrays to deal with non-detections
    masked_mags = np.ma.MaskedArray(multiepoch_mags, mask=np.isnan(multiepoch_mags))
    masked_errs = np.ma.MaskedArray(multiepoch_errs, mask=np.isnan(multiepoch_mags))

    # Compute the weighted mean, standard error, and standard deviation
    # Note the standard deviation is not around the weighted mean
    mean_mags[:,:] = np.ma.average(masked_mags, axis=1, weights=1./masked_errs**2)
    mean_mag_errs[:,:] = 1./np.sqrt(np.ma.sum(1./masked_errs**2, axis=1))
    mean_mags[mean_mags == 0] = np.nan
    mean_mag_errs[mean_mag_errs == 1] = np.nan

    mean_mags[:,0][np.isnan(mean_mags[:,0])] = 99.999
    mean_mags[:,1][np.isnan(mean_mags[:,1])] = 99.999
    mean_mag_errs[:,0][np.isnan(mean_mag_errs[:,0])] = 9.999
    mean_mag_errs[:,1][np.isnan(mean_mag_errs[:,1])] = 9.999

    data_save2 = np.c_[ids, input_list['x'], input_list['y'], input_list['mag_true'],
        mean_mags[:,0], mean_mag_errs[:,0], mean_mags[:,1], mean_mag_errs[:,1], \
        n_found[:,0], n_found[:,1], n_expected[:,0], n_expected[:,1]]
    np.savetxt(f, data_save2, fmt='%9i '+2*'%8.3f '+5*'%6.3f '+4*'%2i ')
