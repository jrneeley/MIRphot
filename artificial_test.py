import numpy as np
import sys
import os
import shutil
#sys.path.insert(0, '/home/jill/python/')
import daophot_tools as dao
#import CRRPphot as crrp
from astropy import wcs

target = sys.argv[1]
xoff = float(sys.argv[2])
yoff = float(sys.argv[3])
n_art_stars = int(sys.argv[4])
n_tests = int(sys.argv[5])

channels = ['I1', 'I2']
exptime = 23.6
frametime = 30
fluxconv = np.array([0.1257, 0.1447])
pixratio = 4.0
phpadu1 = exptime*3.7/(fluxconv[0]*pixratio)
phpadu2 = exptime*3.7/(fluxconv[1]*pixratio)
n_filters = 2
n_epochs = 12

#n_tests = 20
#n_art_stars = 500

os.chdir(target)

##### SETUP DIRECTORY STRUCTURE
print('Setting up directory structure...')

# Make new directory for Artificial Star tests
directory = 'ArtificialStar'
if not os.path.exists(directory): os.makedirs(directory)
os.chdir(directory)

# Copy images into this directory
shutil.copy('../median.fits', 'median.fits')
shutil.copy('../median.psf', 'median.psf')
for i in range(12):
    img = 'I1_{}_e{}_mosaic_dn.fits'.format(target, i+1)
    psf = 'I1_{}_e{}_mosaic_dn.psf'.format(target, i+1)
    shutil.copy('../'+img, 'I1_e{}_orig.fits'.format(i+1))
    shutil.copy('../'+psf, 'I1_e{}_orig.psf'.format(i+1))
    img = 'I2_{}_e{}_mosaic_dn.fits'.format(target, i+1)
    psf = 'I2_{}_e{}_mosaic_dn.psf'.format(target, i+1)
    shutil.copy('../'+img, 'I2_e{}_orig.fits'.format(i+1))
    shutil.copy('../'+psf, 'I2_e{}_orig.psf'.format(i+1))

daophot_setup.get_irac_opt_files(channels, frametime, warm=1, mosaic=1)

##### GENERATE IMAGES WITH ARTIFICIAL STARS

print('Making images with added stars...')

# Generate lists of artifical stars
med_img = 'median.fits'

dao.dao.addstar(med_img, file_stem='median_fake_', num_images=n_tests,
    seed=5, gain=phpadu1, min_mag=12, max_mag=22, num_stars=n_art_stars,
    opt_file='I1-daophot.opt', verbose=0)

for t in range(n_tests):

    # Convert xy coordinates to wcs coord using image header.
    # read in artificial star list
    dt = np.dtype([('x', float), ('y', float)])
    median_pix = np.loadtxt('median_fake_{:02}.add'.format(t+1), dtype=dt,
        skiprows=3, usecols=(1,2))

    # get ra/dec from image
    w = wcs.WCS(med_img)
    ra, dec = w.all_pix2world(median_pix['x']+xoff, median_pix['y']+yoff, 1)


    for i in range(12):

        starlist = 'I1_e{}_fake_{:02}.add'.format(i+1, t+1)
        shutil.copy('median_fake_{:02}.add'.format(t+1), starlist)
        img = 'I1_e{}_orig.fits'.format(i+1)
        w2 = wcs.WCS(img)
        x, y = w2.all_world2pix(ra, dec, 1)

        with open(starlist, 'r') as f:
            all_lines = f.readlines()
        with open(starlist, 'w') as f:
            for j, line in enumerate(all_lines,1):
                if j <= 3:
                    f.writelines(line)
                else:
                    # split line into array
                    temp = line.split()
                    # combine back into string
                    x_new, y_new = float(x[j-4]), float(y[j-4])
                    new_line = '{:6} {:8.3f} {:8.3f} {:8}\n'.format(temp[0], x_new, y_new, temp[3])
                    # save to file
                    f.writelines(new_line)

        dao.dao.addstar(img, file_stem='temp',
            star_list=starlist, seed=5, gain=phpadu1,
            opt_file='I1-daophot.opt', verbose=0)
        shutil.move('temp01.fits', 'I1_e{}_fake_{:02}.fits'.format(i+1, t+1))

        starlist = 'I2_e{}_fake_{:02}.add'.format(i+1, t+1)
        shutil.copy('median_fake_{:02}.add'.format(t+1), starlist)
        img = 'I2_e{}_orig.fits'.format(i+1)
        w2 = wcs.WCS(img)
        x, y = w2.all_world2pix(ra, dec, 1)

        with open(starlist, 'r') as f:
            all_lines = f.readlines()
        with open(starlist, 'w') as f:
            for j, line in enumerate(all_lines,1):
                if j <= 3:
                    f.writelines(line)
                else:
                    # split line into array
                    temp = line.split()
                    # combine back into string
                    x_new, y_new = float(x[j-4]), float(y[j-4])
                    new_line = '{:6} {:8.3f} {:8.3f} {:8}\n'.format(temp[0], x_new, y_new, temp[3])
                    # save to file
                    f.writelines(new_line)

        dao.dao.addstar(img, file_stem='temp',
            star_list=starlist, seed=5, gain=phpadu2,
            opt_file='I2-daophot.opt', verbose=0)
        shutil.move('temp01.fits', 'I2_e{}_fake_{:02}.fits'.format(i+1, t+1))

    ##### RUN PHOTOMETRY PIPELINE
    print('Doing photometry on new median image...')


    img_list = []

    for ch in channels:
        for i in range(12):

            dn_stem = '{}_e{}_fake_{:02}'.format(ch, i+1, t+1)
            if os.path.exists(dn_stem+'.als') == False:
                img_list.append(dn_stem+'.ap')
                # Do first round of initial photometry
                opt_file = '{}-daophot.opt'.format(ch)
                # Run find and phot on each single epoch mosaic
                dao.dao.find(dn_stem, num_frames='5,1', opt_file=opt_file,
                    new_thresh=10.0, verbose=0)
                dao.dao.phot(dn_stem, opt_file=opt_file, verbose=0)
                # copy master PSF to every other frame
                master_psf = '{}_e1_orig.psf'.format(ch)
                shutil.copy(master_psf, dn_stem+'.psf')
                dao.dao.allstar(dn_stem, suppress=1)

    # Skip steps of creating transformations and making the median image.
    # We can just use the files from the initial run.
    shutil.copy('../all-als.mch', 'all-als.mch')
    # Change names of files!
    with open('all-als.mch', 'r') as f:
        all_lines = f.readlines()
    with open('all-als.mch', 'w') as f:
        for j, line in enumerate(all_lines,1):
            # split line into array
            temp = line.split()
            if j <= 12:
                file_new = 'I1_e{}_fake_{:02}.als'.format(j, t+1)
            else:
                file_new = 'I2_e{}_fake_{:02}.als'.format(j-12, t+1)
            # combine back into string
            sep = ' '
            #new_line = '\''+file_new+sep.join(temp[1:])+'\n'
            new_line = '\'{:25} {}\n'.format(file_new, sep.join(temp[1:]))
            # save to file
            f.writelines(new_line)


    # Do photometry on median image
    num_iters = 1

    median_stem = 'median_fake_{:02}'.format(t+1)
    if os.path.exists(median_stem+'.als') == False:
        # Do first round of initial photometry
        opt_file = 'I1-daophot.opt'
        # For the median image roughly 60 frames have been averaged, (5 dithers * 12 aors).
        # Accounting for this adjusts the readout noise and gain from the opt file
        dao.dao.find(median_stem, num_frames='60,1', new_thresh=2, opt_file=opt_file, verbose=1)
        dao.dao.phot(median_stem, opt_file=opt_file, verbose=0)

        # copy master PSF to every other frame
        master_psf = 'median.psf'
        shutil.copy(master_psf, median_stem+'.psf')
        # Run allstar
        dao.dao.allstar(median_stem, verbose=1)

        # Repeat on process on subtracted images
        for jj in range(num_iters):
            print('\n\nIteration {} of {}'.format(jj+2, num_iters+1))
            # make copy of subtracted image for archiving
            sub_img = median_stem+'s.fits'
            shutil.copy(sub_img, '{}_sub_{}.fits'.format(median_stem, jj+1))
            # Find additional stars in subtracted image
            dao.dao.find(sub_img, num_frames='60,1', new_thresh=3, opt_file=opt_file, verbose=0)
            dao.dao.offset(median_stem+'s.coo', id_offset=(jj+1)*100000)
            # Append new stars to list
            cmb_file = '{}_{}.cmb'.format(median_stem, jj+1)
            if jj == 0:
                dao.dao.append(median_stem+'.coo', median_stem+'s.off', \
                    out_file=cmb_file, verbose=1)
            else:
                dao.dao.append('{}_{}.cmb'.format(median_stem, jj), median_stem+'s.off', \
                    out_file=cmb_file, verbose=1)
            # Run daophot/allstar again with new list
            dao.dao.phot(median_stem, opt_file=opt_file, coo_file=cmb_file,
                ap_file=median_stem+'.ap', verbose=0)
            dao.dao.allstar(median_stem, verbose=1)

    starlist = dao.read_dao.read_alf(median_stem+'.als')
    # Change ID in the artificial star lists to new DAO id
    with open(median_stem+'.add', 'r') as f:
        all_lines = f.readlines()
    with open(median_stem+'.add', 'w') as f:
        for j, line in enumerate(all_lines,1):
            if j <= 3:
                f.writelines(line)
            else:
                temp = line.split()
                x_art, y_art = float(temp[1]), float(temp[2])
                d = np.sqrt((starlist['x']-x_art)**2 + (starlist['y']-y_art)**2)
                dao_id = starlist['id'][d == np.min(d)][0]
                # combine back into string
                new_line = '{:7} {:8} {:8} {:8}\n'.format(dao_id, temp[1], temp[2], temp[3])
                # save to file
                f.writelines(new_line)

    dao.other.dao2reg(median_stem+'.add', median_stem, ids=1, radius=1,
        color='red')

    # offset the final .mag file with the xy offset from montage2
    dao.dao.offset(median_stem+'.als', x_offset=xoff, y_offset=yoff)
    dao.dao.allframe('all-als.mch', median_stem+'.off', verbose=1)
