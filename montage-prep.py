import sys
#sys.path.insert(0,'/home/jill/python/')
#import CRRPphot as crrp
from . import daophot_setup
import daophot_tools as dao
import os
import numpy as np
import shutil

targets = ['NGC7078', 'NGC3201', 'NGC6402', 'NGC6121', 'NGC5904']

s = 2
target = targets[s]


channels = ['I1', 'I2']
exptime = 23.6
frametime = 30
fluxconv = [0.1257, 0.1447]
pixratio = 4
starting_dir = os.getcwd()
n_epochs = 12

# Navigate to target directory
os.chdir(target)

img_list = []
for ii, ch in enumerate(channels):

    # copy appropriate daophot options files to current directory
    daophot_setup.get_irac_opt_files(channels, frametime, warm=1, mosaic=1)

    print 'Prepping {} images....'.format(ch)
    for jj in range(n_epochs):
        flux_image = '{}_{}_e{}_mosaic.fits'.format(ch, target, jj+1)
        daophot_setup.spitzer_flux2dn(flux_image, exptime=exptime, \
            fluxconv=fluxconv[ii], pixratio=pixratio)

        print '    Working on single epoch mosaic {} of {}'.format(jj+1, n_epochs)

        dn_stem = '{}_{}_e{}_mosaic_dn'.format(ch, target, jj+1)
        img_list.append(dn_stem+'.ap')
        # Do first round of initial photometry
        #print '    Running daophot...'
        opt_file = '{}-daophot.opt'.format(ch)
        # Run find and phot on each single epoch mosaic
        dao.dao.find(dn_stem, num_frames='5,1', opt_file=opt_file, new_thresh=10.0, verbose=1)
        dao.dao.phot(dn_stem, opt_file=opt_file, verbose=0)
        # copy master PSF to every other frame
        master_psf = '{}{}_mosaic_012919.psf'.format(dao.config.psf_dir, ch)
        shutil.copy(master_psf, dn_stem+'.psf')
        dao.dao.allstar(dn_stem, suppress=1)

# match the aperture photometry files
dao.dao.daomatch(img_list, 'all-ap.mch', force_scale_rot=0)
dao.dao.daomaster('all-ap.mch', frame_num='10,0.5,12', sigma='1')
# A good transformation is more important than a complete catalog at this stage
