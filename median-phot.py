import sys
import daophot_tools as dao
import os
import glob
import numpy as np
import re
import shutil
import matplotlib.pyplot as plt

target = sys.argv[1]
montage_x = sys.argv[2]
montage_y = sys.argv[3]

os.chdir(target)
num_iters = 1


median_stem = 'median'

# Do first round of initial photometry
print('Iteration 1 of {}'.format(num_iters+1))
print('    Running daophot find and phot ...')
opt_file = 'I1-daophot.opt'
# For the median image roughly 60 frames have been averaged, (5 dithers * 12 aors).
# Accounting for this adjusts the readout noise and gain from the opt file
dao.dao.find(median_stem, num_frames='5,1', new_thresh=2, opt_file=opt_file, verbose=1)
dao.dao.phot(median_stem, opt_file=opt_file, verbose=0)

# copy master PSF to every other frame
master_psf = '{}I1_0p6_pixscale_mosaic_dn.psf'.format(dao.config.psf_dir)
shutil.copy(master_psf, median_stem+'.psf')
# Run allstar
print('    Running allstar...')
dao.dao.allstar(median_stem, verbose=1)
dao.other.dao2reg('median.als', 'median-1', ids=0, radius=3)

# Repeat on process on subtracted images
for jj in range(num_iters):
    print('\n\nIteration {} of {}'.format(jj+2, num_iters+1))
    # make copy of subtracted image for archiving
    sub_img = median_stem+'s.fits'
    shutil.copy(sub_img, 'median_sub_{}.fits'.format(jj+1))
    # Find additional stars in subtracted image
    print('    Running daophot find...')
    dao.dao.find(sub_img, num_frames='5,1', new_thresh=3, opt_file=opt_file, verbose=0)
    print('    Running offset and append...')
    dao.dao.offset(median_stem+'s.coo', id_offset=(jj+1)*100000)
    dao.other.dao2reg(median_stem+'s.coo', 'median-{}'.format(jj+2), ids=0)
    # Append new stars to list
    cmb_file = '{}_{}.cmb'.format(median_stem, jj+1)
    if jj == 0:
        dao.dao.append(median_stem+'.coo', median_stem+'s.off', \
            out_file=cmb_file, verbose=1)
    else:
        dao.dao.append('{}_{}.cmb'.format(median_stem, jj), median_stem+'s.off', \
            out_file=cmb_file, verbose=1)
    # Run daophot/allstar again with new list
    print('    Running daophot phot...')
    dao.dao.phot(median_stem, opt_file=opt_file, coo_file=cmb_file,
        ap_file=median_stem+'.ap', verbose=0)
    print('    Running allstar...')
    dao.dao.allstar(median_stem, verbose=1)



# offset the final .mag file with the xy offset from montage2
dao.dao.offset('median.als', x_offset=montage_x, y_offset=montage_y)
dao.other.dao2reg('median.als', 'median-final', ids=0, radius=3)
dao.other.dao2reg('median.off', 'frame1', ids=0, radius=3)
