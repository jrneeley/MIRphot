#!/usr/bin/env python
import re
import shutil
from astropy.io import fits
from daophot_tools import config


def spitzer_flux2dn(image, newname="", exptime=0, fluxconv=0, pixratio=1):
	# exptime should be total (median) integration time
	# fluxconv is the value for a single BCD 

	if (newname == ""):
		newname = re.sub(".fits", "_dn.fits", image)
	shutil.copy(image, newname)
	hdul = fits.open(newname, mode='update')
	prihdr = hdul[0].header
	scidata = hdul[0].data

	# if no exptime or fluxconv is given, get it from the header. 
	# WARNING: this will only be accurate for individual BCDs! 	
	if exptime == 0 : exptime = prihdr['exptime']
	if fluxconv == 0 : fluxconv = prihdr['fluxconv']
	
	# change default fluxconv if not in native pixel ratio
	if pixratio != 1: fluxconv *= pixratio

	scidata *= exptime/fluxconv
	hdul.close()

def get_irac_opt_files(filters, exptime, warm=1, mosaic=1):

        opt_dir = config.opt_dir+'IRAC/'
        if warm == 1:
                opt_dir2 = opt_dir+'warm/'
        if warm == 0:
                opt_dir2 = opt_dir+'cryo/'

        shutil.copy(opt_dir+'allframe.opt', 'allframe.opt')
        shutil.copy(opt_dir+'dummy-daophot.opt', 'daophot.opt')
        if mosaic == 1:
                shutil.copy(opt_dir+'photo-mosaic.opt', 'photo.opt')
                shutil.copy(opt_dir+'allstar-mosaic.opt', 'allstar.opt')
                for filt in filters:
                        opt_file = opt_dir2+filt+'-0p6-pixscale-mosaic.opt'
                        shutil.copy(opt_file, filt+'-daophot.opt')
        if mosaic == 0:
                shutil.copy(opt_dir+'photo.opt', 'photo.opt')
                shutil.copy(opt_dir+'allstar.opt', 'allstar.opt')
                for filt in filters:
                        opt_file = opt_dir2+filt+'-bcd.opt'
                        shutil.copy(opt_file, filt+'-daophot.opt')


def combine_mch_simple(mch_list, output_file='combine.mch'):

        o = open(output_file, 'w')
        for ii, mch in enumerate(mch_list):

                f = open(mch, 'r')
                lines = f.readlines()
                if ii == 0:
                        o.write(lines[0])
                for jj in range(1, len(lines)):
                        o.write(lines[jj])

                f.close()
        o.close()


def copy_mch_artificial(als_list, mch_file, output_file='copy.mch'):

        o = open(output_file, 'w')

        f = open(mch_file, 'r')
        mch_lines = f.readlines()

        o.write(mch_lines[0])

        for ii, als in enumerate(als_list):

                newname = als.ljust(40)
                oldline = mch_lines[1]
                #print newline[1:40]
                #newline[1:40] = newname
                newline = ' \''+newname+oldline[42:]
                o.write(newline)
        f.close()
        o.close()
