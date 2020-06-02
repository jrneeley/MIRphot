import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, '/home/jill/python/')
import daophot_tools as dao
import AstroTools as at
from astropy import wcs
from matplotlib.backends.backend_pdf import PdfPages

target = sys.argv[1]
center_ra = sys.argv[2]
center_dec = sys.argv[3]
xoff = float(sys.argv[4])
yoff = float(sys.argv[5])

os.chdir(target+'/ArtificialStar/')


dt = np.dtype([('id', int), ('x', float), ('y', float), ('mag_true', float),
    ('mag1', float), ('e1', float), ('mag2', float), ('e2', float), ('n1', int),
    ('n2', int), ('ne1', int), ('ne2', int)])
data = np.loadtxt('artificial_stars.txt', dtype=dt)

data['mag1'][data['mag1'] > 90] = np.nan
data['mag2'][data['mag2'] > 90] = np.nan

g = data['ne1'] == 12
g2 = data['ne2'] == 12

fig, ax = plt.subplots(2,1, figsize=(10,8))
bins = np.linspace(10,22,num=25)

ax[0].hist(data['mag_true'][g], color='k', histtype='step', bins=bins,
    linewidth=1.5)
ax[0].hist(data['mag1'][g][~np.isnan(data['mag1'][g])], color='xkcd:steel blue',
    histtype='step', bins=bins, linewidth=1.5)
ax[1].hist(data['mag_true'][g2], color='k', histtype='step', bins=bins,
    linewidth=1.5)
ax[1].hist(data['mag2'][g2][~np.isnan(data['mag2'][g2])], color='xkcd:red',
    histtype='step', bins=bins, linewidth=1.5)
ax[0].set_xlabel('[3.6] mag')
ax[0].set_ylabel('N')
ax[1].set_xlabel('[4.5] mag')
ax[1].set_ylabel('N')
ax[0].text(0.05, 0.9, 'N expected = {}'.format(len(data['mag_true'][g])),
    transform=ax[0].transAxes)
ax[0].text(0.05, 0.85, 'N found = {}'.format(len(data['mag1'][g][~np.isnan(data['mag1'][g])])),
    transform=ax[0].transAxes)
ax[1].text(0.05, 0.9, 'N expected = {}'.format(len(data['mag_true'][g2])),
    transform=ax[1].transAxes)
ax[1].text(0.05, 0.85, 'N found = {}'.format(len(data['mag2'][g2][~np.isnan(data['mag2'][g2])])),
    transform=ax[1].transAxes)
plt.savefig('hist.pdf')

# compute radial distance from center of cluster
c_ra, c_dec = at.coordinates.radec_string2deg(center_ra, center_dec)
w = wcs.WCS('median.fits')
x, y = w.all_world2pix(c_ra, c_dec, 1)
x -= xoff
y -= yoff
r = np.sqrt( (data['x']-x)**2 + (data['y']-y)**2)

fig, ax = plt.subplots(2,1, figsize=(10,8))
cb1 = ax[0].scatter(data['x'], data['y'], c=data['n1'], s=1)
cb2 = ax[1].scatter(data['x'], data['y'], c=data['n2'], s=1)
cbar1 = plt.colorbar(cb1, ax=ax[0])
cbar2 = plt.colorbar(cb2, ax=ax[1])
cbar1.set_label('N Frames')
cbar2.set_label('N Frames')
ax[1].set_xlabel('X')
ax[1].set_ylabel('Y')
ax[0].set_ylabel('Y')
plt.savefig('map.pdf')



diff1 = data['mag1'] - data['mag_true']
diff2 = data['mag2'] - data['mag_true']

fig, ax = plt.subplots(2,1, figsize=(10,8))
c1 = ax[0].scatter(data['mag_true'][g], diff1[g], c=r[g], s=1)
c2 = ax[1].scatter(data['mag_true'][g2], diff2[g2], c=r[g2], s=1)
cb1 = plt.colorbar(c1, ax=ax[0])
cb2 = plt.colorbar(c2, ax=ax[1])
cb1.set_label('R')
cb2.set_label('R')
ax[1].set_xlabel('true mag')
ax[0].set_ylabel('[3.6] - true mag')
ax[1].set_ylabel('[4.5] - true mag')
plt.savefig('mag_diff.pdf')


fig, ax = plt.subplots(2,1, figsize=(10,8))
cb1 = ax[0].scatter(r[g], diff1[g], c=data['mag_true'][g], s=1)
cb2 = ax[1].scatter(r[g2], diff2[g2], c=data['mag_true'][g2], s=1)
cbar1 = plt.colorbar(cb1, ax=ax[0])
cbar2 = plt.colorbar(cb2, ax=ax[1])
cbar1.set_label('True mag')
cbar2.set_label('True mag')
ax[1].set_xlabel('R (pixels)')
ax[0].set_ylabel('[3.6] - true mag')
ax[1].set_ylabel('[4.5] - true mag')
plt.savefig('r_diff.pdf')


s = diff1[g] <= -1.0
s2 = diff2[g2] <= -1.0
fig, ax = plt.subplots(2,1, figsize=(10,8))
cb1 = ax[0].scatter(data['x'][g][s], data['y'][g][s], c=diff1[g][s], s=10)
cb2 = ax[1].scatter(data['x'][g2][s2], data['y'][g2][s2], c=diff2[g2][s2], s=10)
cbar1 = plt.colorbar(cb1, ax=ax[0])
cbar2 = plt.colorbar(cb2, ax=ax[1])
cbar1.set_label('$\Delta$ mag')
cbar2.set_label('$\Delta$ mag')
ax[1].set_xlabel('X')
ax[1].set_ylabel('Y')
ax[0].set_ylabel('Y')
plt.savefig('map-diff.pdf')

fig, ax = plt.subplots(2,2, figsize=(12,8))
s = np.abs(diff1[g]) <= 0.05
s2 = np.abs(diff2[g2]) <= 0.05
n1_e, _, _ = ax[0,0].hist(data['mag_true'][g], color='k', histtype='step', bins=bins,
    linewidth=1.5)
n1_f, _, _ = ax[0,0].hist(data['mag1'][g][s], color='xkcd:steel blue',
    histtype='step', bins=bins, linewidth=1.5)
n2_e, _, _ = ax[1,0].hist(data['mag_true'][g2], color='k', histtype='step', bins=bins,
    linewidth=1.5)
n2_f, _, _ = ax[1,0].hist(data['mag2'][g2][s2], color='xkcd:red',
    histtype='step', bins=bins, linewidth=1.5)
ax[0,0].set_xlabel('[3.6] mag')
ax[0,0].set_ylabel('N')
ax[1,0].set_xlabel('[4.5] mag')
ax[1,0].set_ylabel('N')
ax[0,1].text(0.7, 0.9, 'N expected = {}'.format(len(data['mag_true'][g])),
    transform=ax[0,1].transAxes)
ax[0,1].text(0.7, 0.85, 'N found = {}'.format(len(data['mag1'][g][s])),
    transform=ax[0,1].transAxes)
ax[1,1].text(0.7, 0.9, 'N expected = {}'.format(len(data['mag_true'][g2])),
    transform=ax[1,1].transAxes)
ax[1,1].text(0.7, 0.85, 'N found = {}'.format(len(data['mag2'][g2][s2])),
    transform=ax[1,1].transAxes)
incomp1 = n1_f.astype(float)/n1_e.astype(float)
incomp2 = n2_f.astype(float)/n2_e.astype(float)

ax[0,1].plot(bins[:-1], incomp1, color='xkcd:steel blue')
ax[1,1].plot(bins[:-1], incomp2, color='xkcd:red')
ax[0,1].set_ylim(0,1)
ax[1,1].set_ylim(0,1)
ax[0,1].set_ylabel('Completeness Fraction')
ax[0,1].set_xlabel('[3.6] mag')
ax[1,1].set_ylabel('Completeness Fraction')
ax[1,1].set_xlabel('[4.5] mag')

plt.savefig('completeness-mag.pdf')


bins = np.linspace(0,np.max(r),num=25)
fig, ax = plt.subplots(2,2, figsize=(12,8))
s = np.abs(diff1[g]) <= 0.05
s2 = np.abs(diff2[g2]) <= 0.05
n1_e, _, _ = ax[0,0].hist(r[g], color='k', histtype='step', bins=bins,
    linewidth=1.5)
n1_f, _, _ = ax[0,0].hist(r[g][s], color='xkcd:steel blue',
    histtype='step', bins=bins, linewidth=1.5)
n2_e, _, _ = ax[1,0].hist(r[g2], color='k', histtype='step', bins=bins,
    linewidth=1.5)
n2_f, _, _ = ax[1,0].hist(r[g2][s2], color='xkcd:red',
    histtype='step', bins=bins, linewidth=1.5)
ax[0,0].set_xlabel('R (pixels)')
ax[0,0].set_ylabel('N')
ax[1,0].set_xlabel('R (pixels)')
ax[1,0].set_ylabel('N')
ax[0,1].text(0.7, 0.9, 'N expected = {}'.format(len(r[g])),
    transform=ax[0,1].transAxes)
ax[0,1].text(0.7, 0.85, 'N found = {}'.format(len(r[g][s])),
    transform=ax[0,1].transAxes)
ax[1,1].text(0.7, 0.9, 'N expected = {}'.format(len(r[g2])),
    transform=ax[1,1].transAxes)
ax[1,1].text(0.7, 0.85, 'N found = {}'.format(len(r[g2][s2])),
    transform=ax[1,1].transAxes)
incomp1 = n1_f.astype(float)/n1_e.astype(float)
incomp2 = n2_f.astype(float)/n2_e.astype(float)

ax[0,1].scatter(bins[:-1], incomp1, color='xkcd:steel blue')
ax[1,1].scatter(bins[:-1], incomp2, color='xkcd:red')
ax[0,1].set_ylim(0,1)
ax[1,1].set_ylim(0,1)
ax[0,1].set_ylabel('Completeness Fraction')
ax[0,1].set_xlabel('R (pixels)')
ax[1,1].set_ylabel('Completeness Fraction')
ax[1,1].set_xlabel('R (pixels)')

plt.savefig('completeness-radial.pdf')

plt.show()

# radial completeness in different magnitude bins
bins_mag = np.linspace(12,22,num=51)
bins_rad = np.linspace(0,np.max(r),num=25)
pdf = PdfPages('test.pdf')

plt.rc('font', family='sans-serif', size='8')
plt.rc('ytick', labelsize='x-small')
plt.rc('xtick', labelsize='x-small')

for i in range(len(bins_mag)-1):
    m = (data['mag_true'][g] > bins_mag[i]) & (data['mag_true'][g] <= bins_mag[i+1])
    m2 = (data['mag_true'][g2] > bins_mag[i]) & (data['mag_true'][g2] <= bins_mag[i+1])
    delta_mag1 = diff1[g][m]
    delta_mag2 = diff2[g2][m2]
    rad1 = r[g][m]
    rad2 = r[g2][m2]

    m = (data['mag_true'][g][s] > bins_mag[i]) & (data['mag_true'][g][s] <= bins_mag[i+1])
    m2 = (data['mag_true'][g2][s2] > bins_mag[i]) & (data['mag_true'][g2][s2] <= bins_mag[i+1])
    delta_mag1_clean = diff1[g][s][m]
    delta_mag2_clean = diff2[g2][s2][m2]
    rad1_clean = r[g][s][m]
    rad2_clean = r[g2][s2][m2]

    fig, ax = plt.subplots(2,2, figsize=(12,8))
    for j in range(len(bins_rad)-1):

        r1 = (rad1 < bins_rad[j]) & (rad1 <= bins_rad[j+1])
        r2 = (rad2 < bins_rad[j]) & (rad2 <= bins_rad[j+1])
        mean_bin = (bins_rad[j]+bins_rad[j+1])/2.0
        med1 = np.nanmedian(delta_mag1[r1])
        med2 = np.nanmedian(delta_mag2[r2])
        ax[0,0].scatter(mean_bin, med1, marker='s', color='xkcd:black', alpha=0.7)
        ax[1,0].scatter(mean_bin, med2, marker='s', color='xkcd:black', alpha=0.7)
        ax[0,0].errorbar(mean_bin, np.nanmean(delta_mag1[r1]),
            yerr=np.nanstd(delta_mag1[r1]), fmt='o', color='xkcd:steel blue',
            alpha=0.7)
        ax[1,0].errorbar(mean_bin, np.nanmean(delta_mag2[r2]),
            yerr=np.nanstd(delta_mag2[r2]), fmt='o', color='xkcd:red',
            alpha=0.7)

        r1 = (rad1_clean < bins_rad[j]) & (rad1_clean <= bins_rad[j+1])
        r2 = (rad2_clean < bins_rad[j]) & (rad2_clean <= bins_rad[j+1])
        mean_bin = (bins_rad[j]+bins_rad[j+1])/2.0
        med1 = np.nanmedian(delta_mag1_clean[r1])
        med2 = np.nanmedian(delta_mag2_clean[r2])
        ax[0,1].scatter(mean_bin, med1, marker='s', color='xkcd:black', alpha=0.7)
        ax[1,1].scatter(mean_bin, med2, marker='s', color='xkcd:black', alpha=0.7)
        ax[0,1].errorbar(mean_bin, np.nanmean(delta_mag1_clean[r1]),
            yerr=np.nanstd(delta_mag1_clean[r1]), fmt='o', color='xkcd:steel blue',
            alpha=0.7)
        ax[1,1].errorbar(mean_bin, np.nanmean(delta_mag2_clean[r2]),
            yerr=np.nanstd(delta_mag2_clean[r2]), fmt='o', color='xkcd:red',
            alpha=0.7)



    ax[0,0].set_xlim(0,np.max(r))
    ax[1,0].set_xlim(0,np.max(r))
    ax[0,0].axhline(0, color='k')
    ax[1,0].axhline(0, color='k')
    ax[1,0].set_xlabel('R (pixels)')
    ax[0,0].set_ylabel('[3.6] - true mag')
    ax[1,0].set_ylabel('[4.5] - true mag')
    ax[0,1].set_xlim(0,np.max(r))
    ax[1,1].set_xlim(0,np.max(r))
    ax[0,1].axhline(0, color='k')
    ax[1,1].axhline(0, color='k')
    ax[1,1].set_xlabel('R (pixels)')
    ax[0,1].set_ylabel('[3.6] - true mag')
    ax[1,1].set_ylabel('[4.5] - true mag')
    ax[0,0].set_title('{:.2f} < mag < {:.2f}'.format(bins_mag[i], bins_mag[i+1]))
    ax[0,0].set_ylim(np.nanmin(diff1[g]), np.nanmax(diff1[g]))
    ax[1,0].set_ylim(np.nanmin(diff2[g2]), np.nanmax(diff2[g2]))
    ax[0,1].set_ylim(-0.05, 0.05)
    ax[1,1].set_ylim(-0.05, 0.05)
    pdf.savefig(fig)
    plt.close(fig)
#plt.savefig('test.pdf')

pdf.close()
#plt.show()
