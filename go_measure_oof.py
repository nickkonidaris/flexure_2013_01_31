# Evernote: https://www.evernote.com/shard/s13/sh/2b7f0492-73a3-476a-b232-7bf536694a16/b4f101d108b77553b197db665f6ca7a0

import numpy as np
import pylab as pl
import pyfits as pf

import os
import SEDM
import SEDM.Bias as Bias
import NPK.Centroid as C

path = "/scr2/npk/sedm/raw/2014jan31/"

positions = [[1018, 1020], [1069-2 ,1034-2], [970-2, 1008], [985-2, 1056], 
        [1032, 1058+12], [1039, 993], [1757-3, 320-0], [209-6, 1798]]

positions = map(np.array, positions)

files = os.listdir(path)

results = {}


pl.ion()
cntr = 0


def iterate(img, pos):
    DD = 8

    residual = 5000
    ntry = 0
    prev_centroid = np.array([5e9, 5e9])
    new_pos = pos[:]
    print pos
    while (residual > .2) and (ntry < 10):
        pos = map(np.round, new_pos)
        #SA = filter(dat[pos[0] - DD:pos[0] + DD , pos[1]-DD : pos[1] +DD])
        SA = dat[pos[0] - DD:pos[0] + DD , pos[1]-DD : pos[1] +DD]
        new_centroid = C.wm(SA)
        new_pos[1] = pos[1] + (new_centroid[0] - DD)
        new_pos[0] = pos[0] + (new_centroid[1] - DD)

        residual = np.sqrt(np.sum((new_centroid-prev_centroid)**2))
        prev_centroid = new_centroid

        print new_pos, new_centroid-DD, ntry, residual

        ntry += 1

    return np.array(new_pos), new_centroid, DD, SA
    
def is_out_of_focus(fname):

    nm,h,m,s = fname.rstrip(".fits").split("_")
    h,m,s = map(float, [h,m,s])
    hms = h + m/60. + s/3600.
    low = 20 + 38./60 + 44/3600.
    high= 23 + 32./60 + 01/3600.

    if low <= hms <= high:
        return True
    else:
        return False

def filter(SA):
    t = SA.copy()
    
    cutoff = 500
    ones = t>cutoff
    zeros = t<=cutoff
    t[ones]=1.0
    t[zeros]=0.0

    return t

#files = ['ifu20140131_20_45_06.fits']
for file in files:
    if is_out_of_focus(file):
        pass
    else:
        continue

    print file
    FF = pf.open(os.path.join(path, file))

    if "OBJECT" not in FF[0].header:
        print "skipping"
        continue

    obj = FF[0].header['OBJECT']
    if "Hg" not in obj:
        print "skipping"
        continue

    dat = Bias.remove(FF)

    cms = []
    for i in range(len(positions)):
        p = positions[i].astype(np.float)

        pos, centroid, DD, SA = iterate(dat, p[:])
        cms.append(pos)

        if len(files) == 1:
            '''This is for debugging purposes'''
            print "fig", p
            pl.figure()
            pl.imshow(SA)
            pl.title("Pos: %s" % p)
            pl.plot(centroid[0], centroid[1], 'o')

            if False:
                pl.figure()
                pl.plot(np.mean(SA, axis=0), 'b', drawstyle='steps-mid')
                pl.plot(np.mean(SA, axis=1), 'r', drawstyle='steps-mid')
                pl.axvline(centroid[0],color='green')
                pl.axvline(centroid[1],color='green')

    if True:
        pl.figure(int(round(cntr / 36)))
        pl.subplot(6,6,cntr % 36)
        cntr += 1
        pl.plot(np.mean(SA, axis=0), 'b', drawstyle='steps-mid')
        pl.plot(np.mean(SA, axis=1), 'r', drawstyle='steps-mid')
        pl.axvline(centroid[0], color='green')
        pl.axvline(centroid[1], color='green')
        pl.ylim([0,3000])

        pl.text(-5,1000, file)




    results[file] = {'offsets': cms}
    results[file]['Object'] = obj
    HA = FF[0].header['HA']
    if HA[0] == 'W': sign = -1.0
    else: sign = 1.0

    h,m,s = map(float,HA[1:].split(':'))
    
    deg = sign * h * 15.0 + m * 15.0/60.0 + s*15.0/3600.0
    results[file]['HA'] = deg

    o = FF[0].header['Object']
    try:
        the_ha = float(o.split("[")[1].split(",")[0])*15.0*-1
        print deg, the_ha
    except:
        the_ha = deg

    deg = the_ha

    Dec = FF[0].header['DEC']
    h,m,s = map(float,Dec.split(':'))
    deg = h + m/60.0 + s/3600.0
    results[file]['Dec'] = deg
    
    results[file]['Max-Min'] = np.max(SA) - np.min(SA)
    results[file]['exptime'] = FF[0].header['EXPTIME']
    results[file]['Bgd'] = np.median(dat)

pl.show()

import pickle
pickle.dump( results, open( "results_2014_jan_31_oof.p", "wb" ) )


