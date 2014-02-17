# Evernote: https://www.evernote.com/shard/s13/sh/2b7f0492-73a3-476a-b232-7bf536694a16/b4f101d108b77553b197db665f6ca7a0

import numpy as np
import pylab as pl
import pyfits as pf

import os
import SEDM
import SEDM.Bias as Bias
import NPK.Centroid as C

path = "/scr2/npk/sedm/raw/2014jan31/"

positions = [[1018, 1020], [1069-2 ,1034], [970-2, 1008], [985-2, 1056], 
        [1032, 1058+12], [1039, 993], [1757, 320-2], [209, 1798]]

files = os.listdir(path)

results = {}


pl.ion()
cntr = 0

def is_out_of_focus(fname):

    nm,h,m,s = fname.rstrip(".fits").split("_")
    h,m,s = map(float, [h,m,s])
    hms = h + m/60. + s/3600.
    low = 20 + 38./60 + 44/3600.
    high= 23 + 38./60 + 01/3600.

    if low <= hms <= high:
        return True
    else:
        return False

#files = ['ifu20140131_20_29_41.fits']
#files = ['ifu20140131_18_15_33.fits']
#files = ['ifu20140131_18_30_06.fits']
for file in files:
    print file
    if is_out_of_focus(file):
        print "out of focus"
        continue

    FF = pf.open(os.path.join(path, file))

    if "OBJECT" not in FF[0].header:
        print "skipping"
        continue

    obj = FF[0].header['OBJECT']
    if "Hg" not in obj:
        print "skipping"
        continue

    dat = Bias.remove(FF)

    DD = 8
    cms = []
    for i in range(len(positions)):

        p = positions[i]
        SA = dat[p[0]-DD:p[0]+DD , p[1]-DD:p[1]+DD]
        SA -= np.median(SA)

        centroid = C.wm(SA)

        cms.append([p[1] - centroid[1] + DD, p[0] - centroid[0] + DD])
        if len(files) == 1:
            '''This is for debugging purposes'''
            print "fig", p
            pl.figure()
            pl.imshow(SA)
            pl.plot(centroid[0], centroid[1], 'o')

            if False:
                pl.figure()
                pl.plot(np.mean(SA, axis=0), 'b', drawstyle='steps-mid')
                pl.plot(np.mean(SA, axis=1), 'r', drawstyle='steps-mid')
                pl.axvline(centroid[0],color='blue')
                pl.axvline(centroid[1],color='red')



    if False:
        pl.figure(int(round(cntr / 36)))
        pl.subplot(6,6,cntr % 36)
        cntr += 1
        pl.plot(np.mean(SA, axis=0), 'b', drawstyle='steps-mid')
        pl.plot(np.mean(SA, axis=1), 'r', drawstyle='steps-mid')
        pl.axvline(centroid[0], color='blue')
        pl.axvline(centroid[1], color='red')
        pl.ylim([0,2000])




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
pickle.dump( results, open( "results_2014_jan_31.p", "wb" ) )


