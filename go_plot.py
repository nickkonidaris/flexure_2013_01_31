import numpy as np
import pylab as pl

import pandas
import pickle


results = pickle.load(open("results_2014_jan_31.p", "rb"))


fidname = 'ifu20140131_17_46_54'
fids = results[fidname + '.fits']['offsets']


dxss = [ [], [], [], [] ]
dyss = [ [], [], [], [] ]
HAs = []
Decs = []

repeat_x = [ [], [], [], [] ]
repeat_y = [ [], [], [], [] ]

zs = [ [], [], [], [] ]
ze = [ [], [], [], [] ]

Bgd = []
MM = []
et = []
times = []
center_dx = []
center_dy = []
center_sx = []
center_sy = []

run_starts = []

sorted_files = results.keys()
sorted_files.sort()

counter = 0
for file in sorted_files:
    res = results[file]
    #print file
    r = res['offsets']


    h,m,s =  map(float,file.rstrip(".fits").split("_")[1:4])
    if h > 23: continue

    if "Zenith End" in res['Object']:
        for i in range(len(dxss)):
            ze[i].append(r[i][0] - fids[i][0])
            ze[i].append(r[i][1] - fids[i][1])

        res['HA'] = 0

        run_starts.append(counter+1)
    if "Zenith Start" in res['Object']:
        print file
        for i in range(len(dxss)):
            zs[i].append(r[i][0] - fids[i][0])
            zs[i].append(r[i][1] - fids[i][1])
    
        res['HA'] = 0
        run_starts.append(counter)

    counter += 1
        
    if "Repeat" in res['Object']:
        for i in range(len(dxss)):
            repeat_x[i].append(r[i][0] - fids[i][0])
            repeat_y[i].append(r[i][1] - fids[i][1])
    else:
        HAs.append(res['HA'])
        Decs.append(res['Dec'])
        Bgd.append(res['Bgd'])
        MM.append(res['Max-Min'])
        et.append(res['exptime'])

        dxs = [r[i][0] - fids[i][0] for i in range(4)]
        dys = [r[i][1] - fids[i][1] for i in range(4)]

        center_dx.append(np.median(dxs))
        center_dy.append(np.median(dys))
        center_sx.append(np.std(dxs))
        center_sy.append(np.std(dys))
        
        if np.std(dxs) > .2: print "** ", file
        if h == 0: h =24
        times.append(h + m/60. + s/3600.)
        for i in range(len(dxss)):
            dxss[i].append(r[i][0] - fids[i][0])
            dyss[i].append(r[i][1] - fids[i][1])
    

dxss = np.array(dxss)
dyss = np.array(dyss)
center_dx, center_dy, center_sx, center_sy = map(np.array,
        [center_dx, center_dy, center_sx, center_sy])
HAs, Decs = map(np.array, [HAs, Decs])

pl.ion()
pl.figure(1)
pl.clf()
pl.figure(2)
pl.clf()
for i in range(len(dxss)):
    pl.figure(1)
    pl.plot(dxss[i], dyss[i], 'o')
    pl.figure(2)
    pl.plot(dxss[i], dyss[i]-dyss[0], 'o')




df = pandas.DataFrame({'HA': HAs, 'Decs': Decs, 'Dx0': dxss[0], 'Dy0': dyss[0],
    'Dx1': dxss[1], 'Dy1': dyss[1]})

#pandas.tools.plotting.scatter_matrix(df, alpha=0.3)
fig = pl.figure(3)
ax = fig.add_subplot(111)
ax.set_aspect(1)
pl.clf()
pl.errorbar(center_dx, center_dy, center_sx, center_sy, 'o')
pl.xlim(-4,4)
pl.ylim(-4,4)
pl.xlabel("Pixel Dx")
pl.ylabel("Pixel Dy")
pl.title("Flexure [ZP @ %s]" % fidname)
pl.grid(which='both')

fig = pl.figure(4)
ax = fig.add_subplot(111)
ax.set_aspect(1)
pl.clf()
pl.scatter(center_dx, center_dy, c=HAs)
pl.xlim(-4,4)
pl.ylim(-4,4)
pl.colorbar()
pl.xlabel("Pixel Dx")
pl.ylabel("Pixel Dy")
pl.title("Flexure | HA [ZP @ %s" % fidname)
pl.grid(which='both')

fig = pl.figure(5)
ax = fig.add_subplot(111)
ax.set_aspect('equal')
pl.clf()
pl.scatter(center_dx, center_dy, c=Decs)
pl.xlim(-4,4)
pl.ylim(-4,4)
pl.colorbar()
pl.xlabel("Pixel Dx")
pl.ylabel("Pixel Dy")
pl.title("Flexure | Dec [ZP @ %s" % fidname)
pl.grid(which='both')


starts = []
for i in xrange(len(run_starts)-1):
    if run_starts[i]+1 == run_starts[i+1]:
        pass
    else:
        starts.append(run_starts[i-1])


tracks = []

counter = 10

for ix in xrange(len(starts)-1):
    start = starts[ix]
    end = starts[ix+1]


    files = sorted_files[start:end]

    if len(files)<3:
        continue

    print "**"
    t = {"HA": [],
            "Dec": [],
            "dxs": [],
            "dys": [],
            "cdxs": [],
            "cdys": []}

    fids = results[files[0]]['offsets']

    for file in files:
        print results[file]['Object']
        res = results[file]


        dxs = [res['offsets'][i][0] - fids[i][0] for i in range(4)]
        dys = [res['offsets'][i][1] - fids[i][1] for i in range(4)]

        t["HA"].append(res['HA'])
        t["Dec"].append(res['Dec'])
        t["cdxs"].append(np.median(dxs))
        t["cdys"].append(np.median(dys))

    for k,v in t.items():
        t[k] = np.array(v)

    pl.figure(counter)
    counter += 1
    pl.clf()
    pl.plot(t["HA"], t["cdxs"], 'ob-')
    pl.plot(t["HA"], t["cdys"], 'or-')
    pl.title("Run at Dec %3.1f" % np.mean(t['Dec']))
    pl.xlabel("HA deg")
    pl.ylabel("Offset")
    pl.grid(True)
    pl.legend(["Delta X", "Delta Y"])
    pl.xlim(-100,100)
    pl.ylim(-5,5)
    tracks.append(t)


    pl.figure(99)
    pl.plot(t["cdxs"], t["cdys"], 'o-')



#### FITTING
import core

XS = np.array([center_dy, HAs, Decs]).T
ff = core.multipolyfit(XS, center_dx, 2, model_out=True)
res=map(ff, XS[:,0], XS[:,1], XS[:,2])
ok = np.abs(center_dx - res) < .3

X1 = np.array([center_dy[ok], HAs[ok], Decs[ok]]).T
ff = core.multipolyfit(X1, center_dx[ok], 2, model_out=True)
res=map(ff, XS[:,0], XS[:,1], XS[:,2])
print "Residual: %1.2f" % np.std(center_dx-res)
ff = core.multipolyfit(X1, center_dx[ok], 2, powers_out=True)
print ff

fig = pl.figure(200)
ax = fig.add_subplot(111)
ax.set_aspect(1)
pl.clf()
pl.scatter(center_dx-res, center_dy, c=Decs)
pl.xlim(-4,4)
pl.ylim(-4,4)
pl.colorbar()
pl.xlabel("Residual Dx")
pl.ylabel("Dy")
pl.title("Residual | Dec [ZP @ %s" % fidname)
pl.grid(which='both')


def nff(dy, ha, dec):
    ''' Dy in pixels. HA in degrees, Dec in degrees '''
    
    res = (-9.76170735e-02 
        -5.40857690e-01 * dy
        -1.42770036e-02 * ha
        -8.24304496e-03 * dec
        +9.89255017e-03 * dy**2
        +3.53396833e-04 * dy * ha
        +2.31293635e-03 * dy * dec
        +7.37365997e-05 * ha**2
        +8.83937326e-05 * ha * dec
        -6.34635265e-06 * dec**2
        )
    return res


pl.figure(500)
pl.clf()
pl.scatter(-HAs, Decs,c=(center_dx), s=60)
pl.xlabel("HA [Positive west]")
pl.ylabel("Dec")
cb=pl.colorbar()
cb.set_label("Pixel shift")
pl.grid()
pl.savefig("ha_dec_offset.pdf")



pl.figure(501)
pl.clf()
pl.scatter(-HAs, Decs,c=(center_dx-res), s=60)
pl.xlabel("HA [Positive west]")
pl.ylabel("Dec")
cb=pl.colorbar()
cb.set_label("Pixel shift")
pl.grid()
pl.savefig("ha_dec_residuals.pdf")
