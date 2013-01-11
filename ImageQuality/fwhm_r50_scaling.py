import pylab as pl, numpy as np, binplot as bp, glob as gl, cPickle as p
from DECamCCD import linefit

def plot45line(x,y):
    a = np.array([min(min(x),min(y)),max(max(x),max(y))])
    pl.plot(a,a,'r-')
    return 0



b=np.genfromtxt('r50_whk_whkrms_phi_fwhm_nonwindowed.txt')
ok = b[:,1] != -999.
b = b[ok,:]
r50 = b[:,0]
fwhm = b[:,4]

pl.figure(figsize=(18,6))
pl.subplot(1,3,1)
pl.plot(fwhm,r50,'b.')
pl.xlabel('FWHM [arcsec]')
pl.ylabel('R50 [arcsec]')
bp.bin_scatter(fwhm,r50,nbins=15,fmt='ro')
res = linefit(fwhm,r50)
pl.plot(np.array([0.6,2.5]),res[0]+res[1]*np.array([0.6,2.5]),'r-')
pl.grid()
pl.xlim(0.6,2.5)
pl.ylim(0.4,1.4)
pl.title('R50 ='+str(round(res[1],4))+'*FWHM + ' + str(round(res[0],4)))

low = fwhm <= 0.9
high = fwhm >= 1

pl.subplot(1,3,2)
pl.plot(fwhm[low],r50[low],'b.')
pl.xlabel('FWHM [arcsec]')
pl.ylabel('R50 [arcsec]')
bp.bin_scatter(fwhm[low],r50[low],nbins=8,fmt='ro')
res = linefit(fwhm[low],r50[low])
pl.plot(np.array([0.6,0.92]),res[0]+res[1]*np.array([0.6,0.92]),'r-')
pl.grid()
pl.xlim(0.72,0.92)
pl.ylim(0.45,0.55)
pl.title('R50 ='+str(round(res[1],4))+'*FWHM + ' + str(round(res[0],4)))


pl.subplot(1,3,3)
pl.plot(fwhm[high],r50[high],'b.')
pl.xlabel('FWHM [arcsec]')
pl.ylabel('R50 [arcsec]')
bp.bin_scatter(fwhm[high],r50[high],nbins=13,fmt='ro')
res = linefit(fwhm[high],r50[high])
pl.plot(np.array([0.98,2.5]),res[0]+res[1]*np.array([0.98,2.5]),'r-')
pl.grid()
pl.xlim(0.98,2.5)
pl.ylim(0.5,1.4)
pl.title('R50 ='+str(round(res[1],4))+'*FWHM + ' + str(round(res[0],4)))

pl.figtext(0.4,0.95,'all based on Firstcut catalogs, choose different ranges',color='b')
pl.savefig('fwhm_world_r50_scaling_firstcut.png')
pl.close()


#---from my run of sextractor ----

f = gl.glob('/home/jghao/research/desSV/11242012/allband/desIQ*.p') + gl.glob('/home/jghao/research/desSV/12-21-2012-sv-analysis/desIQ*.p')

nfile = len(f)
fwhm = np.zeros(nfile)
r50 = np.zeros(nfile)
for i in range(nfile):
    r50[i]=p.load(open(f[i],'r'))[8]
    fwhm[i] = p.load(open(f[i],'r'))[9]

pl.figure(figsize=(18,6))
pl.subplot(1,3,1)
pl.plot(fwhm,r50,'b.')
pl.xlabel('FWHM [arcsec]')
pl.ylabel('R50 [arcsec]')
bp.bin_scatter(fwhm,r50,nbins=15,fmt='ro')
res = linefit(fwhm,r50)
pl.plot(np.array([0.6,2.5]),res[0]+res[1]*np.array([0.6,2.5]),'r-')
pl.grid()
pl.xlim(0.6,2.5)
pl.ylim(0.4,1.4)
pl.title('R50 ='+str(round(res[1],4))+'*FWHM + ' + str(round(res[0],4)))

low = fwhm <= 0.9
high = fwhm >= 1

pl.subplot(1,3,2)
pl.plot(fwhm[low],r50[low],'b.')
pl.xlabel('FWHM [arcsec]')
pl.ylabel('R50 [arcsec]')
bp.bin_scatter(fwhm[low],r50[low],nbins=8,fmt='ro')
res = linefit(fwhm[low],r50[low])
pl.plot(np.array([0.6,0.92]),res[0]+res[1]*np.array([0.6,0.92]),'r-')
pl.grid()
pl.xlim(0.72,0.92)
pl.ylim(0.45,0.55)
pl.title('R50 ='+str(round(res[1],4))+'*FWHM + ' + str(round(res[0],4)))


pl.subplot(1,3,3)
pl.plot(fwhm[high],r50[high],'b.')
pl.xlabel('FWHM [arcsec]')
pl.ylabel('R50 [arcsec]')
bp.bin_scatter(fwhm[high],r50[high],nbins=13,fmt='ro')
res = linefit(fwhm[high],r50[high])
pl.plot(np.array([0.98,2.5]),res[0]+res[1]*np.array([0.98,2.5]),'r-')
pl.grid()
pl.xlim(0.98,2.5)
pl.ylim(0.5,1.4)
pl.title('R50 ='+str(round(res[1],4))+'*FWHM + ' + str(round(res[0],4)))

pl.figtext(0.4,0.95,'all based on my run of sextractor, choose different ranges',color='b')

pl.savefig('fwhm_world_r50_scaling_myrun_sextractor.png')
pl.close()


#whkrms ----

f = gl.glob('/home/jghao/research/desSV/11242012/allband/desIQ*.p') + gl.glob('/home/jghao/research/desSV/12-21-2012-sv-analysis/desIQ*.p')

nfile = len(f)
whk=np.zeros(nfile)
whkSex=np.zeros(nfile)
whkrms = np.zeros(nfile)
whkrmsSex = np.zeros(nfile)
r50 = np.zeros(nfile)
r50Sex = np.zeros(nfile)
for i in range(nfile):
    whk[i]=p.load(open(f[i],'r'))[1]
    whkSex[i]=p.load(open(f[i],'r'))[5]
    whkrms[i]=p.load(open(f[i],'r'))[3]
    whkrmsSex[i] = p.load(open(f[i],'r'))[7]
    r50Sex[i] = p.load(open(f[i],'r'))[8]
    r50[i] = p.load(open(f[i],'r'))[4]


#compare my run with firstcut ----

f = gl.glob('/home/jghao/research/desSV/firstcut_compare/desIQ*.p')
f.sort()
nfile = len(f)
whk=np.zeros(nfile)
whkSex=np.zeros(nfile)
whkrms = np.zeros(nfile)
whkrmsSex = np.zeros(nfile)
r50 = np.zeros(nfile)
r50Sex = np.zeros(nfile)
fwhmSex = np.zeros(nfile)
expid = np.zeros(nfile).astype('i')
for i in range(nfile):
    expid[i] = p.load(open(f[i],'r'))[0]
    whk[i]=p.load(open(f[i],'r'))[1]
    whkrms[i]=p.load(open(f[i],'r'))[3]
    r50[i] = p.load(open(f[i],'r'))[4]
    whkSex[i]=p.load(open(f[i],'r'))[5]
    whkrmsSex[i] = p.load(open(f[i],'r'))[7]
    r50Sex[i] = p.load(open(f[i],'r'))[8] 
    fwhmSex[i] = p.load(open(f[i],'r'))[9]

fwhm = r50*2.
b=np.genfromtxt('/home/jghao/research/desSV/firstcut_compare/firstcut_image_quality_final_forcompare.txt',delimiter=',')

fwhk = b[:,2]
fwhkrms = b[:,3]
fr50 = b[:,1]
ffwhm = b[:,5]

pl.figure(figsize=(15,15))
pl.subplot(3,3,1)
pl.plot(whk,whkSex,'b.')
plot45line(whk,whkSex)
pl.grid()
pl.xlabel('whk  [arcsec]')
pl.ylabel('whk_sextractor [arcsec]')
pl.subplot(3,3,2)
pl.plot(whk,fwhk,'b.')
plot45line(whk,fwhk)
pl.grid()
pl.xlabel('whk  [arcsec]')
pl.ylabel('whk_firstcut [arcsec]')
pl.subplot(3,3,3)
pl.plot(whkSex,fwhk,'b.')
plot45line(whkSex,fwhk)
pl.grid()
pl.xlabel('whk_sextractor  [arcsec]')
pl.ylabel('whk_firstcut [arcsec]')

pl.subplot(3,3,4)
pl.plot(whkrms,whkrmsSex,'b.')
plot45line(whkrms,whkrmsSex)
pl.grid()
pl.xlabel('whkrms  [arcsec]')
pl.ylabel('whkrms_sextractor [arcsec]')
pl.subplot(3,3,5)
pl.plot(whkrms,fwhkrms,'b.')
plot45line(whkrms,fwhkrms)
pl.grid()
pl.xlabel('whkrms  [arcsec]')
pl.ylabel('whkrms_firstcut [arcsec]')
pl.subplot(3,3,6)
pl.plot(whkrmsSex,fwhkrms,'b.')
plot45line(whkrmsSex,fwhkrms)
pl.grid()
pl.xlabel('whkrms_sextractor  [arcsec]')
pl.ylabel('whkrms_firstcut [arcsec]')

pl.subplot(3,3,7)
pl.plot(r50,r50Sex,'b.')
plot45line(r50,r50Sex)
pl.grid()
pl.xlabel('r50  [arcsec]')
pl.ylabel('r50_sextractor [arcsec]')

pl.subplot(3,3,8)
pl.plot(r50,fr50,'b.')
plot45line(r50,fr50)
pl.grid()
pl.xlabel('r50  [arcsec]')
pl.ylabel('r50_firstcut [arcsec]')
pl.subplot(3,3,9)
pl.plot(r50Sex,fr50,'b.')
plot45line(r50Sex,fr50)
pl.grid()
pl.xlabel('r50_sextractor  [arcsec]')
pl.ylabel('r50_firstcut [arcsec]')

pl.figtext(0.35,0.95,'comparison between my measurement and firstcut',color='b',fontsize=19)
pl.savefig('firstcut_myrun_compare.png')
pl.close()


# ---all firstcut distribution -----
b=np.genfromtxt('/home/jghao/research/ggsvn/des-sv/firstcut_image_quality_12282012.txt',delimiter=',')
ok = b[:,2] != -999.
b = b[ok,]
fwhk = b[:,2]
fwhkrms = b[:,3]
fr50 = b[:,1]
ffwhm = b[:,5]

pl.figure(figsize=(18,6))
pl.subplot(1,3,1)
pl.hist(fr50,bins=20)
pl.xlabel('r50_firstcut [arcsec]')

pl.subplot(1,3,2)
pl.hist(fwhk,bins=20)
pl.xlabel('whk_firstcut [arcsec]')

pl.subplot(1,3,3)
pl.hist(fwhkrms,bins=20)
pl.xlabel('whkrms_firstcut [arcsec]')

pl.figtext(0.4,0.95,'based on firstcut catalogs as of 12/28/2012',color='b',fontsize=19)

pl.savefig('firstcut_cat_statistics.png')
pl.close()

# ---weighted measurement distribution -----
f = gl.glob('/home/jghao/research/desSV/11242012/allband/desIQ*.p') + gl.glob('/home/jghao/research/desSV/12-21-2012-sv-analysis/desIQ*.p')+gl.glob('/home/jghao/research/desSV/firstcut_compare/desIQ*.p')

nfile = len(f)
whk=np.zeros(nfile)
whkrms = np.zeros(nfile)
r50 = np.zeros(nfile)
for i in range(nfile):
    whk[i]=p.load(open(f[i],'r'))[1]
    whkrms[i]=p.load(open(f[i],'r'))[3]
    r50[i] = p.load(open(f[i],'r'))[4]

pl.figure(figsize=(18,6))
pl.subplot(1,3,1)
pl.hist(r50,bins=20)
pl.xlabel('r50_weighted [arcsec]')

pl.subplot(1,3,2)
pl.hist(whk,bins=20)
pl.xlabel('whk_weighted [arcsec]')

pl.subplot(1,3,3)
pl.hist(whkrms,bins=20)
pl.xlabel('whkrms_weighted [arcsec]')

pl.figtext(0.2,0.95,'based on my weighted moments from 147 exposures from 11/24,12/09,12/16, 12/18',color='b',fontsize=19)

pl.savefig('weighted_cat_statistics.png')
pl.close()


