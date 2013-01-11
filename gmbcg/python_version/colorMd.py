"""
This script calcualte the color model based on the DES mock. 
"""
import pyfits as pf
import numpy as np
import pylab as pl
import binplot as bp


data0 = pf.getdata('/home/jghao/research/data/des_mock/v3.04/truthCat/DES_Mock_v3.04_Baseline_truth_00.fit')
data1 = pf.getdata('/home/jghao/research/data/des_mock/v3.04/truthCat/DES_Mock_v3.04_Baseline_truth_01.fit')
data2 = pf.getdata('/home/jghao/research/data/des_mock/v3.04/truthCat/DES_Mock_v3.04_Baseline_truth_02.fit')
data3 = pf.getdata('/home/jghao/research/data/des_mock/v3.04/truthCat/DES_Mock_v3.04_Baseline_truth_03.fit')

bg0 = data0[data0.field('central') == 1]
bg1 = data1[data1.field('central') == 1]
bg2 = data2[data2.field('central') == 1]
bg3 = data3[data3.field('central') == 1]

gmr = np.concatenate((bg0.field('omag')[:,0]-bg0.field('omag')[:,1],bg1.field('omag')[:,0]-bg1.field('omag')[:,1],bg2.field('omag')[:,0]-bg2.field('omag')[:,1],bg3.field('omag')[:,0]-bg3.field('omag')[:,1]))

rmi = np.concatenate((bg0.field('omag')[:,1]-bg0.field('omag')[:,2],bg1.field('omag')[:,1]-bg1.field('omag')[:,2],bg2.field('omag')[:,1]-bg2.field('omag')[:,2],bg3.field('omag')[:,1]-bg3.field('omag')[:,2]))

imz = np.concatenate((bg0.field('omag')[:,2]-bg0.field('omag')[:,3],bg1.field('omag')[:,2]-bg1.field('omag')[:,3],bg2.field('omag')[:,2]-bg2.field('omag')[:,3],bg3.field('omag')[:,2]-bg3.field('omag')[:,3]))

zmy = np.concatenate((bg0.field('omag')[:,3]-bg0.field('omag')[:,4],bg1.field('omag')[:,3]-bg1.field('omag')[:,4],bg2.field('omag')[:,3]-bg2.field('omag')[:,4],bg3.field('omag')[:,3]-bg3.field('omag')[:,4]))

z = np.concatenate((bg0.field('z'), bg1.field('z'),bg2.field('z'),bg3.field('z')))

idxgr = z < 0.4
idxri = (z >=0.4) * (z <0.75)
idxiz = (z >= 0.75) * (z< 1.15)
idxzy = z >= 1.15

pl.figure(figsize=(14,10))
pl.subplot(2,2,1)
resgmr = np.polyfit(z[idxgr],gmr[idxgr],1)
pl.plot(z[idxgr],gmr[idxgr],'b.',alpha=0.5)
bp.bin_scatter(z[idxgr],gmr[idxgr],binsize = 0.05,fmt = 'ro',scatter=True)
pl.plot(z[idxgr],np.polyval(resgmr,z[idxgr]),'r,')
pl.xlabel('z')
pl.ylabel('g - r')
pl.title('slope, intercept:'+str(resgmr))

pl.subplot(2,2,2)
resrmi = np.polyfit(z[idxri],rmi[idxri],1)
pl.plot(z[idxri],rmi[idxri],'b.',alpha=0.5)
bp.bin_scatter(z[idxri],rmi[idxri],binsize = 0.05,fmt = 'ro',scatter=True)
pl.plot(z[idxri],np.polyval(resrmi,z[idxri]),'r,')
pl.xlabel('z')
pl.ylabel('r - i')
pl.title('slope, intercept:'+str(resrmi))

pl.subplot(2,2,3)
resimz = np.polyfit(z[idxiz],imz[idxiz],1)
pl.plot(z[idxiz],imz[idxiz],'b.',alpha=0.5)
bp.bin_scatter(z[idxiz],imz[idxiz],binsize = 0.05,fmt = 'ro',scatter=True)
pl.plot(z[idxiz],np.polyval(resimz,z[idxiz]),'r,')
pl.xlabel('z')
pl.ylabel('i - z')
pl.title('slope, intercept:'+str(resimz))

pl.subplot(2,2,4)
reszmy = np.polyfit(z[idxzy],zmy[idxzy],1)
pl.plot(z[idxzy],zmy[idxzy],'b.',alpha=0.5)
bp.bin_scatter(z[idxzy],zmy[idxzy],binsize = 0.05,fmt = 'ro',scatter=True)
pl.plot(z[idxzy],np.polyval(reszmy,z[idxzy]),'r,')
pl.xlabel('z')
pl.ylabel('z - y')
pl.title('slope, intercept:'+str(reszmy))

pl.savefig('/home/jghao/research/data/des_mock/v3.04/truthCat/des_mockV3.04_color.png')
#-----the inversion -------
# y = a x + b -> x = 1/a * y - b/a

def gmrz(gmr):
    z = 1./2.681*gmr - 0.665/2.681
    return z

def rmiz(rmi):
    z = 1./2.063*rmi - (-0.266/2.063)
    return z

def imzz(imz):
    z = 1./1.667*imz - (-0.708/1.667)
    return z

def zmyz(zmy):
    z = 1./0.192 * zmy - (0.223/0.192)
    return z


# verify ---
pl.plot(gmr[idxgr],z[idxgr],'b,')
pl.plot(gmr[idxgr],gmrz(gmr[idxgr]),'r,')

pl.plot(rmi[idxri],z[idxri],'b,')
pl.plot(rmi[idxri],rmiz(rmi[idxri]),'r,')

pl.plot(imz[idxiz],z[idxiz],'b,')
pl.plot(imz[idxiz],imzz(imz[idxiz]),'r,')

pl.plot(zmy[idxzy],z[idxzy],'b,')
pl.plot(zmy[idxzy],zmyz(zmy[idxzy]),'r,')
