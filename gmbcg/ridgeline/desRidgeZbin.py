#! /usr/bin/env python
# This code calculate the ridgeline colors using ecgmm for the des mock----
# J. Hao @ Fermilab, 4/1/2012
# This is copied from mockRidgeZbin.py and modified to apply to the real des data. J. Hao, 12/2/2012 

import numpy as np
import pyfits as pf
import pylab as pl
import esutil as es
import ecgmmPy as gmm
import rwecgmmPy as rwgmm
import scipy.stats as sts
import glob as gl
import time
import binplot as bp


def zgmr(z):
    gmr = 2.681*z + 0.665
    return gmr

def zrmi(z):
    rmi = 2.063*z -0.266
    return rmi

def zimz(z):
    imz = 1.667*z-0.708
    return imz

def zzmy(z):
    zmy = 0.192 * z + 0.223
    return zmy

#-----0.4L* in z-band ---------
def limz(x):
    """
    corresponding i_absmag <= -20.5 
    """
    A=np.exp(3.17)
    k=0.15
    lmz=A*x**k
    return(lmz)


def ecgmmRidgeZbin(z=None,color=None,colorErr=None,mag=None):
    #--define some quantity to be returned ----
    ok = (color>=-1) * (color <= 5.)
    color = color[ok]
    colorErr = colorErr[ok]
    mag = mag[ok]
    z = z[ok]
    alpha0 = []
    alpha1 = []
    mu0 = []
    mu1 = []
    sigma0 = []
    sigma1 = []
    ntot = []
    amp = []
    aic1 = []
    aic2 = []
    zmd=[]
    #-----------------------------------------
    startTime=time.time()
    h,rev = es.stat.histogram(z,binsize=0.01,rev=True)
    for i in range(h.size):
        if rev[i] != rev[i+1]:
            print i
            indices = rev[rev[i]:rev[i+1]]
            zmed = np.mean(z[indices])
            ok = mag[indices] <= limz(zmed)
            indices = indices[ok]
            if len(indices) <=10:
                continue
            alpha=np.array([0.5,0.5])
            mu=np.array([sts.scoreatpercentile(color[indices],per=70),sts.scoreatpercentile(color[indices],per=40)])
            sigma=np.array([0.04,0.3])
            a2 = gmm.aic_ecgmm(color[indices],colorErr[indices],alpha,mu,sigma)           
            a1 = gmm.wstat(color[indices],colorErr[indices])[2]
            if a2 < a1:
                srt=np.argsort(sigma)
                alpha0.append(alpha[srt[0]])
                alpha1.append(alpha[srt[1]])
                mu0.append(mu[srt[0]])
                mu1.append(mu[srt[1]])
                sigma0.append(sigma[srt[0]])
                sigma1.append(sigma[srt[1]])
                aic1.append(a1)
                aic2.append(a2)
                amp.append(len(indices)*alpha[srt[0]])
                zmd.append(zmed)
                print a2, a1
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)
    return np.array(alpha0),np.array(alpha1),np.array(mu0),np.array(mu1),np.array(sigma0),np.array(sigma1),np.array(aic1),np.array(aic2),np.array(amp),np.array(zmd)





#------catalogs -------------------

if __name__ == '__main__':
    from desRidgeZbin import *
    catdir = '/home/jghao/research/data/DES_photoz_challenge/'
    #gal = pf.getdata(catdir+'vvds_des_mag_detmodel.fits',1)
    gal = pf.getdata(catdir+'vvds_des_mag_auto.fits',1)
    ra = gal.field('ra')
    dec = gal.field('dec')
    z = gal.field('spz')
    mag=gal.field('mag_z')-10
    gmr=gal.field('mag_g') - gal.field('mag_r')
    gmrErr=np.sqrt(gal.field('magerr_g')**2+gal.field('magerr_r')**2)
    rmi=gal.field('mag_r') - gal.field('mag_i')
    rmiErr=np.sqrt(gal.field('magerr_r')**2+gal.field('magerr_i')**2)
    imz=gal.field('mag_i') - gal.field('mag_z')
    imzErr=np.sqrt(gal.field('magerr_i')**2+gal.field('magerr_z')**2)

    #----g - r color -----
    pl.figure(figsize=(10,10))
    pl.subplot(2,2,1)
    idx_gr = z < 0.4
    resGR = ecgmmRidgeZbin(z=z[idx_gr],color=gmr[idx_gr],colorErr=gmrErr[idx_gr],mag=mag[idx_gr])
    color = resGR[2]
    redshift = resGR[9]
    sgm0 = resGR[4]
    sgm1 = resGR[6]
    zz = np.arange(0.05,0.5,0.1)
    bp.bin_scatter(redshift,color,binsize=0.05,fmt='bo',label='measurement')
    pl.plot(zz,zgmr(zz),'r-',label='color model')
    pl.title('des_vvds_mag_detmodel')
    pl.legend(loc='best')
    pl.xlabel('redshift')
    pl.ylabel('Ridgeline Color: g-r')

    # --- r - i color -------
    pl.subplot(2,2,2)
    idx_ri = (z >= 0.4)*(z < 0.75 )
    resRI = ecgmmRidgeZbin(z=z[idx_ri],color=rmi[idx_ri],colorErr=rmiErr[idx_ri],mag=mag[idx_ri])
    color = resRI[2]
    redshift = resRI[9]
    sgm0 = resRI[4]
    sgm1 = resRI[6]
    zz = np.arange(0.4,0.9,0.1)
    bp.bin_scatter(redshift,color,binsize=0.05,fmt='bo',label='measurement')
    pl.plot(zz,zrmi(zz),'r-',label='color model')
    pl.title('des_vvds_mag_detmodel')
    pl.legend(loc='best')
    pl.xlabel('redshift')
    pl.ylabel('Ridgeline Color: r-i')

    #-----i - z color -----------
    pl.subplot(2,2,3)
    idx_iz = (z >= 0.75)*(z < 1.15)
    resIZ = ecgmmRidgeZbin(z=z[idx_iz],color=imz[idx_iz],colorErr=imzErr[idx_iz],mag=mag[idx_iz])
    color = resIZ[2]
    redshift = resIZ[9]
    sgm0 = resIZ[4]
    sgm1 = resIZ[6]
    zz = np.arange(0.75,1.3,0.1)
    bp.bin_scatter(redshift,color,binsize=0.05,fmt='bo',label='measurement')
    pl.plot(zz,zimz(zz),'r-',label='color model')
    pl.title('des_vvds_mag_detmodel')
    pl.legend(loc='best')
    pl.xlabel('redshift')
    pl.ylabel('Ridgeline Color: i-z')
    
    pl.figtext(0.5,0.95,'mag_auto')
    #pl.savefig(catdir+'des_vvds_ridgeline_mag_detmodel_zbin.png')
    pl.savefig(catdir+'des_vvds_ridgeline_mag_auto_zbin.png')
    
    pl.close()

    #------gmr ridgeline distribution ------

    redshift = resGR[9]
    alp0 = resGR[0]
    alp1 = resGR[1]
    mu0 = resGR[2]
    mu1 = resGR[3]
    sgm0 = resGR[4]
    sgm1 = resGR[5]
    z,alp0m,yerr = bp.bin_scatter(redshift,alp0,binsize=0.05,plot=False)
    z,alp1m,yerr = bp.bin_scatter(redshift,alp1,binsize=0.05,plot=False)
    z,mu0m,yerr = bp.bin_scatter(redshift,mu0,binsize=0.05,plot=False)
    z,mu1m,yerr = bp.bin_scatter(redshift,mu1,binsize=0.05,plot=False)
    z,sgm0m,yerr = bp.bin_scatter(redshift,sgm0,binsize=0.05,plot=False)
    z,sgm1m,yerr = bp.bin_scatter(redshift,sgm1,binsize=0.05,plot=False)
    
    alpha = np.array(zip(alp0m,alp1m))
    mu = np.array(zip(mu0m,mu1m))
    sigma = np.array(zip(sgm0m,sgm1m))
    
    x=np.arange(-0.5,2.5,0.01)
    pl.figure(figsize=(15,15))
    nn = len(z)
    for i in range(nn):
        pl.subplot(3,3,i+1)
        gmm.ecgmmplot(x,alpha[i],mu[i],sigma[i])
        pl.title('zbin: '+str(round(z[i],3))+' +/- 0.025')
        pl.xlabel('g - r')
    pl.figtext(0.5,0.95,'detmodel color')
    pl.savefig(catdir+'des_vvds_gmr_detmodel_distribution_obs_zbin.png')
    pl.close()
    
    #------rmi ridgeline distribution ------

    redshift = resRI[9]
    alp0 = resRI[0]
    alp1 = resRI[1]
    mu0 = resRI[2]
    mu1 = resRI[3]
    sgm0 = resRI[4]
    sgm1 = resRI[5]
    z,alp0m,yerr = bp.bin_scatter(redshift,alp0,binsize=0.05,plot=False)
    z,alp1m,yerr = bp.bin_scatter(redshift,alp1,binsize=0.05,plot=False)
    z,mu0m,yerr = bp.bin_scatter(redshift,mu0,binsize=0.05,plot=False)
    z,mu1m,yerr = bp.bin_scatter(redshift,mu1,binsize=0.05,plot=False)
    z,sgm0m,yerr = bp.bin_scatter(redshift,sgm0,binsize=0.05,plot=False)
    z,sgm1m,yerr = bp.bin_scatter(redshift,sgm1,binsize=0.05,plot=False)
    
    alpha = np.array(zip(alp0m,alp1m))
    mu = np.array(zip(mu0m,mu1m))
    sigma = np.array(zip(sgm0m,sgm1m))
    
    x=np.arange(-0.5,1.5,0.01)
    pl.figure(figsize=(15,15))
    nn = len(z)
    for i in range(nn):
        pl.subplot(3,3,i+1)
        gmm.ecgmmplot(x,alpha[i],mu[i],sigma[i])
        pl.title('zbin: '+str(round(z[i],3))+' +/- 0.025')
        pl.xlabel('r - i')
    pl.figtext(0.5,0.95,'detmodel color')
    pl.savefig(catdir+'des_vvds_rmi_detmodel_distribution_obs_zbin.png')
    pl.close()
 

    #------imz ridgeline distribution ------

    redshift = resIZ[9]
    alp0 = resIZ[0]
    alp1 = resIZ[1]
    mu0 = resIZ[2]
    mu1 = resIZ[3]
    sgm0 = resIZ[4]
    sgm1 = resIZ[5]
    z,alp0m,yerr = bp.bin_scatter(redshift,alp0,binsize=0.05,plot=False)
    z,alp1m,yerr = bp.bin_scatter(redshift,alp1,binsize=0.05,plot=False)
    z,mu0m,yerr = bp.bin_scatter(redshift,mu0,binsize=0.05,plot=False)
    z,mu1m,yerr = bp.bin_scatter(redshift,mu1,binsize=0.05,plot=False)
    z,sgm0m,yerr = bp.bin_scatter(redshift,sgm0,binsize=0.05,plot=False)
    z,sgm1m,yerr = bp.bin_scatter(redshift,sgm1,binsize=0.05,plot=False)
    
    alpha = np.array(zip(alp0m,alp1m))
    mu = np.array(zip(mu0m,mu1m))
    sigma = np.array(zip(sgm0m,sgm1m))
    
    x=np.arange(-0.5,1.5,0.01)
    pl.figure(figsize=(15,15))
    nn = len(z)
    for i in range(nn):
        pl.subplot(3,3,i+1)
        gmm.ecgmmplot(x,alpha[i],mu[i],sigma[i])
        pl.title('zbin: '+str(round(z[i],3))+' +/- 0.025')
        pl.xlabel('i - z')
    pl.figtext(0.5,0.95,'detmodel color')
    pl.savefig(catdir+'des_vvds_rmi_detmodel_distribution_obs_zbin.png')
    pl.close()
 
    #-----check color -----
    idx_iz = (z >= 0.75)*(z < 1.15)
    idx_gr = z < 0.4
    idx_ri = (z >= 0.4)*(z < 0.75 )
    pl.figure(figsize=(10,10))
    pl.subplot(2,2,1)
    #bp.bin_scatter(z,gmr,nbins=20,fmt='bo',scatter=True)
    pl.plot(z[idx_gr],gmr[idx_gr],'b.')
    pl.xlabel('z')
    pl.ylabel('g-r')
    pl.ylim(-0.5,3)
    pl.plot(z[idx_gr],zgmr(z[idx_gr]),'r.')
    pl.subplot(2,2,2)
    #bp.bin_scatter(z,rmi,nbins=20,fmt='bo',scatter=True)
    pl.plot(z[idx_ri],rmi[idx_ri],'b.')
    pl.xlabel('z')
    pl.ylabel('r-i')
    pl.ylim(-0.5,2.)
    pl.plot(z[idx_ri],zrmi(z[idx_ri]),'r.')
    pl.subplot(2,2,3)
    #bp.bin_scatter(z,imz,nbins=20,fmt='bo',scatter=True)
    pl.plot(z[idx_iz],imz[idx_iz],'b.')
    pl.xlabel('z')
    pl.ylabel('i-z')
    pl.ylim(-0.5,2.)
    pl.plot(z[idx_iz],zimz(z[idx_iz]),'r.')

    pl.figtext(0.45,0.95,'MAG_AUTO')
    pl.savefig(catdir+'des_vvds_color_mag_auto.png')

    
    pl.figtext(0.45,0.95,'MAG_DETMODEL')
    pl.savefig(catdir+'des_vvds_color_mag_detmodel.png')

    #---check the same color using sdss coadd ---
    b=pf.getdata('/home/jghao/research/photoz_nnp/data/coadd_photoz_training_ribarmar_jghao.fit')
    gmr = b.DERED_G - b.DERED_R
    rmi = b.DERED_R - b.DERED_I
    imz = b.DERED_I - b.DERED_Z
    z = b.Z

    idx_iz = (z >= 0.75)*(z < 1.15)
    idx_gr = z < 0.4
    idx_ri = (z >= 0.4)*(z < 0.75 )
    pl.figure(figsize=(10,10))
    pl.subplot(2,2,1)
    #bp.bin_scatter(z,gmr,nbins=20,fmt='bo',scatter=True)
    pl.plot(z[idx_gr],gmr[idx_gr],'b,')
    pl.xlabel('z')
    pl.ylabel('g-r')
    pl.ylim(-0.5,3)
    pl.plot(z[idx_gr],zgmr(z[idx_gr]),'r.')
    pl.subplot(2,2,2)
    #bp.bin_scatter(z,rmi,nbins=20,fmt='bo',scatter=True)
    pl.plot(z[idx_ri],rmi[idx_ri],'b,')
    pl.xlabel('z')
    pl.ylabel('r-i')
    pl.ylim(-0.5,2.)
    pl.plot(z[idx_ri],zrmi(z[idx_ri]),'r.')
    pl.subplot(2,2,3)
    #bp.bin_scatter(z,imz,nbins=20,fmt='bo',scatter=True)
    pl.plot(z[idx_iz],imz[idx_iz],'b,')
    pl.xlabel('z')
    pl.ylabel('i-z')
    pl.ylim(-0.5,2.)
    pl.plot(z[idx_iz],zimz(z[idx_iz]),'r.')

    pl.figtext(0.5,0.95,'SDSS Coadd dered Model_MAG')
    pl.savefig(catdir+'SDSS_coadd_dered_model_mag.png')
