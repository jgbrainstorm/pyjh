#! /usr/bin/env python
# This code calculate the ridgeline colors using ecgmm for the des mock----
# J. Hao @ Fermilab, 4/1/2012

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


def ecgmmRidge(ra_c=None, dec_c=None,photoz_c=None,r200_c=None,m200_c=None,ra=None, dec=None,color=None,colorErr=None,mag=None,candidateIdx=None):
    #--define some quantity to be returned ----
    rac = ra_c[candidateIdx]
    decc = dec_c[candidateIdx]
    photozc = photoz_c[candidateIdx]
    r200c = r200_c[candidateIdx]
    m200c = m200_c[candidateIdx]
    ok = (color>=-1) * (color <= 5.)
    color = color[ok]
    colorErr = colorErr[ok]
    #colorErr = np.zeros(len(color))
    mag = mag[ok]
    ra = ra[ok]
    dec = dec[ok]
    BCGalpha0 = []
    BCGalpha1 = []
    BCGmu0 = []
    BCGmu1 = []
    BCGsigma0 = []
    BCGsigma1 = []
    BCGntot = []
    BCGamp = []
    BCGaic1 = []
    BCGaic2 = []
    BCGphotoz=[]
    BCGm200c = []
    #-----------------------------------------
    Ncandidates = len(photozc)
    ridgeZ = np.zeros(Ncandidates)
    depth = 10
    h=es.htm.HTM(depth)
    Cosmo = es.cosmology.Cosmo(h=0.7)
    DA=Cosmo.Da(0,photozc)
    srad=np.rad2deg(r200c/DA)
    m1,m2,d12 = h.match(rac,decc,ra,dec,srad,maxmatch=5000)
    r12=np.deg2rad(d12)*DA[m1]
    indices=(mag[m2]<=limz(photozc[m1])) # no bcg assumed
    m1 = m1[indices]
    m2 = m2[indices]
    h,rev = es.stat.histogram(m1, binsize=1, rev=True)
    startTime=time.time()
    for i in range(h.size):
        if rev[i] != rev[i+1]:
            print i
            indx = rev[ rev[i]:rev[i+1]]
            alpha=np.array([0.5,0.5])
            mu=np.array([sts.scoreatpercentile(color[m2[indx]],per=70),sts.scoreatpercentile(color[m2[indx]],per=40)])
            sigma=np.array([0.04,0.3])
            aic2 = gmm.aic_ecgmm(color[m2[indx]],colorErr[m2[indx]],alpha,mu,sigma)           
            aic1 = gmm.wstat(color[m2[indx]],colorErr[m2[indx]])[2]
            if aic2 < aic1:
                srt=np.argsort(sigma)
                BCGalpha0.append(alpha[srt[0]])
                BCGalpha1.append(alpha[srt[1]])
                BCGmu0.append(mu[srt[0]])
                BCGmu1.append(mu[srt[1]])
                BCGsigma0.append(sigma[srt[0]])
                BCGsigma1.append(sigma[srt[1]])
                BCGaic1.append(aic1)
                BCGaic2.append(aic2)
                BCGamp.append(len(indx)*alpha[srt[0]])
                BCGphotoz.append(photozc[m1[indx[0]]])
                BCGm200c.append(m200c[m1[indx[0]]])
                print aic2, aic1
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)
    return np.array(BCGalpha0),np.array(BCGalpha1),np.array(BCGmu0),np.array(BCGmu1),np.array(BCGsigma0),np.array(BCGsigma1),np.array(BCGaic1),np.array(BCGaic2),np.array(BCGamp),np.array(BCGphotoz),np.array(BCGm200c)





#------catalogs -------------------

if __name__ == '__main__':
    version='v4.00'
    catdir = '/home/jghao/research/data/des_mock/'+version+'/'
    #galtype = 'truth'
    galtype = 'obs'
    halo = pf.getdata(catdir+'haloCat/DES_Mock_'+version+'_halos.fit')
    if galtype == 'truth':
        gal = pf.getdata(catdir+'truthCat/DES_Mock_'+version+'_Baseline_truth_00.fit')
    else:
        gal = pf.getdata(catdir+'obsCat/DES_Mock_'+version+'_Baseline_04.fit')
        
    #----choose halos by mass threshold ---
    mass_threshold = 1e14
    halo = halo[halo.field('m200')>= mass_threshold]

    #---centers from halo catalog -------
    ra_c = halo.field('ra')
    dec_c = halo.field('dec')
    photoz_c = halo.field('z')
    r200_c = halo.field('r200')
    m200_c = halo.field('m200')
    #----catalog to be searched -------
    if galtype == 'truth':
        ra = gal.field('ra')
        dec = gal.field('dec')
        mag=gal.field('omag')[:,3]
        ngal = len(mag)
        gmr=gal.field('omag')[:,0] - gal.field('omag')[:,1]
        gmrErr=np.sqrt(gal.field('omagerr')[:,0]**2+gal.field('omagerr')[:,1]**2)
        rmi=gal.field('omag')[:,1] - gal.field('omag')[:,2]
        rmiErr=np.sqrt(gal.field('omagerr')[:,1]**2+gal.field('omagerr')[:,2]**2)
        imz=gal.field('omag')[:,2] - gal.field('omag')[:,3]
        imzErr=np.sqrt(gal.field('omagerr')[:,2]**2+gal.field('omagerr')[:,3]**2)
    else:
        ra = gal.field('ra')
        dec = gal.field('dec')
        mag=gal.field('mag_z')
        gmr=gal.field('mag_g') - gal.field('mag_r')
        gmrErr=np.sqrt(gal.field('magerr_g')**2+gal.field('magerr_r')**2)
        rmi=gal.field('mag_r') - gal.field('mag_i')
        rmiErr=np.sqrt(gal.field('magerr_r')**2+gal.field('magerr_i')**2)
        imz=gal.field('mag_i') - gal.field('mag_z')
        imzErr=np.sqrt(gal.field('magerr_i')**2+gal.field('magerr_z')**2)

    #----g - r color -----
    pl.figure(figsize=(10,10))
    pl.subplot(2,2,1)
    idx_gr = photoz_c < 0.4
    resGR= ecgmmRidge(ra_c=ra_c, dec_c=dec_c,photoz_c=photoz_c,r200_c=r200_c,m200_c = m200_c,ra=ra, dec=dec,color=gmr,colorErr=gmrErr,mag=mag,candidateIdx=idx_gr)
    color = resGR[2]
    redshift = resGR[9]
    sgm0 = resGR[4]
    sgm1 = resGR[6]
    zz = np.arange(0.05,0.5,0.1)
    bp.bin_scatter(redshift,color,binsize=0.05,fmt='bo',label='measurement')
    pl.plot(zz,zgmr(zz),'r-',label='color model')
    if galtype == 'truth':
        pl.title('DES Mock '+version+' : Truth')
    else:
        pl.title('DES Mock'+version+': Obs Catalog')
    pl.legend(loc='best')
    pl.xlabel('redshift')
    pl.ylabel('Ridgeline Color: g-r')

# --- r - i color -------
    pl.subplot(2,2,2)
    idx_ri = (photoz_c >= 0.4)*(photoz_c < 0.75 )
    resRI= ecgmmRidge(ra_c=ra_c, dec_c=dec_c,photoz_c=photoz_c,r200_c=r200_c,m200_c = m200_c,ra=ra, dec=dec,color=rmi,colorErr=rmiErr,mag=mag,candidateIdx=idx_ri)
    color = resRI[2]
    redshift = resRI[9]
    sgm0 = resRI[4]
    sgm1 = resRI[6]
    zz = np.arange(0.4,0.9,0.1)
    bp.bin_scatter(redshift,color,binsize=0.05,fmt='bo',label='measurement')
    pl.plot(zz,zrmi(zz),'r-',label='color model')
    if galtype == 'truth':
        pl.title('DES Mock '+version+' : Truth')
    else:
        pl.title('DES Mock '+version+' : Obs Catalog')
    pl.legend(loc='best')
    pl.xlabel('redshift')
    pl.ylabel('Ridgeline Color: r-i')

#-----i - z color -----------
    pl.subplot(2,2,3)
    idx_iz = (photoz_c >= 0.75)*(photoz_c < 1.15)
    resIZ= ecgmmRidge(ra_c=ra_c, dec_c=dec_c,photoz_c=photoz_c,r200_c=r200_c,m200_c = m200_c,ra=ra, dec=dec,color=imz,colorErr=imzErr,mag=mag,candidateIdx=idx_iz)
    color = resIZ[2]
    redshift = resIZ[9]
    sgm0 = resIZ[4]
    sgm1 = resIZ[6]
    zz = np.arange(0.75,1.3,0.1)
    bp.bin_scatter(redshift,color,binsize=0.05,fmt='bo',label='measurement')
    pl.plot(zz,zimz(zz),'r-',label='color model')
    if galtype == 'truth':
        pl.title('DES Mock '+version+' : Truth')
    else:
        pl.title('DES Mock '+version+' : Obs Catalog')
    pl.legend(loc='best')
    pl.xlabel('redshift')
    pl.ylabel('Ridgeline Color: i-z')

#-----Mass - richness color -----------
    pl.subplot(2,2,4)
    mass = np.concatenate([resGR[10],resRI[10],resIZ[10]])
    rich = np.concatenate([resGR[8],resRI[8],resIZ[8]])
    pl.loglog(mass,rich,'b.')
    if galtype == 'truth':
        pl.title('DES Mock '+version+' : Truth')
    else:
        pl.title('DES Mock '+version+' : Obs Catalog')
    pl.legend(loc='best')
    pl.xlabel('M200')
    pl.ylabel('Richness')
    pl.ylim(10,200)
    if galtype == 'truth':
        pl.savefig(catdir+'DES_Mock_'+version+'_ridgeline_mass_rich_truthCat.png')
    else:
        pl.savefig(catdir+'DES_Mock_'+version+'_ridgeline_mass_rich_obsCat.png')
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
        if galtype == 'truth':
            pl.ylabel('DES Mock '+version+' : Truth')
        else:
            pl.ylabel('DES Mock '+version+' : Obs catalog')
    if galtype == 'truth':
        pl.savefig(catdir+'DES_Mock_'+version+'_gmr_distribution_truth.png')
    else:
        pl.savefig(catdir+'DES_mock_'+version+'_gmr_distribution_obs.png')
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
        if galtype == 'truth':
            pl.ylabel('DES Mock '+version+' : Truth')
        else:
            pl.ylabel('DES Mock '+version+' : obs catalog')
    if galtype == 'truth':
        pl.savefig(catdir+'DES_Mock_'+version+'_rmi_distribution_truth.png')
    else:
        pl.savefig(catdir+'DES_Mock_'+version+'_rmi_distribution_obs.png')
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
        if galtype == 'truth':
            pl.ylabel('DES Mock '+version+' : Truth')
        else:
            pl.ylabel('DES Mock '+version+' : obs catalog')
    if galtype == 'truth':
        pl.savefig(catdir+'DES_Mock_'+version+'_imz_distribution_truth.png')
    else:
        pl.savefig(catdir+'DES_Mock_'+version+'_imz_distribution_obs.png')
    pl.close()
 
