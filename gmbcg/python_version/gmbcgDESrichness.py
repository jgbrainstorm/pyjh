"""
This is the python implementation of richness remeasurement using the radially weighted ECGMM.
J. Hao @ Fermilab
9/17/2012
"""

from gmbcgFinder import *
import cPickle as p

def BCGidx(bcgRA,bcgDEC,galRA,galDEC):
    depth = 10
    h=es.htm.HTM(depth)
    srad = 7.5e-05  # pixel = 0.27 arcsec, in degrees
    m1,m2,d12 = h.match(bcgRA,bcgDEC,galRA,galDEC,srad)
    return m2

    


#-----setup catalog directories---
BCGCatDir = '/home/jghao/research/data/des_mock/v4.00/clusterCat/gaussian_photoz/'
InputCatDir = '/home/jghao/research/data/des_mock/v4.00/obsCat/'

BCGF = gl.glob(BCGCatDir+'*.fit')
galF = gl.glob(InputCatDir+'*.fit')
galF.sort()
BCGF.sort()
NF = len(galF)

Truth = 0
#----read in file and prepare the input variable -----
for i in range(10):
    i=0
    cat=pf.getdata(galF[i],1)
    hdu = pf.open(BCGF[i])
    bcg = pf.getdata(BCGF[i],1)
    ra=cat.field('ra')
    dec=cat.field('dec')
    photoz=cat.field('PHOTOZ_GAUSSIAN')
    objID=cat.field('id')
    
    mag=cat.field('mag_z')
    gmr=cat.field('mag_g') - cat.field('mag_r')
    gmrErr=np.sqrt(cat.field('magerr_g')**2+cat.field('magerr_r')**2)
    rmi=cat.field('mag_r') - cat.field('mag_i')
    rmiErr=np.sqrt(cat.field('magerr_r')**2+cat.field('magerr_i')**2)
    imz=cat.field('mag_i') - cat.field('mag_z')
    imzErr=np.sqrt(cat.field('magerr_i')**2+cat.field('magerr_z')**2)
    zmy=cat.field('mag_z') - cat.field('mag_y')
    zmyErr=np.sqrt(cat.field('magerr_z')**2+cat.field('magerr_y')**2)

    gr = bcg[bcg.field('photoz') < 0.4]
    
    grCandIdx = BCGidx(gr.field('RA'),gr.field('DEC'),cat.field('RA'),cat.field('DEC'))
    resGR= gmbcgFinder(objID=objID,ra=ra, dec=dec, photoz=photoz,color=gmr,colorErr=gmrErr,mag=mag,bcgCandidateIdx=grCandIdx)

    ri= bcg[(bcg.field('photoz') >= 0.4)*(bcg.field('photoz') < 0.65)]
    riCandIdx = BCGidx(ri.field('RA'),ri.field('DEC'),cat.field('RA'),cat.field('DEC'))
    resRI= gmbcgFinder(objID=objID,ra=ra, dec=dec, photoz=photoz,color=rmi,colorErr=rmiErr,mag=mag,bcgCandidateIdx=riCandIdx)

    bcg.LH = 0.
    idgr = np.arange(len(bcg.LH))[bcg.field('photoz') < 0.4]
    bcg[idgr].LH = resGR[10]
    idri = np.arange(len(bcg.LH))[(bcg.field('photoz') >= 0.4)*(bcg.field('photoz') < 0.65)]
    bcg[idri].LH = resRI[10]
    bcg = bcg[(bcg.field('lh') > 5)*(bcg.field('lh')<500)]
    
"""
BCGIdx, BCGobjID, BCGalpha0, BCGalpha1,BCGmu0, BCGmu1, BCGsigma0,BCGsigma1,BCGaic1,BCGaic2,BCGamp


"""
