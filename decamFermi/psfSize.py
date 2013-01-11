#! /usr/bin/env python
import numpy as np
import pyfits as pf
import scipy.ndimage as nd
import pylab as pl
import sys


def getxy(datas):
    ok=datas>=datas.mean()+1.5*datas.std()
    good=nd.binary_opening(ok,structure=np.ones((6,6)))    
    datagood=datas*good
    structuring_element=np.ones((3,3))
    segmentation,segments=nd.label(good,structure=structuring_element)    
    coords=nd.center_of_mass(datagood,segmentation,range(1,segments+1))
    xcoords=np.array([x[1] for x in coords])
    ycoords=np.array([x[0] for x in coords])
    return xcoords,ycoords

def fwhm(data,obj,nobj):
    fwhmres = np.zeros(nobj)
    for i in range(nobj):
        fwhmres[i] = data[obj[i]]

def findbstr(filename):
    """
    find the bright stars on the image
    """
    data = pf.getdata(filename)
    hdr = pf.getheader(filename)
    saturate = hdr['saturate']
    bsIDX = (data >= 0.5*saturate)* (data <= 0.8*saturate)
    good=nd.binary_opening(bsIDX,structure=np.ones((3,3)))  
    objData = data*good
    seg,nseg=nd.label(good,structure=np.ones((3,3)))  
    coords=nd.center_of_mass(objData,seg,range(1,nseg+1))
    xcoord=np.array([x[1] for x in coords])
    ycoord=np.array([x[0] for x in coords])
    obj = nd.find_objects(seg, nseg)
    

dir = '/home/jghao/research/data/des_optics_psf/dc6b_image/'

imgname = dir+'decam-34-0-r-0_01.fits'
bkgname = dir+'decam-34-0-r-0_01_bkg.fits'
