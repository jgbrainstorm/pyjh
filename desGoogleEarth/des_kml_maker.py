#! /usr/bin/env python

#---this is the final one by coming previous pieces --
from fits2color import *
from des_kml_creation import *


def sliceImg(filename,band,dirr):
    # note that the update of wcs not correct yet, only vliad in stripe 82
    nimg = 4
    npix = 14000/4.
    for i in range(nimg):
        for j in range(nimg):
            print i,j
            imghdu= pf.open(filename)
            imghdu[0].data = imghdu[0].data[i*npix:(i+1)*npix,j*npix:(j+1)*npix]
            xdiff = 1750+j*npix - 7000
            ydiff = 1750+i*npix - 7000
            imghdu[0].header.update('CRPIX1',1750)
            imghdu[0].header.update('CRPIX2',1750)
            imghdu[0].header.update('NAXIS1',3500)
            imghdu[0].header.update('NAXIS2',3500)
            crval1 = imghdu[0].header['CRVAL1']
            crval2 = imghdu[0].header['CRVAL2']
            imghdu[0].header.update('CRVAL1',crval1-xdiff*0.27/3600.)
            imghdu[0].header.update('CRVAL2',crval2+ydiff*0.27/3600.)
            imghdu.writeto(dirr+'medium/image_'+band+'_'+str(i)+'_'+str(j)+'.fits')
    return '---done!---'

def testcontrast(dirr,contrast=None,bright = None,size=None,nonlinearity=None):
    i = 1
    j = 2
    redF = dirr+'medium/image_z_'+str(i)+'_'+str(j)+'.fits'
    greenF = dirr+'medium/image_i_'+str(i)+'_'+str(j)+'.fits'
    blueF = dirr+'medium/image_r_'+str(i)+'_'+str(j)+'.fits'
    img=colorImg(redF,greenF,blueF,scale=[0.07*size,0.05*size,0.04*size],nonlinearity=nonlinearity,smooth=1.)
    img = ImageEnhance.Brightness(img)
    img = img.enhance(bright)
    ehImg=ImageEnhance.Contrast(img)
    newImg=ehImg.enhance(contrast)
    #newImg.show()
    newImg.save('test.png')


#----main program --------

if __name__ == '__main__': 
    from des_kml_maker import *
    dirr = '/home/jghao/research/data/des_realimage/des-google/elgordo/'
    imghdr = 'elgordo'

    redF=dirr+imghdr+'_i.fits'
    greenF=dirr+imghdr+'_r.fits'
    #blueF=dirr+'sptw_g.fits'
    infraredF=dirr+imghdr+'_z.fits'

    pngDIR=dirr
    fitsDIR=dirr
    kmlDIR=dirr+imghdr+'kml/'
    pngName = 'des_riz_'+imghdr+'.png'
    kmlName = pngName[0:-4]+'.kml'

    #-----slice image -------
    sliceImg(redF,'i',dirr)
    sliceImg(greenF,'r',dirr)
    #sliceImg(blueF,'g',dirr)
    sliceImg(infraredF,'z',dirr)

    #-----make g,r,z color image ---
    """
    nimg = 4
    for i in range(nimg):
        for j in range(nimg):
            print i,j
            redF = dirr+'medium/image_z_'+str(i)+'_'+str(j)+'.fits'
            greenF = dirr+'medium/image_r_'+str(i)+'_'+str(j)+'.fits'
            blueF = dirr+'medium/image_g_'+str(i)+'_'+str(j)+'.fits'
            size = 0.45
            img=colorImg(redF,greenF,blueF,scale=[0.048*size,0.044*size,0.046*size],nonlinearity=20,smooth=1.)
            img = ImageEnhance.Brightness(img)
            img = img.enhance(1.4)
            img=ImageEnhance.Contrast(img)
            img=img.enhance(1.9)
            img.save(dirr+'medium/image_rgb_'+str(i)+'_'+str(j)+'.png')
    """
    #-----make r,i,z color image ---
    nimg = 4
    for i in range(nimg):
        for j in range(nimg):
            print i,j
            redF = dirr+'medium/image_z_'+str(i)+'_'+str(j)+'.fits'
            greenF = dirr+'medium/image_i_'+str(i)+'_'+str(j)+'.fits'
            blueF = dirr+'medium/image_r_'+str(i)+'_'+str(j)+'.fits'
            size = 0.25
            img=colorImg(redF,greenF,blueF,scale=[0.07*size,0.05*size,0.04*size],nonlinearity=20,smooth=1.)
            img = ImageEnhance.Brightness(img)
            img = img.enhance(1.4)
            img=ImageEnhance.Contrast(img)
            img=img.enhance(1.9)
            img.save(dirr+'medium/image_rgb_'+str(i)+'_'+str(j)+'.png')



    #----combine back to a large color image-----
    blank_image = Image.new('RGB',(14000,14000))
    nimg = 4
    npix=3500
    for i in range(nimg):
        for j in range(nimg):
            print i,j
            colimg = dirr+'medium/image_rgb_'+str(i)+'_'+str(j)+'.png'
            im = Image.open(colimg)
            blank_image.paste(im,(j*npix,i*npix))
            blank_image.save(dirr+pngName)
 
    # ---- make kml regioneated images -------

    step2=make_kml_image(pngDIR,fitsDIR,kmlDIR)

    # ---- make final kml ------
    step3=make_root_kml(kmlDIR,kmlName)

