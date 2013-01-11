from fits2color import *
from wcs2kmlHao import *
import os



dir = '/home/jghao/research/data/des_realimage/des-google/stripe82_coadd/'
voiletF=dir+'image_g.fits'
redF=dir+'image_z.fits'
greenF=dir+'image_i.fits'
blueF=dir+'image_r.fits'


sliceImg(redF,'z')
sliceImg(greenF,'i')
sliceImg(blueF,'r')
sliceImg(voiletF,'g')


#-----make color image ---
#r,i,z
nimg = 4
for i in range(nimg):
    for j in range(nimg):
        print i,j
        redF = dir+'medium/image_z_'+str(i)+'_'+str(j)+'.fits'
        greenF = dir+'medium/image_i_'+str(i)+'_'+str(j)+'.fits'
        blueF = dir+'medium/image_r_'+str(i)+'_'+str(j)+'.fits'
        size = 0.1
        img=colorImg(redF,greenF,blueF,scale=[0.054*size,0.04*size,0.04*size],nonlinearity=20,smooth=1.2)
        img = ImageEnhance.Brightness(img)
        img = img.enhance(1.7)
        img=ImageEnhance.Contrast(img)
        img=img.enhance(2.5)
        img.save(dir+'medium/riz/image_rgb_'+str(i)+'_'+str(j)+'.png')

#-----make color image ---
#g,r,i
nimg = 4
for i in range(nimg):
    for j in range(nimg):
        print i,j
        redF = dir+'medium/image_i_'+str(i)+'_'+str(j)+'.fits'
        greenF = dir+'medium/image_r_'+str(i)+'_'+str(j)+'.fits'
        blueF = dir+'medium/image_g_'+str(i)+'_'+str(j)+'.fits'
        size = 0.1
        img=colorImg(redF,greenF,blueF,scale=[0.042*size,0.04*size,0.045*size],nonlinearity=10,smooth=1.)
        img = ImageEnhance.Brightness(img)
        img = img.enhance(1.4)
        img=ImageEnhance.Contrast(img)
        img=img.enhance(1.8)
        img.save(dir+'medium/gri/image_rgb_'+str(i)+'_'+str(j)+'.png')



def testcontrast(contrast=None,bright = None,size=None,nonlinearity=None):
    i = 1
    j = 2
    redF = dir+'medium/image_i_'+str(i)+'_'+str(j)+'.fits'
    greenF = dir+'medium/image_r_'+str(i)+'_'+str(j)+'.fits'
    blueF = dir+'medium/image_g_'+str(i)+'_'+str(j)+'.fits'
    img=colorImg(redF,greenF,blueF,scale=[0.042*size,0.040*size,0.045*size],nonlinearity=nonlinearity,smooth=1.)
    img = ImageEnhance.Brightness(img)
    img = img.enhance(bright)
    ehImg=ImageEnhance.Contrast(img)
    newImg=ehImg.enhance(contrast)
    newImg.show()

            
#----combine back --riz----
blank_image = Image.new('RGB',(14000,14000))
nimg = 4
npix=3500
for i in range(nimg):
    for j in range(nimg):
        print i,j
        colimg = dir+'medium/riz/image_rgb_'+str(i)+'_'+str(j)+'.png'
        im = Image.open(colimg)
        blank_image.paste(im,(j*npix,i*npix))
        
 
#----combine back --gri----
blank_image = Image.new('RGB',(14000,14000))
nimg = 4
npix=3500
for i in range(nimg):
    for j in range(nimg):
        print i,j
        colimg = dir+'medium/gri/image_rgb_'+str(i)+'_'+str(j)+'.png'
        im = Image.open(colimg)
        blank_image.paste(im,(j*npix,i*npix))
        blank_image.save('des_coadd_1deg_gri.png')
 


#------kml files---
nimg = 20
for i in range(nimg):
    for j in range(nimg):
        print i,j
        fitsName = dir+'small/image_z_'+str(i)+'_'+str(j)+'.fits'
        imgName=dir+'small/image_rgb_'+str(i)+'_'+str(j)+'.png'
        genkml(fitsfile=fitsName,imagefile=imgName)

make_root_kml(kmlDIR = dir+'small/')
