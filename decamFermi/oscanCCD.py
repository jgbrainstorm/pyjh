#! /usr/bin/env python                                                                                                              
import sys
sys.path.append('/home/jghao/fmccd')
if len(sys.argv) == 1:
    print 'syntax: oscanCCD imagefile extension'
else:
    from fermiMCCD import *
    hdu=pf.open(sys.argv[1])
    i=int(sys.argv[2])
    col0=hdu[i].header['biasseca'].split('[')[1].split(']')[0].split(',')[0].split(':')[0]
    col1=hdu[i].header['biasseca'].split('[')[1].split(']')[0].split(',')[0].split(':')[1]
    row0=hdu[i].header['biasseca'].split('[')[1].split(']')[0].split(',')[1].split(':')[0]
    row1=hdu[i].header['biasseca'].split('[')[1].split(']')[0].split(',')[1].split(':')[1]
    oscanA=hdu[i].data[int(row0)-1:int(row1),int(col0)-1:int(col1)]
    mdA=np.median(oscanA,axis=1)
    rowA=np.arange(0,mdA.shape[0])     
    col0=hdu[i].header['biassecb'].split('[')[1].split(']')[0].split(',')[0].split(':')[0]
    col1=hdu[i].header['biassecb'].split('[')[1].split(']')[0].split(',')[0].split(':')[1]
    row0=hdu[i].header['biassecb'].split('[')[1].split(']')[0].split(',')[1].split(':')[0]
    row1=hdu[i].header['biassecb'].split('[')[1].split(']')[0].split(',')[1].split(':')[1]
    oscanB=hdu[i].data[int(row0)-1:int(row1),int(col0)-1:int(col1)]        
    mdB=np.median(oscanB,axis=1)
    rowB=np.arange(0,mdB.shape[0])
    pl.plot(rowA,mdA,'r-',label='A')
    pl.plot(rowB,mdB,'b-',label='B')
    pl.xlabel('# of row')
    pl.ylabel('Median Pixel Value (ADU)')
    pl.title('Overscan of '+hdu[i].header['detpos'])
    pl.ylim(-100,5000)
    pl.legend()
    pl.show()

