#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 15:35:13 2018

@author: ben
"""

from osgeo import gdal, gdalconst
import numpy as np
import matplotlib.pyplot as plt


class mapData(object):
    def __init__(self):
        self.x=None
        self.y=None
        self.z=None
        self.projection=None
        self.extent=None
        
    def from_geotif(self, file, bands=None, bounds=None, skip=1):
        """
            Read a raster from a DEM file
        """
        ds=gdal.Open(file, gdalconst.GA_ReadOnly)
        GT=ds.GetGeoTransform()
        proj=ds.GetProjection()
        if bands is None:
            n_bands=ds.RasterCount
            bands=np.arange(n_bands, dtype=int)
        if not isinstance(bands, (list, tuple, np.ndarray)):
            bands=[bands]
        # get geolocation info, allocate outputs    
        band=ds.GetRasterBand(1)
        nodataValue=band.GetNoDataValue()
        # ii and jj are the pixel center coordinates.  0,0 in GDAL is the upper-left
        # corner of the first pixel.
        ii=np.arange(0, band.XSize)+0.5
        jj=np.arange(0, band.YSize)-0.5
        x=GT[0]+GT[1]*ii
        y=GT[3]+GT[5]*jj
        if bounds is not None:
            cols = np.where(( x>=bounds[0][0] ) & ( x<= bounds[0][1] ))[0]
            rows = np.where(( y>=bounds[1][0] ) & ( y<= bounds[1][1] ))[0]
        else:
            rows=np.arange(band.YSize, dtype=int)
            cols=np.arange(band.XSize, dtype=int)
        z=list()
        for band_num in bands:
            band=ds.GetRasterBand(int(band_num+1))
            z.append(band.ReadAsArray(int(cols[0]), int(rows[0]), int(cols[-1]-cols[0]+1), int(rows[-1]-rows[0]+1)))
            if skip > 1:
                z[-1]=z[-1][::skip, ::skip]
        if len(bands)==1:
            z=z[0]
        else:
            z=np.stack(z, axis=2)
        ds=None
    
        if skip >1:
            cols=cols[::skip]
            rows=rows[::skip]
        if nodataValue is not None and np.isfinite(nodataValue):
            bad = z==np.array(nodataValue).astype(z.dtype)
            z = np.float64(z)
            z[bad] = np.NaN
        else:
            z = np.float64(z)
        x=x[cols]
        y=y[rows]
        self.x=x
        self.y=y
        self.z=z
        self.projection=proj
        self.extent=[np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)]
        return self
    
    def add_alpha_band(self, alpha=None, nodata_vals=None):
        if alpha is not None:
            self.z=np.concatenate([self.z, alpha], axis=2);
            return
        if nodata_vals is not None:
            alpha=np.ones_like(self.z[:,:,0])
            if hasattr(nodata_vals, 'len') and len(nodata_vals)==3:
                for ii in range(3):
                    alpha[~np.isfinite(self.z[:,:,ii]) | (self.z[:,:,ii]==nodata_vals[ii])]=0
            elif nodata_vals is not None:
                alpha[np.all(~np.isfinite(self.z) | (self.z==nodata_vals), axis=2)]=0
        else:
            alpha=np.any(~np.isfinite(self.z), axis=2)
        if len(alpha.shape)<3:
            alpha.shape=(alpha.shape[0], alpha.shape[1], 1)
        self.z=np.concatenate([self.z, alpha], axis=2)
        return self

    def normalize(self, z0=[0., 255.], z1=[0., 1.], truncate=True, dtype=np.float64):    
        
        self.z=((self.z.astype(np.float64))-z0[0])/(z0[1]-z0[0])*(z1[1]-z1[0])+z1[0]
        if truncate:
            self.z[self.z < z1[0]]=z1[0]
            self.z[self.z > z1[1]]=z1[1]
        self.z=self.z.astype(dtype)
        return self
        
    def show(self, ax=None):
        if ax is None:
            h_im=plt.imshow(self.z, extent=self.extent)
        else:
            h_im=ax.imshow(self.z, extent=self.extent)
        return h_im
        
        
#thefile='/Volumes/ice2/ben/DMS/20180504/DMS_1842637_02317_20180405_12465616.tif'
#M=mapData().from_geotif(thefile, skip=3).normalize().add_alpha_band(nodata_vals=0)
#M.show()