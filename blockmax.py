# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 15:44:31 2018

@author: ben
"""

import numpy as np
import matplotlib.pyplot as plt
def blockmax(x, y, z, delta, xy0=[0.,0.]):
    yscale=np.ceil((np.nanmax(y)-np.nanmin(y))/delta*1.1)
    zscale=(np.nanmax(z)-np.nanmin(z))*1.1
    xr=np.floor((x-xy0[0])/delta)
    yr=np.floor((y-xy0[1])/delta)
    xyind=xr*yscale+(yr-np.min(yr))
    ind=np.argsort(xyind+(z-np.nanmin(z))/zscale)
    xs=x[ind]
    ys=y[ind]
    zs=z[ind]
    xyind=xyind[ind]
    ux, ix=np.unique(xyind, return_index=True)
    xm=np.zeros_like(ux)+np.NaN
    ym=np.zeros_like(ux)+np.NaN
    zm=np.zeros_like(ux)+np.NaN
    ix=np.concatenate(ix, xyind.size)
    for  count, i0 in enumerate(ix[:-1]):
        iM=ix[count+1]-1
        zm[count]=zs[iM]
        xm[count]=xs[iM]
        ym[count]=ys[iM]
    return xm, ym, zm
