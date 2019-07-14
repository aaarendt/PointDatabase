# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 21:40:42 2018

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
def in_axes(x, y, ax=plt.gca()):
    xl=ax.get_xlim()
    yl=ax.get_ylim()
    print(xl)
    print(yl)
    return np.where(np.logical_and(np.logical_and(x>=xl[0], x<=xl[1]), np.logical_and(y>=yl[0], y<=yl[1])))
    