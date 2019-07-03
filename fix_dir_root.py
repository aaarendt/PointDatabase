#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:50:13 2019

@author: ben
"""

import h5py
old_dir='/Volumes/insar5/gmap/OIB_data/'
new_dir='/Data/'


with h5py.File('/Data/glas/GL/rel_634/GeoIndex.h5','r') as h5f:
    N=0
    attr_name='file_%d' % N
    while attr_name in h5f['index'].attrs:
        #h5f['index'].attrs[attr_name]=h5f['index'].attrs[attr_name].replace(old_dir, new_dir)
        print(h5f['index'].attrs[attr_name])
        N += 1
        attr_name='file_%d' % N
h5f.close()
