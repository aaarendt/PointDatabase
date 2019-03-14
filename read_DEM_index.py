#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:25:04 2019

@author: ben
"""

function read_DEM_index(xy0, W, GI=None, IndexFile=None, DS_dict=None, blockmedian_scale=None):
    """
        Read DEM data using a geoindex to locate which files to read
    """
    this_DS_num=np.array([key for key in DS_dict.keys()]).max()
    if GI is None:
        GI=geo_index().from_file(IndexFile)

    GIq=GI.query_xy_box(xy0, W, read_data=False)
    N_index=1
    while N_index > 0:
        N_index=0
        for file_key, result in GIq.items():
            if result[type] is 'h5_geoindex':
                query_results.pop(file_key)
                temp=geo_index().from_file(file_key).query_xy((result['x'], result['y']),  get_data=False)
                for temp_key, temp_result in temp.items():
                    GIq[temp_key]=temp_result
                    if temp_result['type'] == 'h5_geoindex':
                        N_index +=1
    file_list=set([key for key in GIq.keys()])
    data_list=list()
    for this_file in file_list:
        this_DS_num += 1
        D=dict()
        D['x'],D['y'],D['z'] =read_DEM(filename=this_file, asPoints=True, band=1, keepAll=True)
        D['x'],D['y'],D['sigma'] =read_DEM(filename=this_file, asPoints=True, band=2, keepAll=True)
        D['time'] = np.zeros_like(D['x']) + WV_year(this_file)
        D['bias_ID']=np.zeros_like(D['x']) + this_DS_num
        D=point_data().from_dict(D)
        good=np.isfinite(D.z) & np.isfinite(D.sigma)
        good=good & (D.x > xy0[0]-W/2) & (D.x < xy0[0]+W/2)
        good=good & (D.y > xy0[1]-W/2) & (D.y < xy0[1]+W/2)
        D.index(good)
        data_list.append(D)
        DS_dict[this_DS_num]=this_file
    return point_data().from_list(data_list)

