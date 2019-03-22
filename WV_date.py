# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:51:44 2019

@author: ben
"""
import numpy as np
import re
from datetime import date
def WV_date(filename):
    date_re=re.compile('\d\d.*_(2\d\d\d)(\d\d)(\d\d)_')
    m=date_re.search(filename)
    if m is None:
        return np.NaN
    return date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    
def WV_year(filename):
    this_date=WV_date(filename)
    this_delta=this_date-date(2000, 1, 1)    
    return 2000+this_delta.days/365.25
    