#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 21:14:09 2021

@author: jonathan
"""

import os
import numpy as np
from netCDF4 import Dataset

test_folder = "/data2/jonathan/work/code/educational_VMEC/test"

inout = "out" # "in" or "out"
iteration = 7
runId = "BETA_5_ICUR_5K"





ref_fname = None
if inout == "out":
    ref_fname = os.path.join(test_folder, "vac_ref", "vacout_ref_%s_%06d.nc"%(runId, iteration))
else:
    ref_fname = os.path.join(test_folder, "vac_ref", "vacin_%s_%06d.nc"%(runId, iteration))
tst_fname = os.path.join(test_folder, "vac", "vac%s_%s_%06d.nc"%(inout, runId, iteration))

print("ref: "+ref_fname)
print("tst: "+tst_fname)

if not os.path.isfile(ref_fname):
    raise RuntimeError("reference %s not found"%(ref_fname,))
if not os.path.isfile(tst_fname):
    raise RuntimeError("test %s not found"%(tst_fname,))

ref_data = {}
d = Dataset(ref_fname, "r")
for key in d.variables:
    ref_data[key] = d[key][()]
d.close()

tst_data = {}
d = Dataset(tst_fname, "r")
for key in d.variables:
    tst_data[key] = d[key][()]
d.close()

# compare data
for key in ref_data:
    if not key in tst_data:
        raise RuntimeError("key %s not found in %s", key, tst_fname)
    
    r = ref_data[key]
    t = tst_data[key]
    
    if (np.shape(r) != np.shape(t)):
        raise RuntimeError("shape mismatch in %s: ref has %s, tst has %s"%(key, str(np.shape(r)), str(np.shape(t))))
    
    try:    
        d = np.sum(np.abs(r-t))
        if d != 0.0:
            print("value mismatch in %s: %g"%(key, d))
    except:
        pass
