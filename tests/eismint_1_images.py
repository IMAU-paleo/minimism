#!/usr/bin/env python

import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("eismint_1/help_fields_grid_ANT.nc")
halfy = int(len(ds.y)/2)
Hi = ds.Hi.isel(time=-1, y = halfy).where((ds.Hi.x > 0) & (ds.Hi.x < 7e5))
plot = Hi.plot.line(x='x')
plt.savefig('eismint_1_hi.png')
