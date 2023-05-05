# Rewritten by Christian Kragh Jespersen, @astrockragh, to work in Python 3.8 (Jan-2023)
# Dependencies are numpy, scipy, and pandas
# If anyone (including future me) wants to optimize for speed, the bottleneck is the intpk4 integral


import numpy as np
import pandas as pd
from .quickcv import quickcv
import os

#### Sample parameter ####

# ########### Input parameters ###########
# s1 = 3.4/60 #side 1, in degrees 
# s2 = 10.2/60 #side 2, in degrees
# zarr = np.array([3,4,5,6.5]) #redshift bin edges
# ############################################

# ########### Cosmology ###########
# OmegaM = 0.308
# OmegaL = 0.692
# sigma8 = 0.82
# acc = 'low' # low takes about 20 seconds per calculation, high takes about 8 minutes per calculation
############################################

# The bias values for 7.0, 7.5, and 8.0 were incorrectly extrapolated, here they are set to be the same as for 8.5 in this version by astrockragh 19/11-2022

b070=0.062
b170=2.59
b270=1.025

b075=0.062
b175=2.59
b275=1.025

b080=0.062
b180=2.59
b280=1.025

b085=0.062
b185=2.59
b285=1.025

b090=0.074
b190=2.58
b290=1.039

b095=0.042
b195=3.17
b295=1.147

b0100=0.053
b1100=3.07
b2100=1.225

b0105=0.069
b1105=3.19
b2105=1.269

b0110=0.173
b1110=2.89
b2110=1.438

def get_cv(side1, side2, zarr, name = None, dz = None, acc = 'low', OmegaM = 0.308, OmegaL = 0.692, OmegaBaryon = 0.022/(0.678)**2, sigma8 = 0.82, ns = 0.96, h = 0.678, verbose = False):
   ''' Method to calculate the cosmic variance for a given survey area and redshift bins.
      Give side lengths in degrees
      side1: side length 1 [degrees]
      side2: side length 2 [degrees]
      zarr: array of redshift bin edges
      name: name of the output file. File will be saved at dfs/name.csv if name is not None
      dz: redshift bin size. If None, the redshift bin size will be the difference between the redshift bin edges. If dz is given, zarr is the redshift bin centers.
      acc: accuracy of the calculation. 'low' or 'high'. 'low' is faster (20 sec per calculation), 'high' is more accurate (8 mins per calculation)
      OmegaM, OmegaL, OmegaBaryon, sigma8, ns, h: cosmological parameters
      verbose: whether to print out the time taken for the integral
      
      Output is dataframe with columns:
      zmid: redshift bin center
      dz: redshift bin size
      cv_dm: cosmic variance for dark matter
      cv_xx: cosmic variance for x.x - 0.25 < log(M_*/M_{sun}) < x.x + 0.25

      if given name, the dataframe will be saved to dfs/name.csv, along with a df with metadata, dfs/name_meta.csv so that the runs don't have to be repeated
      '''

   metadict = {'side1': side1,
            'side2': side2,
            'zarr': np.array(zarr),
            'dz': dz,
            'acc': acc,
            'OmegaM': OmegaM,
            'OmegaL': OmegaL,
            'OmegaBaryon': OmegaBaryon, 
            'sigma8': sigma8,
            'ns': ns,
            'h': h,
            'verbose': verbose}

   df_meta = pd.DataFrame(metadict)
   if dz is None:
      nz = len(zarr)-1
   else:
      nz = len(zarr)

   columns = ['zmid', 'dz', 'cv_dm', 'cv_70', "cv_75", "cv_80", "cv_85", "cv_90", "cv_95", "cv_100", "cv_105", "cv_110"]
   df_values = pd.DataFrame(index = np.arange(nz), columns = columns)

   for i in range(nz):
      if dz is None:
         dz      = zarr[i+1]-zarr[i]
         zmid    = zarr[i+1]+zarr[i]
         zmid    = zmid/2
      else:
         if i == 0:
            print('dz is given, zarr is redshift bin centers')
         zmid    = zarr[i]
      
      s       = quickcv.quickcv(side1, side2, np.array([zmid]), np.array([dz]), OmegaM = OmegaM, OmegaL = OmegaL, OmegaBaryon = OmegaBaryon, \
                                sigma8 = sigma8, ns = ns, h = h, acc = acc, verbose = verbose)
      
      ## These are the linear bias values
      bias70  = b070  * (zmid+1)**b170  + b270
      bias75  = b075  * (zmid+1)**b175  + b275
      bias80  = b080  * (zmid+1)**b180  + b280
      bias85  = b085  * (zmid+1)**b185  + b285
      bias90  = b090  * (zmid+1)**b190  + b290
      bias95  = b095  * (zmid+1)**b195  + b295
      bias100 = b0100 * (zmid+1)**b1100 + b2100
      bias105 = b0105 * (zmid+1)**b1105 + b2105
      bias110 = b0110 * (zmid+1)**b1110 + b2110

      ## All gg variables are galaxy-galaxy cosmic variance
      sgg70  = bias70  * s
      sgg75  = bias75  * s
      sgg80  = bias80  * s
      sgg85  = bias85  * s
      sgg90  = bias90  * s
      sgg95  = bias95  * s
      sgg100 = bias100 * s
      sgg105 = bias105 * s
      sgg110 = bias110 * s
      print(zmid, dz, s, sgg70, sgg75, sgg80, sgg85, sgg90, sgg95, sgg100, sgg105, sgg110)
      df_values.loc[i] = [zmid, dz, s, sgg70, sgg75, sgg80, sgg85, sgg90, sgg95, sgg100, sgg105, sgg110]
      #prints the mean redshift, redshift bins size, sigma_dm and the sigma_gg for nine stellar mass bins

   # Save the dataframes to csv files so that we don't have to run the code again
   if name:
      if not os.path.isdir('dfs'):
         os.mkdir('dfs')
      print('Saving dataframes of values to dfs/' + name )
      print(df_values)
      df_values.to_csv('dfs/' + name + '.csv', index = False)
      df_meta.to_csv('dfs/' + name + '_meta' + '.csv', index = False)
   return df_values

## Sample use
# if __name__ == '__main__':
#    get_cv(3.4/60, 10.2/60, np.array([3,4,5,6.5]), name = 'test')