===============================
Cosmic Variance Calculator
===============================

Package to calculate cosmic variance in pencil-beam surveys
---------------------------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/cosmic_variance.svg
        :target: https://pypi.python.org/pypi/cosmic_variance

.. image:: https://img.shields.io/travis/astrockragh/cosmic_variance.svg
        :target: https://travis-ci.com/astrockragh/cosmic_variance

.. image:: https://readthedocs.org/projects/cosmic-variance/badge/?version=latest
        :target: https://cosmic-variance.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


Python package based on the IDL code released with the Cosmic Variance Cookbook of Moster et al. (2010)

The code is based on galaxy stellar mass bins (as described in https://arxiv.org/pdf/1001.1737.pdf), scaled to dark matter cosmic variance (as described in https://arxiv.org/pdf/astro-ph/0109130.pdf). 

This is significantly more useful than dark matter - only variance, since the empirical galaxy variance is significantly higher.

Free software: MIT license

Install and Use
-------------------

To install the package, simply run:

.. code-block:: bash

        pip install cosmic-variance

Then in your script/notebook, import the package as:

.. code-block:: python

        import cosmic_variance as cv

The main use of the package is through the get_cv function, which takes in a rectangular survey geometry with side lengths side1 and side2 (in degrees), and an array of redshift bin edges, and returns a pandas dataframe with the cosmic variance for 0.5 dex galaxy stellar mass bins for each redshift bin.

.. code-block:: python

        import cosmic_variance as cv
        import numpy as np
        # Example of using the main function, get_cv to calculate
        # cosmic variance for a single JWST pointing

        #### these arguments are required ####
        side1 = 2.2/60. # /60 to convert from arcmin to degrees
        side2 = 2*2.2/60. # /60 to convert from arcmin to degrees
        zarray = np.array([7,8,9,11,13]) # redshift bin edges

        #### these arguments are optional ####
        name = 'JWST' # name of the survey, if provided, the output file will be saved as dfs/{name}.csv along with a meta file.
        # Default is None, in which case the output will not be saved

        acc = 'low' # accuracy of the calculation, 'low' or 'high, low is default, faster and sufficient for almost all applications

        verbose = False # if True, will print out the progress of the calculation, default is False

        #If you want to use a different cosmology, you can specify it by the following in the get_cv call
        # OmegaM = 0.308, OmegaL = 0.692, OmegaBaryon = 0.022/(0.678)**2 sigma8 = 0.82, ns = 0.96, h = 0.678

        cv_df = cv.get_cv(side1, side2, zarray, name = name, acc=acc, verbose = verbose)

This will calculate the cosmic variance for a 2.2 arcmin x 4.4 arcmin survey in redshifts bin [7, 8], [8,9], [9,11], [11,13] and save the output.