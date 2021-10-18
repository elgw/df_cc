# Correction for shift and low order chromatic aberrations

This repository contains code to correct for shifts and chromatic aberrations
between wide field microscopy channels. The implementation is based on a
[paper](https://doi.org/10.1046/j.1365-2818.2000.00754.x) by Michal Kozubek,
although any bugs or misconceptions are likely to be mine.

Typically these functions are used on calibration beads to set up
a mapping between a reference channel (pick the best one) and the other
channels.

These files are slightly out of context here as they are part of a bigger
toolbox (DOTTER) that might be released at some time point. To use them with
MATLAB, it might be a good idea to add the corresponding folder to the path,
i.e.

``` MATLAB
addpath('~/where/i/keep/df_cc/')
```

# Files
## Create a correction (`.cc`) file:
 - `df_cc` -- Cluster dots to find correspondences, and find correction.
Uses the following functions:
 - `df_cc_cluster` -- find correspondences between points.
 - `df_cc_create` -- find the transformation between two or more channels.
 - `df_cc_view` -- metrics about the correction results.

## Correct dots
 - `df_cc_apply_dots` -- transform dots to match a reference channels

## Correct images:
 - `df_cc_apply_image` -- correct a single image
 - `df_cc_apply_image_folder` -- correct all images in a folder

## Demos/tests
 - `df_cc_demo` -- demo using a constant shift / dot correction.
 - `df_cc_ut` -- unit tests

## Utility files
These files are used by the other scripts:
 - `df_cc_eudist` euclidean distance
 - `df_cc_poly2mat` -- create matrix for finding coefficients.

## Depreciated
 - `df_cc_apply_m` -- Correct `M.dots` in a DOTTER `.NM` file.
 - `df_cc_apply_n` -- Correct all `N.userDots` in a DOTTER `.NM` file.

# TODO
 - [ ] Use a polynomial also in the axial direction?
 - [ ] Add some more demos and demo images.
