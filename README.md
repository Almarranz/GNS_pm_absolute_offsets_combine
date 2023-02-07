# GNS_pm_absolute_offsets_combine
Alignment with Gaia using lists of combained chips in order to get more Gaia´ stars

0. gns1_B_lists.py. Combines GNS1 H and d ks lists. IN this version we do not
cutout the overlapping fields. The idea is to use more Gaia stars for the 
GNS2 fiel and see if that improves the alignemnt

1. gaia_refstars_offsets_gns_pixels.py. Give Gaia offsets cordinates and move
then to the corresponding epoch. Then uses astroaling to alingn the *gns xy 
coordenates* with the Gaia offsets.
>NOTE1:We use foregroud stars in GNS1 for m
>matching Gaia stars and use the matching stars in astroalign. If we use 
>astroalign directly with the whole lists (fore and background), aa cant find
>the transformation. The same thing we have to do with GNS2, but, since GNS2 
>only has H band, we have to match with GNS1 for foreground, so we are limited
>to the Gaia star overlapping with GNS1. Once aa has found the tranformation
>for GNS2 we aplay it to the whole GNS2 lists.

2. IDL_gaia_offset_gns_xy_align.pro Align GNS to Gaia with a polinomal of degree
1 or 2

3. IDL_GaiaRf_align_gns_offsets.pro Find matching stars and substract for pm.









# # GNS_pm_absolute_offsets
# ## Obtain absolute proper motions from GNS1 and GNS2, using GAIA as reference frame

# ### GRF (ICRS): OFFSETS(xy GNS and  deproyected coordinates GAIA)
# 0. gns1_B_lists.py. Combines GNS1 H and d ks lists

# 1. gaia_refstars_offsets_gns_pixels.py. Give Gaia offsets cordinates and move
# then to the corresponding epoch. Then uses astroaling to alingn the *gns xy 
# coordenates* with the Gaia offsets.

# 2. IDL_gaia_offset_gns_xy_align.pro Align GNS to Gaia with a polinomal of degree
# 1 or 2

# 2.1 (A) IDL_BSgaiaRF_offsets_align_gns.pro. genareted lists using bootstrapping
# estimate for the alignemnt uncertainties.
# 2.1 (B) IDL_JKgaiaRF_offsets_align_gns. genareted lists using jackinnifing
# estimate for the alignemnt uncertainties.

# 2.2 (A). bs_uncert_offsets.py. Give the ra and dec bs uncertainties for each star
# 2.2 (B). jk_uncert_offsets.py. Give the ra and dec jk uncertainties for each star


# 3. IDL_GaiaRf_align_gns_offsets.pro Find matching stars and substract for pm.
# Also incorporates the aignemnet errors to the pm ones.
# >Note: this can be easyly made in python.

# 3.1 gaia_stars_offsets_comp.py. Computes residuals with Gaia catalog. The point
# is to get rid of the 3sigma stars and repit the alignment without then.
# 4.pm_analysis_offsets_radec.py. Find and plot clusters using dbscan.