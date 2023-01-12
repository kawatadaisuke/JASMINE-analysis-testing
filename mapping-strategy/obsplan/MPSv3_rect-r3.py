
# MPSv3_rect-r3.py
# finer dithering version of MPSv3r2_rect-r2.py which covers
# the whole JASMINE region similar number of times and at the different 
# position in the detector.

import numpy as np
from telescope_baseline.mapping.aperture import lb_detector_unit, four_square_convexes, inout_four_square_convexes, ang_detector_unit, lb2ang, ang2lb
import h5py

from matplotlib import patches

if __name__ == "__main__":
    import pkg_resources           
    from telescope_baseline.mapping.read_catalog import read_jasmine_targets
    from telescope_baseline.mapping.mapset import ditheringmap, inout_convexesset
    from telescope_baseline.mapping.plot_mapping import plot_targets, plot_n_targets, hist_n_targets, plot_ae_targets, hist_ae_targets, convert_to_convexes, plot_convexes
    import matplotlib.pyplot as plt
    import numpy as np
    import tqdm
    # size of one detector
    each_width_mm = 19.52
    # the size of one detector + gap. 
    width_mm = 22.4
    EFL_mm = 4370.0
    l_center = 0.7
    b_center = 0.6
    # size of the Galactic centre survey
    dl_gcs = 0.7+1.4
    db_gcs = 1.2
    PA_deg = 0.0
#    Ndither = [26, 17]
# from the output of "number of dither l, b below.  
    Ndither = [264,174]
    gap_width_mm=width_mm - each_width_mm
    # make gap_width_mm smaller
    # gap_width_mm *= 0.8
    print(' gap width (mm)=', gap_width_mm)
    # deg vs mm
    mm_deg = EFL_mm*np.pi/180.0
    deg_mm = 1.0/mm_deg
    print(' 1 deg is ', mm_deg, ' mm')
    # dithering width (deg) -> (mm)
    dithering_width_mm = 0.01*mm_deg
    print('dithering width (mm) =', dithering_width_mm)

    # shift the starting point further by the size of one detector
    l_center += (each_width_mm+0.5*gap_width_mm-dithering_width_mm)*deg_mm
    b_center += (each_width_mm+0.5*gap_width_mm-dithering_width_mm)*deg_mm
    print(' starting central point (l, b) =', l_center, b_center)
    print(' size of one detector and gap (deg) =', each_width_mm*deg_mm, \
          gap_width_mm*deg_mm)
    add_dither_deg = (2.0*width_mm-gap_width_mm-dithering_width_mm)*deg_mm
    print('additional size to dither =', add_dither_deg)
    print('number of dither l, b=', \
         (add_dither_deg+dl_gcs)*mm_deg/dithering_width_mm, \
         (add_dither_deg+db_gcs)*mm_deg/dithering_width_mm)


    hdf = pkg_resources.resource_filename('telescope_baseline', 'data/cat_hw14.5.hdf')
    targets, l, b, hw = read_jasmine_targets(hdf)
    Nstar = 10**(hw/-2.5)/10**(12.5/-2.5)
    # square region    
    convexesset = ditheringmap(l_center, b_center, PA_deg, \
                     dithering_width_mm=dithering_width_mm, Ndither=Ndither, \
                     width_mm=width_mm, each_width_mm=each_width_mm, \
                     EFL_mm=EFL_mm, left=1.0, top=-0.75)

    print(' shape of the final convexesset =', np.shape(convexesset))

    pos = convert_to_convexes(convexesset)
    plot_convexes(l, b, pos)

    #check inout dithering map
    ans = inout_convexesset(targets,convexesset)

    #count number 
    nans = np.sum(ans,axis=0) #number of observation for each star after covering the whole fields  
    plot_n_targets(l, b, nans, cmap="CMRmap_r")    
    hist_n_targets(nans)

    # S/N computing, Kawata-san needs to correct the below.
    scale = 1.0
    # number of small frames per the Galactic centre field.
    nsf_per_gcf = np.shape(convexesset)[0]
    print(' number of small fields to cover the whole filed =', nsf_per_gcf)
    ac = 6000 # astrometric accuracy per exposure (micro arcsec per frame)
    norbits = 6000  # number of orbits used for the GC survey in 3 years
    nex_per_sf = 20 # number of exposure per small frame
    nsf_per_orb = 8 # number of small frames per orbit.
    norb_gcf = nsf_per_gcf/nsf_per_orb   # number of orbit required to cover the whole Gf survey area.
    print(' number of orbits to cover the whole GC filed =', norb_gcf)
    scale = nex_per_sf*(norbits/norb_gcf)
    final_ac=ac/np.sqrt(nans*scale)/np.sqrt(Nstar)    

    plot_ae_targets(l, b, final_ac, cmap="CMRmap_r")
    hist_ae_targets(final_ac)

    # save the data
    with h5py.File("MPSv3r3_rect.h5", "w") as f:
      dsetl = f.create_dataset("l", data=l)
      dsetb = f.create_dataset("b", data=b)
      dsethw = f.create_dataset("hw", data=hw)
      dsetnans = f.create_dataset("nans", data=nans)
      dsetfac = f.create_dataset("final_ac", data=final_ac)      


