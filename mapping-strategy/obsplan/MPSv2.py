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
    each_width_mm=19.52
    width_mm=22.4
    EFL_mm=4370.0
    l_center=0.7
    b_center=0.7
    PA_deg=0.0
    Ndither=[37,37]
    gap_width_mm=width_mm - each_width_mm
    # make gap_width_mm smaller
    # gap_width_mm *= 0.8
    print(' gap width (mm)=', gap_width_mm)

    hdf=pkg_resources.resource_filename('telescope_baseline', 'data/cat_hw14.5.hdf')
    targets,l,b,hw=read_jasmine_targets(hdf)
    Nstar=10**(hw/-2.5)/10**(12.5/-2.5)
    # square region    
    convexesset_sq=ditheringmap(l_center,b_center,PA_deg, dithering_width_mm=gap_width_mm, Ndither=Ndither,
                             width_mm=width_mm, each_width_mm=each_width_mm, EFL_mm=EFL_mm, left=1.0,top=-0.75)

    # rectoangular region
    l_center_ra=-0.7
    b_center_ra=0.3
    PA_deg=0.0
    Ndither_ra=[21,16]
    convexesset_ra=ditheringmap(l_center_ra, b_center_ra,PA_deg, dithering_width_mm=gap_width_mm, Ndither=Ndither_ra,
                             width_mm=width_mm, each_width_mm=each_width_mm, EFL_mm=EFL_mm, left=1.0,top=-0.75)

    print(' shape of convexesset square=', np.shape(convexesset_sq))
    print(' shape of convexesset rectoangular=', np.shape(convexesset_ra))    
    convexesset = np.concatenate([convexesset_sq, convexesset_ra])
    print(' shape of the final convexesset =', np.shape(convexesset))

    pos=convert_to_convexes(convexesset)
    plot_convexes(l,b,pos)

    #check inout dithering map
    ans=inout_convexesset(targets,convexesset)

    #count number 
    nans=np.sum(ans,axis=0) #number of obs    
    plot_n_targets(l,b,nans,cmap="CMRmap_r")    
    hist_n_targets(nans)


    # S/N computing, Kawata-san needs to correct the below.
    scale=1.0
    # number of small frames per the Galactic centre field.
    nsf_per_gcf = np.shape(convexesset)[0]
    print(' number of small fields to cover the whole filed =', nsf_per_gcf)
    ac=6000 # astrometric accuracy per exposure (micro arcsec per frame)
    norbits = 6000  # number of orbits used for the GC survey in 3 years
    nex_per_sf = 20 # number of exposure per small frame
    nsf_per_orb = 8 # number of small frames per orbit.
    norb_gcf = nsf_per_gcf/nsf_per_orb   # number of orbit required to cover the whole Gf survey area.
    print(' number of orbits to cover the whole GC filed =', norb_gcf)
    scale = nex_per_sf*(norbits/norb_gcf)
    final_ac=ac/np.sqrt(nans*scale)/np.sqrt(Nstar)    

    plot_ae_targets(l,b,final_ac,cmap="CMRmap_r")
    hist_ae_targets(final_ac)

    # save the data
    with h5py.File("MPSv2.h5", "w") as f:
      dsetl = f.create_dataset("l", data=l)
      dsetb = f.create_dataset("b", data=b)
      dsethw = f.create_dataset("hw", data=hw)
      dsetnans = f.create_dataset("nans", data=nans)
      dsetfac = f.create_dataset("final_ac", data=final_ac)      


