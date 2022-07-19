if __name__ == "__main__":
    import pkg_resources           
    from telescope_baseline.mapping.read_catalog import read_jasmine_targets
    from telescope_baseline.mapping.plot_mapping import plot_targets, plot_n_targets, hist_n_targets, plot_ae_targets, hist_ae_targets
    from telescope_baseline.mapping.mapset import obsn_MPSv1
    import matplotlib.pyplot as plt
    import numpy as np
    import tqdm
    each_width_mm=19.52
    width_mm=22.4
    EFL_mm=4370.0
    l_center=-1.2
    b_center=0.0
    PA_deg=0.0

    hdf=pkg_resources.resource_filename('telescope_baseline', 'data/cat_hw14.5.hdf')
    targets,l,b,hw=read_jasmine_targets(hdf)
    Nstar=10**(hw/-2.5)/10**(12.5/-2.5)
    
    Ng=11
    dizLmax=1.5
    grid=np.linspace(-dizLmax,dizLmax,Ng)
    dummy_array=np.ones((Ng,Ng))
    gx=(grid[:,np.newaxis]*dummy_array).flatten()
    gy=(grid[np.newaxis,:]*dummy_array).flatten()
    
    nans=np.zeros_like(hw)
    for i in tqdm.tqdm(range(0,Ng*Ng)):
        nans_each,pos_each=obsn_MPSv1(targets,l_center,b_center,PA_deg, width_mm=width_mm, each_width_mm=each_width_mm, EFL_mm=EFL_mm,left=gx[i],top=gy[i])
        nans=nans+nans_each
        if i==0:
            pos_all=pos_each
        else:
            pos_all=np.vstack([pos_all,pos_each])

    scale=50.0*6000/12.0/(Ng*Ng)
    
    print("Nobs=",len(nans[nans>0]))
    ac=6000 # micro arcsec per frame
    
    final_ac=ac/np.sqrt(nans*scale)/np.sqrt(Nstar)    

    plot_ae_targets(l,b,final_ac,cmap="CMRmap_r")    
    plot_n_targets(l,b,scale*nans,pos=pos_all,cmap="CMRmap_r")    
    plot_n_targets(l,b,scale*nans,cmap="CMRmap_r",outfile="nno.png")    
    hist_ae_targets(final_ac)
    hist_n_targets(scale*nans)
