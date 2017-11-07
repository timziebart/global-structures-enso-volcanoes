# Introduction

This repository was created to publish the code for the analysis in [1]. It is meant for recreation of the results in that work and a starting point for further functional climate network analysis.

*Disclaimer: As the main goal of this code publication is the reproduction of the results in [1], the code is not well-documented. If you wish to use this work as a basis for your own work, please contact me (Tim.Kittel@pik-potsdam.de) .*

# Prerequesites

# The Repository

Please clone the repository and then go inside (for the following commands to work).
```
git clone https://github.com/timkittel/global-structures-enso-volcanoes.git
cd global-structures-enso-volcanoes
```

## Data

As detailed in [1] this project uses the **daily averaged surface air temperature data (SAT at sigma=0.995)** of the [NCEP/NCAR Reanalysis I Project](https://www.esrl.noaa.gov/psd/data/reanalysis/reanalysis.shtml). (They provide [ftp download servers](ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/) and you basically need all the files `air.sig995.1948.nc` until `air.sig995.2015.nc`.)

The easiest is if create a subfolder called `data` and place all the `.nc` files inside. If you prefer to keep them somewhere else, you can use the `--data-directory` flag of `fullrun.py`

## Python & Dependencies

The Code was developed on an `64bit Ubuntu 16.04 LTS` using `Python 3.6.1`. Furthermore, the following additional Python packages are necessary. The version numbers correspond to the ones I used during the development.

* `numpy 1.12.1`
* `scipy 0.19.0`
* `matplotlib 2.0.2`
* `basemap 1.1.0`
* `pandas 0.20.1`
* `argcomplete 1.8.2`
* `igraph 0.7.1.post6`
* `h5py 2.7.0`
* `netCDF4 1.2.8`

Further, please install `simple-mpi` which can be found [here on github](https://github.com/timkittel/simple-mpi). (I recommend to clone it in a different folder, install it with `pip` there and then return to the `global-structures-enso-volcanoes` folder.)

# Running the Code

*Comment: Running the below analysis is computationally very expensive. If the following code detects `mpi` (e.g. because it was called via `mpirun`) it automatically parallelizes the computation. I strongly recommend to use this option! I am happy to provide details/help if necessary.*

*Comment: `./fullrun.py -h` gives a help for the usage. Please check it before running the lines below.*

*Comment: `./fullrun.py` creates (possibly many) temporary files when using `mpi`. If you want these to be save in a different location, use the `--scratch-directory` flag.*

First, create the cache file for the icosahedral grid `.icosahedral-grid.cache.npy`. *In this step, do not use `mpi`!*
```
./icosahedral_grid.py
```
Then run the functional climate network analysis (except modularity). A file called `Output.FullRun.daily-paper.icosahedral.hdf5` will be created. *I strongly recommend the usage of `mpi`, e.g. using `mpirun`.!*
```
./fullrun.py daily paper
```
The computation of the modularity is not done daily, because the computation time would increase drastically. Hence, it is done in a separate step. A file called `Output.FullRun.normal-modularity.icosahedral.hdf5` will be created.
```
./fullrun.py normal modularity
```
The analysis for the volcanoes is done with the icosahedral grid where the ENSO-big region is removed, using the following line. A file called `Output.FullRun.daily-paper.icosahedral_without_ENSO_big.hdf5` will be created.
```
./fullrun.py daily paper --grid icosahedral_without_ENSO_big
```
Finally, if you want to recreate the comparison of the community detection algorithms, run the following line. A file called `Output.FullRun.normal-cmp-modularity.icosahedral.hdf5` will be created.
```
./fullrun.py normal comparison-modularity
```


# Plotting the Results

*Comment: `./paper-pix.py -h` gives a help for the usage. Please check it before running the lines below.*

regions of interest (Fig. 1 in [1]):
```
./paper-pix.py Output.FullRun.daily-paper.icosahedral.hdf5 icosahedral regions-of-interest --save
```

transitivity and global average link distance (Fig. 2(a) and (c) in [1]):
```
./paper-pix.py Output.FullRun.daily-paper.icosahedral.hdf5 icosahedral ENSO-global --save
```

modularity (Fig. 2(b) in [1]):
```
/paper-pix.py Output.FullRun.normal-modularity.icosahedral.hdf5 icosahedral ENSO-global --save
```

composites (Fig. 3 in [1]):
```
./paper-pix.py Output.FullRun.daily-paper.icosahedral.hdf5 icosahedral composites --save
```

regionalized climate network properties (Fig. 4 in [1]):
```
./paper-pix.py Output.FullRun.daily-paper.icosahedral.hdf5 icosahedral ENSO-local --save
```

analysis of volcanoes (Fig. 5 in [1]):
```
./paper-pix.py Output.FullRun.daily-paper.icosahedral_without_ENSO_big.hdf5 icosahedral_without_ENSO_big volcanoes --save
```

comparison of the community detection algorithms
```
./paper-pix.py Output.FullRun.normal-cmp-modularity.icosahedral.hdf5 icosahedral cmp-modularity
```

# References

[1] T. Kittel, C. Ciemer, N. Lotfi, T. Peron, F. Rodrigues,
J. Kurths, and R. V. Donner. “Global teleconnectivity structures of the
El Niño–Southern Oscillation and large volcanic eruptions: An evolving
network perspective”. in prep. 2017















