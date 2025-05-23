A variant of FORTRAN code QSSP, for calculating co-seismic static deformation of a self-gravitating, spherically symmetric, isotropic and elastic earth.

This code package consists of two program codes: The first program "QSSPgrnco" generates co-seismic Green's function database, and the second program "Getcoseis" provides Green's functions or complete co-seismic deformation field requested by the user.

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../QSSPgrnco2018/SourceCode

~>gfortran -o qsspgrnco *.f -O3

to get the excutable code qsspgrnco.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

References

Wang, R., S. Heimann, Y. Zhang, H. Wang, and T. Dahm (2017). Complete synthetic seismograms based on a spherical self-gravitating Earth model with an atmos-phere-ocean-mantle-core structure. Geophysical Journal International, doi: 10.1093/gji/ggx259.
