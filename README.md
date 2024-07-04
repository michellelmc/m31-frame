# m31-frame
Python code for converting RA and dec to M31-centric galactic coordinates. Written in collaboration with Paula Gherghinescu.

# Description

m31-frame allows you to convert coordinates in RA and dec to an M31-centric galactic coordinate system ($l_{M31}, b_{M31}$). In m31_frame.py, we provide functions that will produce these coordinates using transformations from either [Metz et al 2007](https://ui.adsabs.harvard.edu/abs/2007MNRAS.374.1125M/abstract) (radec2m31_metz) or [Conn et al. 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...758...11C/abstract) and [2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...766..120C/abstract) (radec2m31_conn), also used in [Savino et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...938..101S/abstract).

Functionality is also provide for calculating tangent plane coordinates for M31 (project_tp). 

Finally, we provide functions for calculating uncertainty tracks to use in an Aitoff plot (aitoff_error_metz, aitoff_error_conn).

# Documentation

We provide worked examples in m31_frame.ipynb.
