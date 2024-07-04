import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pylab as plt
import astropy.units as u
from gal_radii_pb import correct_rgc_mc

def project_tp(centre, coord):
    """
    Function to project RA and dec onto the tangent plane given a specific centre.

    Inputs: 
    -------
    centre - Astropy coordinate for your object 
    coord  - Astropy coordinate for your object
   

    Returns:
    --------
    xi, eta - The angular coordinates projected onto a tangent plane centred on centre
    """
    sep = centre.separation(coord)
    eta = np.cos(coord.dec.rad)*np.sin(coord.ra.rad-centre.ra.rad)
    xi = (np.cos(centre.dec.rad)*np.sin(coord.dec.rad)-np.sin(centre.dec.rad)*np.cos(coord.dec.rad)
        * np.cos(coord.ra.rad-centre.ra.rad))/np.cos(sep.rad)
    return xi,eta

def aitoff_error_conn(coord, dist, d_eu, d_ed):
    """
    Function that samples along the line of sight distance uncertainty
    to create an uncertainty track for an aitoff plot. Conn projection.
    
    Inputs: 
    -------
    coord - Astropy coordinate for your object
    dist  - Distance to object in kpc
    d_eu, d_ed - upper and lower distance uncertainties in kpc.

    Returns:
    --------
    l_track, b_track - numpy arrays with positions of uncertainty in l and b
    """
    dist_samp = np.linspace(dist-d_eu, dist+d_ed,100)
    l_track = np.zeros(100)
    b_track = np.zeros(100)
    for i in range(100):
        l_track[i], b_track[i], r_tmp = radec2m31_conn(coord,dist_samp[i])
 
    return l_track, b_track


def aitoff_error_metz(coord, dist, d_eu, d_ed):
    """
    Function that samples along the line of sight distance uncertainty
    to create an uncertainty track for an aitoff plot. Metz projection.
    Inputs: 
    -------
    coord - Astropy coordinate for your object
    dist  - Distance to object in kpc
    d_eu, d_ed - upper and lower distance uncertainties in kpc.

    Returns:
    --------
    l_track, b_track - numpy arrays with positions of uncertainty in l and b
    """
    dist_samp = np.linspace(dist-d_eu, dist+d_ed,100)
    l_track = np.zeros(100)
    b_track = np.zeros(100)
    for i in range(100):
        l_track[i], b_track[i], r_tmp = radec2m31_metz(coord,dist_samp[i])
 
    return l_track, b_track


# Define rotation matrices 

def rotate_x(angle):
    """
    Rotate points around the x-axis by the given angle in radians.
    NB Metz has the -ve signs flipped in their example rotation.

    Inputs: 
    -------
    angle - angle to rotate around, defined in radians

    Returns:
    --------
    rotation_matrix - calculated rotation matrix

    """
    rotation_matrix = np.array([[1, 0, 0],
                                [0, np.cos(angle), -np.sin(angle)],
                                [0, np.sin(angle), np.cos(angle)]])
    return rotation_matrix

def rotate_y(angle):
    """
    Rotate points around the y-axis by the given angle in radians.

    Inputs: 
    -------
    angle - angle to rotate around, defined in radians

    Returns:
    --------
    rotation_matrix - calculated rotation matrix

    """
    rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle)],
                                [0, 1, 0],
                                [-np.sin(angle), 0, np.cos(angle)]])
    return rotation_matrix

def rotate_z(angle):
    """
    Rotate points around the z-axis by the given angle in radians.

    Inputs: 
    -------
    angle - angle to rotate around, defined in radians

    Returns:
    --------
    rotation_matrix - calculated rotation matrix

    """
    rotation_matrix = np.array([[np.cos(angle), np.sin(angle), 0],
                                [-np.sin(angle), np.cos(angle), 0],
                                [0, 0, 1]])
    return rotation_matrix

def rotate_rpq(ra,dec):

	Rrpq = np.array([[np.cos(dec)*np.cos(ra),-np.sin(ra),-np.sin(dec)*np.cos(ra)],
					[np.cos(dec)*np.sin(ra),np.cos(ra),-np.sin(dec)*np.sin(ra)],
					[np.sin(dec),0.,np.cos(dec)]])
	return Rrpq




def radec2m31_metz(coord,dist): # heliocentric not galactocentric 
    """
    Convert equatorial coordinates (RA, Dec) to Andromeda-centric coordinates.
    Following Metz et al 2007 
    (https://ui.adsabs.harvard.edu/abs/2007MNRAS.374.1125M/abstract)

    Important note: I have followed Metz et al but added in:
    (1) An extra rotation around z of np.pi to match projections in other papers 
    (Savino et al. 2022, Conn et al. 2013, probably due to AItoff axes in matplotlib)
    (2) This is flipped/reflected in orientation along the l axis cf. Conn/Savino. 

    Also, as noted by myself and Paula, there is an error in the Metz appendix with where their
    matrices are transposed. The below seems to output the Metz results (particularly their
    values for R_M31) correctly. 

    Inputs:
    --------

    astropy coordinates object to transform
    distance to objects
    
    Returns:
    --------

    l_m31, b_m31, r_new: M31 galactic coordinates (degrees, degrees, kpc)

    """

    # Define M31 frame (coordinates, inclination and position angle)
    theta = np.pi/2-np.deg2rad(39.5)
    inc1 = np.pi/2-np.deg2rad(77.5)
    # angle beta=-0.6219604 corrects the offset between the Sun and MW Galactic Center
    # i.e. rotation of the x and y planes about the z axis in a counter-clockwise direction, for a left-handed coordinate system
    beta = 0.6219604*(np.pi/180.)
    c_m31 = SkyCoord(10.68470833, 41.26875, unit = (u.deg,u.deg))
    r_m31 = 785


    x = dist * np.cos(coord.ra.rad) * np.cos(coord.dec.rad)
    y = dist * np.sin(coord.ra.rad) * np.cos(coord.dec.rad)
    z = dist * np.sin(coord.dec.rad)
    cart = np.array([x,y,z])
    xM31 = r_m31 * np.cos(c_m31.ra.rad) * np.cos(c_m31.dec.rad)
    yM31 = r_m31 * np.sin(c_m31.ra.rad) * np.cos(c_m31.dec.rad)
    zM31 = r_m31 * np.sin(c_m31.dec.rad)
    M31cart = np.array([xM31,yM31,zM31])

    # Added an additional rotation around the z axis of 180 degrees to line up
    #  with Savino/Conn/matplotlib projection (albeit reflected)
    result = np.matmul(rotate_z(np.pi).T,np.matmul(rotate_z(beta).T, np.matmul(rotate_y(inc1).T, 
                       np.matmul(rotate_x(theta).T, 
                       rotate_rpq(c_m31.ra.rad,c_m31.dec.rad).T))))

    
    R_M31 = np.matmul(result,(cart-M31cart)) # Matrix for getting from MW-centric to M31-centric as: x_M31 = R_M31*x_MW
    r_new = np.sqrt(np.sum(R_M31**2))
 

    l_m31 = np.rad2deg(np.arctan2(R_M31[1],R_M31[0]))
    # Correct for l<0 
    if l_m31 < 0:
        l_m31+=360
    
    b_m31 = np.rad2deg(np.arcsin(R_M31[2]/r_new))

    return l_m31, b_m31, r_new

def radec2m31_conn(coord,dist):
    """
    Convert equatorial coordinates (RA, Dec) to Andromeda-centric coordinates.
    Following Conn et al 2013 (and used in Savino et al. 2022) 
    (https://ui.adsabs.harvard.edu/abs/2013ApJ...766..120C/abstract
    https://ui.adsabs.harvard.edu/abs/2012ApJ...758...11C/abstract
    https://ui.adsabs.harvard.edu/abs/2022ApJ...938..101S/abstract)
    
    Inputs:
    -------
    coord - astropy coordinates object to transform 
    dist - distance to objects in kpc
    Returns:
    --------
    l_m31, b_m31: M31 galactic coordinates (degrees, degrees, kpc)

    """
    # Define M31 frame (coordinates, inclination and position angle)
    theta = np.deg2rad(39.5)
    inc1 = -np.deg2rad(77.5)
    c_m31 = SkyCoord(10.68470833, 41.26875, unit = (u.deg,u.deg))
    r_m31 = 785

    xi, eta = project_tp(c_m31,coord)

    # X, Y, Z in M31 coordinates
    ang_d = coord.separation(c_m31)    
    x = dist * np.cos(ang_d.rad)*np.tan(xi)    
    y = dist * np.sin(eta)
    z = dist * np.cos(ang_d.rad)-785
    M31cart = np.array((x,y,z))

    # Conn 2013 rotation

    R_M31 = np.matmul(rotate_z(np.pi/2),np.matmul(rotate_x(inc1),np.matmul(rotate_z(theta),M31cart)))

    r_new = np.sqrt(np.sum(R_M31**2))
 
    l_m31 = np.rad2deg(np.arctan2(R_M31[1],R_M31[0])) +180
    # Correct for l<0 
    if l_m31 < 0:
        l_m31+=360
    
    b_m31 = np.rad2deg(np.arcsin(R_M31[2]/r_new))

    return l_m31, b_m31, r_new



