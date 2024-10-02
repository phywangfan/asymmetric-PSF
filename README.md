# asymmetric-PSF
This repository is a MATLAB-based open source code corresponding to the manuscript tilted by "Asymmetric point-spread function in the tilted plane". 

There are two main programs for implementing designed points and a 3D object, named by 'main_aPSF_points.m' and 'main_aPS_pointCGHF.m', respectively.

The program 'main_aPSF_points.m' is performed to generate Fig. 2 and Fig. 4 by changing parameters in the code.
The program 'main_aPS_pointCGHF.m' is performed to generate Fig. 6.

# object
'chess_two_2.ply',  is two pieces of pawns cosisting of point cloud.

# Sub-functions
1. 'intersection_LinePlane.m', is to solve for the four intersections of the four radial lines emanating from the source point with the tilted plane.
2. 'rotmatrix.m', is to obtain the rotation matrix by inputing the given rotated angle. Noting that this is not the matrix for calculating aPSF, but for obtaining the normal vector of tilted plane.
3. 'RotMatrix_Rodrigues.m', is to obtain the rotation matrix by inputing the normal vector, which is used for calculating aPSF.
4. 'rot_poly2D.m', is to rotate the wavefield from the global system to the local system by assigning a value of ‘true’ to the last variable, while rotate from the local system to the global system by assigning a value of "false" to the last variable.
5. 'ASM_prop_W.m', is to propagate wavefield forward based on angular spectrum. 


