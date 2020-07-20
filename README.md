
## IsoFlow 0.5.3
- Fixed issue with inverse metric indexing in residual calculation
- Added plotting of wall temperature and pressures
- Tweaked capsule meshing to provide less skewed meshes

## IsoFlow 0.5.2
- Separated copy function of object pressure coefficient

## IsoFlow 0.5.1
- Residual plotting no longer resets on every simulation run, only resets on initialize

###########################################################################
## IsoFlow 0.5.0
###########################################################################

- Added Roe FVS scheme
- Added pressure coefficient plotting
- Added ability to copy pressure coefficient data to system clipboard
- Implemented Fortran subroutine for calculating interface sound speeds

###########################################################################
## IsoFlow 0.4.0
###########################################################################

- Added MUSCL interpolation
- Higher order scheme menu and limiters added
- Overhauled boundary condition handling, added dynamic outflow switching for subsonic vs. supersonic outflow

###########################################################################
## IsoFlow 0.3.0
###########################################################################

-Added capability for NACA and biconvex airfoils
-Updated boundary conditions for use with airfoil geometries
-Added ability to change angle of attack of inlet flow

-Added SLAU solver 