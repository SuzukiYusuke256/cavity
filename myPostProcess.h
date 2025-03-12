#ifndef _MY_POST_PROCESS_H_
#define _MY_POST_PROCESS_H_

#include "myIO.h"

/**
 * Function to interpolate velocities from staggered grid to cell centers
 * 
 * @param u          x-direction velocity field (on staggered grid)
 * @param v          y-direction velocity field (on staggered grid)
 * @param uCenter    array to store x-direction velocity at cell centers
 * @param vCenter    array to store y-direction velocity at cell centers
 * @param xn         number of cells in x-direction (internal domain)
 * @param yn         number of cells in y-direction (internal domain)
 * @return           0 on success, -1 on failure
 */
int calcCellCenterVelocity(double* u, double* v, 
                           double* uCenter, double* vCenter,
                           int xn, int yn);

/**
 * Function to calculate velocity gradients at cell centers
 * 
 * @param uCenter    x-direction velocity at cell centers
 * @param vCenter    y-direction velocity at cell centers
 * @param dudx       array to store ∂u/∂x at cell centers
 * @param dudy       array to store ∂u/∂y at cell centers
 * @param dvdx       array to store ∂v/∂x at cell centers
 * @param dvdy       array to store ∂v/∂y at cell centers
 * @param xn         number of cells in x-direction (internal domain)
 * @param yn         number of cells in y-direction (internal domain)
 * @param rdx        inverse of grid spacing in x-direction (1/dx)
 * @param rdy        inverse of grid spacing in y-direction (1/dy)
 * @return           0 on success, -1 on failure
 */
int calcCellCenterVelocityGradients(double* uCenter, double* vCenter, 
                                   double* dudx, double* dudy, 
                                   double* dvdx, double* dvdy,
                                   int xn, int yn, 
                                   double rdx, double rdy);

int calcSurfaceVelocityGradients(double* uCenter, double* vCenter, double* wallDudy, int xn, int yn, double rdx, double rdy);

/**
 * Integrated function to calculate velocity gradients at cell centers from staggered grid velocities
 * 
 * @param u          x-direction velocity field (on staggered grid)
 * @param v          y-direction velocity field (on staggered grid)
 * @param dudx       array to store ∂u/∂x at cell centers
 * @param dudy       array to store ∂u/∂y at cell centers
 * @param dvdx       array to store ∂v/∂x at cell centers
 * @param dvdy       array to store ∂v/∂y at cell centers
 * @param xn         number of cells in x-direction (internal domain)
 * @param yn         number of cells in y-direction (internal domain)
 * @param rdx        inverse of grid spacing in x-direction (1/dx)
 * @param rdy        inverse of grid spacing in y-direction (1/dy)
 * @return           0 on success, -1 on failure
 */
int calcStaggeredToCellCenterVelocityGradients(double* u, double* v, 
                                              double* dudx, double* dudy, 
                                              double* dvdx, double* dvdy,
                                              int xn, int yn, double rdx, double rdy);

#endif