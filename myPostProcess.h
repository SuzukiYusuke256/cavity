#ifndef _MY_POST_PROCESS_H_
#define _MY_POST_PROCESS_H_

#include "myIO.h"

// Function to interpolate velocities from staggered grid to cell centers
int calcCellCenterVelocity(double* u, double* v, 
                           double* uCenter, double* vCenter,
                           int xn, int yn);

// Function to calculate velocity gradients at cell centers
int calcCellCenterVelocityGradients(double* uCenter, double* vCenter, 
                                   double* dudx, double* dudy, 
                                   double* dvdx, double* dvdy,
                                   int xn, int yn, 
                                   double rdx, double rdy);

int calcSurfaceVelocityGradients(double* uCenter, double* vCenter, double* wallDudy, int xn, int yn, double rdx, double rdy);

// Integrated function to calculate velocity gradients at cell centers 
// from staggered grid velocities
int calcStaggeredToCellCenterVelocityGradients(double* u, double* v, 
                                              double* dudx, double* dudy, 
                                              double* dvdx, double* dvdy,
                                              int xn, int yn, double rdx, double rdy);

// Function to calculate viscous dissipation at cell centers.
int calcViscousDissipation(double nu, double* dudx, double* dudy, double* dvdx, double* dvdy, int xn, int yn, double dx, double dy, double* dissip, double* totalDissip);

// Function to calculate wall work
int calcWallWork(double nu, double* wallDudy, int xn, int yn, double dx, double dy, double* wallWork, double* totalWallWork);

// Function to calculate kinetic energy
int calcKineticEnergy(double* u, double* v, int xn, int yn, double dx, double dy, double* ke, double* totalKe);

void print_usage(const char* programName);

#endif