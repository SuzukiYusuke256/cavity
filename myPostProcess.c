#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "myConst.h"
#include "config.h"
#include "myIO.h"
#include "myPostProcess.h"

// calculate gradient
int main()
{
    // const char* caseName = "cavity_smac_01-2";
    // int timeStep         = 0;

    const char* caseName = "test_02";
    int timeStep         = 0;
    
    // read config
    Config cfg;
    char cfgName[128] = "";
    strcpy(cfgName,caseName); // cfgName = caseName
    strcat(cfgName,"/config"); // cfgName = caseName/config

    readConfig(cfgName,&cfg);

    // if the case name does not match the config file, return error
    if (strcmp(caseName, cfg.caseName) != 0) 
    {
        fprintf(stderr, "Error: caseName does not match cfg.caseName\n");
        return -1;
    }

    const int nx = cfg.nx;
    const int ny = cfg.ny;

    const double rdx = (double)nx; // 0.25 dx = 1/xn. rdx = 1/(1/xn) = xn;
    const double rdy = (double)ny;

    double* u          = (double*)calloc((nx+3)*(ny+2), sizeof(double));
    double* v          = (double*)calloc((nx+2)*(ny+3), sizeof(double));
    double* uCenter    = (double*)calloc((nx+2)*(ny+2), sizeof(double));
    double* vCenter    = (double*)calloc((nx+2)*(ny+2), sizeof(double));
    double* dudx       = (double*)calloc(nx*ny, sizeof(double));
    double* dudy       = (double*)calloc(nx*ny, sizeof(double));
    double* dvdx       = (double*)calloc(nx*ny, sizeof(double));
    double* dvdy       = (double*)calloc(nx*ny, sizeof(double));
    double* wallDudy   = (double*)calloc(nx, sizeof(double));

    double* dissip     = (double*)calloc(nx*ny, sizeof(double));
    double totalDissip = 0.0;

    readData(u,nx+3,ny+2,caseName,timeStep,"U");
    readData(v,nx+2,ny+3,caseName,timeStep,"V");

    // Calculate velocities at cell centers
    calcCellCenterVelocity(u, v, uCenter, vCenter, nx, ny);

    // Calculate velocity gradients at cell centers
    calcCellCenterVelocityGradients(uCenter, vCenter, dudx, dudy, dvdx, dvdy, nx, ny, rdx, rdy);
    calcSurfaceVelocityGradients(uCenter,vCenter,wallDudy,nx,ny,rdx,rdy);
    calcViscousDissipation(0.01, dudx, dudy, dvdx, dvdy, dissip, &totalDissip, nx, ny);

    writeData(uCenter, nx+2,ny+2,caseName,timeStep,"vU",           nx,ny); // U at the cell center
    writeData(vCenter, nx+2,ny+2,caseName,timeStep,"vV",           nx,ny); // V at the cell center
    writeData(dudx,    nx,  ny,  caseName,timeStep,"vDUdx",        nx,ny);
    writeData(dudy,    nx,  ny,  caseName,timeStep,"vDUdy",        nx,ny);
    writeData(dvdx,    nx,  ny,  caseName,timeStep,"vDVdx",        nx,ny);
    writeData(dvdy,    nx,  ny,  caseName,timeStep,"vDVdy",        nx,ny);
    writeData(wallDudy,nx,  1,   caseName,timeStep,"sDUdy",        nx,ny);
    writeData(dissip,  nx,  ny,  caseName,timeStep,"viscousDissip",nx,ny);

    printf("Total viscous dissipation: %lf\n", totalDissip);
    // writeData(totalDissip, 1,1,  caseName,0,"totalViscousDissip",nx,ny);
    
    
    // Free memory
    free(dudx);
    free(dudy);
    free(dvdx);
    free(dvdy);
    free(u);
    free(v);
    free(uCenter);
    free(vCenter);
    free(wallDudy);
    free(dissip);
    
    return 0;
}

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
int calcCellCenterVelocity(double* u, double* v, double* uCenter, double* vCenter, int xn, int yn) 
{
    // Error check
    if (u == NULL || v == NULL || uCenter == NULL || vCenter == NULL) {
        fprintf(stderr, "Error: Null pointer provided to calcCellCenterVelocity\n");
        return -1;
    }
    
    // Interpolate velocities from staggered grid to cell centers
    for (int jj = 0; jj < yn+2; jj++) {
        for (int ii = 0; ii < xn+2; ii++) {
            // u is offset in x-direction (located at (i+0.5,j)), average adjacent u values
            uCenter[ii + jj*(xn+2)] = 0.5 * (u[ii + jj*(xn+3)] + u[ii+1 + jj*(xn+3)]);
            // v is offset in y-direction (located at (i,j+0.5)), average adjacent v values
            vCenter[ii + jj*(xn+2)] = 0.5 * (v[ii + jj*(xn+2)] + v[ii + (jj+1)*(xn+2)]);

            // printf("%d %d %lf\n",ii,jj,u[ii + jj*(xn+3)]);
            // printf("%d %d %lf\n",ii,jj,uCenter[ii + jj*xn]);
        }
    }
    
    return 0;
}

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
                                   double rdx, double rdy) {
    
    if (uCenter == NULL || vCenter == NULL || dudx == NULL || dudy == NULL || dvdx == NULL || dvdy == NULL) {
        fprintf(stderr, "Error: Null pointer provided to calcCellCenterVelocityGradients\n");
        return -1;
    }
    
    // Calculate gradients for internal cells
    // the size of array is (xn x yn) (ex) 128x128
    for (int jj = 1; jj <= yn; jj++) {
        for (int ii = 1; ii <= xn; ii++) {

            // Cell index (internal domain)
            // the size of uCenter is (xn+2 x yn+2) including ghost cells
            int idx = (ii-1) + (jj-1)*xn;
            
            // Calculate gradients using central difference method
            dudx[idx] = (uCenter[(ii+1) + jj*(xn+2)] - uCenter[(ii-1) + jj*(xn+2)]) * 0.5 * rdx;
            dudy[idx] = (uCenter[ii + (jj+1)*(xn+2)] - uCenter[ii + (jj-1)*(xn+2)]) * 0.5 * rdy;
            dvdx[idx] = (vCenter[(ii+1) + jj*(xn+2)] - vCenter[(ii-1) + jj*(xn+2)]) * 0.5 * rdx;
            dvdy[idx] = (vCenter[ii + (jj+1)*(xn+2)] - vCenter[ii + (jj-1)*(xn+2)]) * 0.5 * rdy;
        }
    }
    
    // Calculate gradients for boundary cells using one-sided differences
    // Left and right boundaries
    // for (int jj = 1; jj < yn-1; jj++) {
    //     // Left boundary (i=0)
    //     dudx[0 + jj*xn] = (uCenter[1 + jj*xn] - uCenter[0 + jj*xn]) * rdx;
    //     dudy[0 + jj*xn] = (uCenter[0 + (jj+1)*xn] - uCenter[0 + (j-1)*xn]) * r2dy;
    //     dvdx[0 + jj*xn] = (vCenter[1 + jj*xn] - vCenter[0 + jj*xn]) * rdx;
    //     dvdy[0 + jj*xn] = (vCenter[0 + (j+1)*xn] - vCenter[0 + (j-1)*xn]) * r2dy;
        
    //     // Right boundary (i=xn-1)
    //     dudx[(xn-1) + jj*xn] = (uCenter[(xn-1) + jj*xn] - uCenter[(xn-2) + jj*xn]) * rdx;
    //     dudy[(xn-1) + jj*xn] = (uCenter[(xn-1) + (j+1)*xn] - uCenter[(xn-1) + (j-1)*xn]) * r2dy;
    //     dvdx[(xn-1) + jj*xn] = (vCenter[(xn-1) + jj*xn] - vCenter[(xn-2) + jj*xn]) * rdx;
    //     dvdy[(xn-1) + jj*xn] = (vCenter[(xn-1) + (j+1)*xn] - vCenter[(xn-1) + (j-1)*xn]) * r2dy;
    // }
    
    // // Top and bottom boundaries
    // for (int i = 1; i < xn-1; i++) {
    //     // Bottom boundary (j=0)
    //     dudx[i + 0*xn] = (uCenter[(i+1) + 0*xn] - uCenter[(i-1) + 0*xn]) * r2dx;
    //     dudy[i + 0*xn] = (uCenter[i + 1*xn] - uCenter[i + 0*xn]) * rdy;
    //     dvdx[i + 0*xn] = (vCenter[(i+1) + 0*xn] - vCenter[(i-1) + 0*xn]) * r2dx;
    //     dvdy[i + 0*xn] = (vCenter[i + 1*xn] - vCenter[i + 0*xn]) * rdy;
        
    //     // Top boundary (j=yn-1)
    //     dudx[i + (yn-1)*xn] = (uCenter[(i+1) + (yn-1)*xn] - uCenter[(i-1) + (yn-1)*xn]) * r2dx;
    //     dudy[i + (yn-1)*xn] = (uCenter[i + (yn-1)*xn] - uCenter[i + (yn-2)*xn]) * rdy;
    //     dvdx[i + (yn-1)*xn] = (vCenter[(i+1) + (yn-1)*xn] - vCenter[(i-1) + (yn-1)*xn]) * r2dx;
    //     dvdy[i + (yn-1)*xn] = (vCenter[i + (yn-1)*xn] - vCenter[i + (yn-2)*xn]) * rdy;
    // }
    
    return 0;
}

// velocity gradient at the moving wall. only x-direction
int calcSurfaceVelocityGradients(double* uCenter, double* vCenter, double* wallDudy, int xn, int yn, double rdx, double rdy) 
{

    // Calculate gradients for internal cells
    // the size of array is (xn x yn) (ex) 128x128
    for (int ii = 0; ii < xn; ii++) {

        // Calculate gradients using central difference method with inverse multiplication
        wallDudy[ii] = (uCenter[(ii+1) + (yn+1)*(xn+2)] - uCenter[(ii+1) + yn*(xn+2)]) * rdy;
    }

    return 0;
}


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
                                              int xn, int yn, double rdx, double rdy) {
    // Allocate arrays for cell center velocities
    double* uCenter = (double*)calloc(xn * yn, sizeof(double));
    double* vCenter = (double*)calloc(xn * yn, sizeof(double));
    
    if (uCenter == NULL || vCenter == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in calcStaggeredToCellCenterVelocityGradients\n");
        return -1;
    }
    
    // Calculate velocities at cell centers
    calcCellCenterVelocity(u, v, uCenter, vCenter, xn, yn);
    
    // Calculate inverse of 2*dx and 2*dy for central differences
    double r2dx = 0.5 * rdx;
    double r2dy = 0.5 * rdy;
    
    // Calculate velocity gradients at cell centers
    calcCellCenterVelocityGradients(uCenter, vCenter, dudx, dudy, dvdx, dvdy, xn, yn, rdx, rdy);
    
    // Free memory
    free(uCenter);
    free(vCenter);
    
    return 0;
}

/**
 * Function to calculate viscous dissipation at cell centers
 * 
 * @param nu          kinematic viscosity
 * @param dudx        array of ∂u/∂x at cell centers
 * @param dudy        array of ∂u/∂y at cell centers
 * @param dvdx        array of ∂v/∂x at cell centers
 * @param dvdy        array of ∂v/∂y at cell centers
 * @param dissip      array to store viscous dissipation
 * @param totalDissip pointer to store total dissipation
 * @param nx          number of cells in x-direction
 * @param ny          number of cells in y-direction
 * @return            0 on success, -1 on failure
 */
int calcViscousDissipation(double nu, double* dudx, double* dudy, double* dvdx, double* dvdy, double* dissip, double* totalDissip, int nx, int ny) {

    if (dudx == NULL || dudy == NULL || dvdx == NULL || dvdy == NULL || dissip == NULL) {
        fprintf(stderr, "Error: Null pointer provided to calcViscousDissipation\n");
        return -1;
    }

    // Volume of each cell
    double dV = 1.0 / (double)(nx*ny); 

    printf("dV: %lf\n",dV);

    // Initialize total dissipation
    *totalDissip = 0.0;

    // Calculate viscous dissipation for each cell
    for (int jj = 0; jj < ny; jj++) {
        for (int ii = 0; ii < nx; ii++) {

            int idx = ii + jj*nx;

            // Calculate squared terms
            double dudx2 = dudx[idx] * dudx[idx];  // tensor[0]^2
            double dvdx2 = dvdx[idx] * dvdx[idx];  // tensor[1]^2
            double dudy2 = dudy[idx] * dudy[idx];  // tensor[2]^2
            double dvdy2 = dvdy[idx] * dvdy[idx];  // tensor[3]^2

            // Calculate cross terms
            double cross_term = dvdx[idx] * dudy[idx];  // tensor[1]*tensor[3]

            // Calculate dissipation using simplified formula
            dissip[idx] = nu * dV *(2.0*dudx2 + dvdx2 + dudy2 + 2.0*dvdy2 + 2.0*cross_term);
            *totalDissip += dissip[idx];
        }
    }

return 0;
}