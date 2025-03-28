#ifndef CONFIG_H_
#define CONFIG_H_

// Configuration structure for simulation parameters
typedef struct Config {
    char caseName[50];         // Name of the test case
    int stepNum;               // Number of simulation steps
    int outputInterval;        // Interval for output generation
    double deltaT;             // Time step size
    int nx;                    // Grid size in x direction
    int ny;                    // Grid size in y direction
    double re;                 // Reynolds number
    int withInitialCondition;  // Flag for initial condition (0 or 1)
    double convergenceThreshold; // Threshold for convergence check
    int writePrecision;        // Precision for output data
} Config;

#endif 