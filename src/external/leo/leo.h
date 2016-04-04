/*
 * runEM.h
 *
 * Code generation for function 'runEM4'
 *
 * C source code generated on: Thu Jun  5 13:30:30 2014
 *
 */

#ifndef __LEO_H__
#define __LEO_H__
/* Include files */

#include <iostream>
#include <math.h>
//#include <random>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <armadillo>
#include <list>
#include <omp.h>

#define EPSILON 70000
#define ITERATION_LIMIT 100
#define INF 1000000000
/* Function Declarations */

typedef struct{
    arma::vec bandwidthPredictions;
    arma::vec powerPredictions;
    double bandwidthError;
    double powerError;
}PredictionResults;

/**
 * @param appId Between [0, N - 1] (N = number of applications in the file).
 **/
PredictionResults compute(int appId, char* powerFile, char* bandwidthFile, arma::vec* sampledPower, arma::vec* sampledPerformance);

#endif
/* End of code generation (runEM.h) */
