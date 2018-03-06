//
// Created by Lei Ma on 9/6/17.
//

/*
Parameters and intializations are set here
*/

#ifndef HALO_PARALLEL_INITIALIZER_H
#define HALO_PARALLEL_INITIALIZER_H

#include <iostream>
#include <algorithm>
#include <iterator>
#include <complex>
#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>

using namespace std;

/*
State contains three elements. It can be used to describe a 2 by 2 matrix with determined trace.
In the case of traceless density matrix, state_type[0] is the 1 1 element which has to be real, 
state_type[1] is the real part of 12 element, state_type[2] is imaginary part of 12 element.
*/
typedef array< double , 3 > state_type;

// Define lambda and omega

// const double s2theta = 0.916515138991; // sin^2 2theta;theta_{12} with sin^2(theta)=0.3
const double s2theta = 0.0; // sin^2 2theta;theta_{12} with sin^2(theta)=0.3
//const double s2theta = 0.01; // sin^2 2theta;theta_{12} with sin^2(theta)=0.3
//const double s2theta = 0.2; // sin^2 2theta;theta_{12} with sin^2(theta)=0.3
const double c2theta = sqrt( 1.0 - pow( s2theta, 2.0) );

const double omegav = -1.0; // scale all quantities using omega, now distance is x omega
//const double omegav = -1.0; // scale all quantities using omega, now distance is x omega
//const double omegav = 0.0; // scale all quantities using omega, now distance is x omega
const double refl = 0.0; // 0.01; // reflected neutrinos, the backward beam is always reduced by this number, thus the actual mu we use is multiplied by this number
const double mu = 5.0; // 50.0 * omegav;// * refl; // in unit of omegav where omegav = 1 // 2017-10-13: Versions before this one was using the convention mu=50.0 * omegav * refl; 
// const double mu = 4.0; // 50.0 * omegav;// * refl; // in unit of omegav where omegav = 1 // 2017-10-13: Versions before this one was using the convention mu=50.0 * omegav * refl; 


#endif
