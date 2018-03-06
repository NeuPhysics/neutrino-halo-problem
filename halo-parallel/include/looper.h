//
// Created by Lei Ma on 9/6/17.
//

#ifndef HALO_PARALLEL_LOOPER_H
#define HALO_PARALLEL_LOOPER_H

#include "initializer.h"
#include "stepper.h"
#include "helper.h"
#include <omp.h> //For openmp

// #include "recorder.h"// For test

namespace Looper
{

   void vacuum_euler_forward( state_type rho_self_array [], const double dt, const int length){ // loop through for iter iterations



      for(int i =0; i<length-1; i++){
         
        rho_self_array[i+1] =  Stepper::vacuum_euler_forward(rho_self_array[i],  dt) ;

      }
   }

   void interaction_euler_forward( state_type rho_self_array [], state_type rho_counter_array [], const double dt, const int length){ // loop through for iter iterations


      // state_type rhs;

      for(int i =0; i<length-1; i++){
         
        Stepper::euler_forward(rho_self_array[i+1], rho_self_array[i], rho_counter_array[length - 1 - i], dt) ;
        // rho_self_array[i+1] = rhs;

      }
   }
   
   void halo_euler_forward( StateArray rho_forward_array, StateArray rho_backward_array, const double dt, const int length){ // loop through for iter iterations

      
      state_type rhs;

      for(int i =0; i<length-1; i++){
         
        Stepper::euler_forward(rhs, rho_forward_array[i], rho_backward_array[length - 1 - i], dt) ;
        rho_forward_array[i+1] =  rhs; 
        Stepper::euler_forward(rhs, rho_backward_array[i], rho_forward_array[length - 1 - i], dt) ;
        rho_backward_array[i+1] = rhs;

      }
   }




   void halo_euler_forward_one(StateArray* rho_array_ptr, StateArray* rho_array_store_ptr, const double dt, const int totallength, const double reflection = 1, const double muf = 5.0, const double costheta = -1.0){ // loop through for iter iterations
     
      int length = totallength/2;

      #pragma omp parallel for
      for(int i =0; i<length-1; i++){
         
        Stepper::euler_forward_one( (*rho_array_ptr)[i+1], (*rho_array_store_ptr)[i], (*rho_array_store_ptr)[totallength - 2  - i], dt, reflection, muf, costheta) ;
        Stepper::euler_forward_one( (*rho_array_ptr)[length + i], (*rho_array_store_ptr)[length -1 + i], (*rho_array_store_ptr)[length -1 - i], dt, 1.0, muf, costheta) ;
      }
   }

   void halo_euler_forward_one_avg(StateArray* rho_array_ptr, StateArray* rho_array_store_ptr, const double dt, const int totallength, const double alpha, const double reflection = 1, const double muf = 5.0, const double costheta = -1.0){ // loop through for iter iterations
     
      int length = totallength/2;

      #pragma omp parallel for
      for(int i =0; i<length-1; i++){
         
        Stepper::euler_forward_one( (*rho_array_ptr)[i+1], (*rho_array_store_ptr)[i], (*rho_array_store_ptr)[totallength - 2  - i], dt, reflection, muf, costheta) ;
        Stepper::euler_forward_one( (*rho_array_ptr)[length + i], (*rho_array_store_ptr)[length -1 + i], (*rho_array_store_ptr)[length -1 - i], dt, 1.0, muf, costheta) ;

        double sumrecpf = 0.0;
        double sumrecpb = 0.0;
        double elef = 0.0;
        double eleb = 0.0;
        // Average the new results with old results
        for(int j=0; j<3;j++){
            elef = alpha * (*rho_array_store_ptr)[i][j] + (1 - alpha) * (*rho_array_ptr)[i][j];
            (*rho_array_ptr)[i][j] = elef;
            eleb= alpha * (*rho_array_store_ptr)[length + i][j] + (1 - alpha) * (*rho_array_ptr)[length + i][j];
            (*rho_array_ptr)[length + i][j] = eleb;
            sumrecpf = sumrecpf + elef*elef;
            sumrecpb = sumrecpb + eleb*eleb;
        }
           sumrecpf = 1/( std::sqrt(sumrecpf) ); 
           sumrecpb = 1/( std::sqrt(sumrecpb) ); 
            for(int j=0;j < 3;j++) {
                (*rho_array_ptr)[i][j] = (*rho_array_ptr)[i][j] * sumrecpf;
                (*rho_array_ptr)[length + i][j] = (*rho_array_ptr)[length + i][j] * sumrecpb;
            }

      }
   }

   void halo_euler_forward_one_incline(StateArray *rho_array_ptr, StateArray *rho_array_store_ptr, const double dt, const int totallength, const double alpha, const double reflection = 1, const double muf = 5.0, const double costheta = -1.0)
   { // loop through for iter iterations

     int length = totallength / 2;
     double alpha_rescaled = alpha/length;

     #pragma omp parallel for
     for (int i = 0; i < length - 1; i++)
     {

      //  state_type hamilf;
      //  state_type hamilb;

       // Stepper::euler_forward_one_w_h has been validated and compared to previous results.
       Stepper::euler_forward_one_incline(alpha_rescaled, (*rho_array_ptr)[i + 1], (*rho_array_store_ptr)[i], (*rho_array_store_ptr)[totallength - 2 - i], dt, reflection, muf, costheta);
       Stepper::euler_forward_one_incline(alpha_rescaled, (*rho_array_ptr)[length + i], (*rho_array_store_ptr)[length - 1 + i], (*rho_array_store_ptr)[length - 1 - i], dt, 1.0, muf, costheta);

      //  int ipfsign = Helper::sgnf(innerproductf);
      //  int ipbsign = Helper::sgnf(innerproductb);
     
     }
   }

   void halo_evolution_op_one(StateArray* rho_array_ptr, StateArray* rho_array_store_ptr, const double dt, const int totallength, const double reflection = 1.0, const double muf = 5.0, const double costheta = -1.0){ // loop through for iter iterations
     
      int length = totallength/2;

      #pragma omp parallel for
      for(int i =0; i<length-1; i++){
         
        Stepper::evolution_op_one( (*rho_array_ptr)[i+1], (*rho_array_store_ptr)[i], (*rho_array_store_ptr)[totallength - 2  - i], dt, reflection, muf, costheta) ;
        Stepper::evolution_op_one( (*rho_array_ptr)[length + i], (*rho_array_store_ptr)[length -1 + i], (*rho_array_store_ptr)[length -1 - i], dt, 1.0, muf, costheta) ;
      }
   }
   
   void halo_euler_forward_one_nunubar(StateArray* rho_array_ptr, StateArray* rho_array_store_ptr, StateArray* rho_another_array_ptr, StateArray* rho_another_array_store_ptr, const double dt, const int totallength, const double spectrum[2], const double reflection, const double mu_arr[2], const double costheta[4]){ // loop through for iter iterations
    
    // spectrum[2] = {left beam, right beam}
     // mu_arr = {mu_left, mu_right} or {mu_1,mu_2} for the two beams
     // I define {mu_self, mu_the_other} when doing the calculations since the Hamiltonian takes in such parameters.
     // costheta[4] = { cos(2theta_left), cos(2theta_right) , cos(theta_right- theta_left), cos (theta_right+theta_left) }

      int length = totallength/2;
      // when calculating the right beam, the order of spectrum is reversed so I define the reversed spectrum
      double spectrum_r[2]; 
      double mu_arr_r[2]; 
      spectrum_r[0] = spectrum[1]; // spectrum_r is for the calculation of the right beam
      spectrum_r[1] = spectrum[0];
      mu_arr_r[0] = mu_arr[1];
      mu_arr_r[1] = mu_arr[0];
      // define the costheta's needed for each beam
      double costheta_l[3];
      double costheta_r[3];

      costheta_l[0] = costheta[0];
      costheta_l[1] = costheta[3];
      costheta_l[2] = costheta[2];

      costheta_r[0] = costheta[1];
      costheta_r[1] = costheta[3];
      costheta_r[2] = costheta[2];

      // build the reflection array
      double refl_arr_f[2];
      refl_arr_f[1] = reflection;
      refl_arr_f[0] = 1.0;
      double refl_arr_b[2];
      refl_arr_b[1] = 1.0;
      refl_arr_b[0] = reflection;

      // interaction_nunubar( state_type &h_store, const state_type &rho_counter, const state_type &rho_ya_counter, const state_type &rho_same_direction, const double spectrum[2], const double reflection[2], const double muf[2], const double costheta[3] )
      #pragma omp parallel for
      for(int i =0; i<length-1; i++){
         
        // the left beam forward
        Stepper::euler_forward_one_nunubar( (*rho_array_ptr)[i+1], (*rho_array_store_ptr)[i], (*rho_array_store_ptr)[totallength - 2  - i], (*rho_another_array_store_ptr)[totallength - 2  - i], (*rho_another_array_store_ptr)[i], dt, spectrum, refl_arr_f, mu_arr, costheta_l) ;
        // right beam forward
        Stepper::euler_forward_one_nunubar( (*rho_another_array_ptr)[i+1], (*rho_another_array_store_ptr)[i],  (*rho_another_array_store_ptr)[totallength - 2  - i], (*rho_array_store_ptr)[totallength - 2  - i], (*rho_array_store_ptr)[i], dt, spectrum_r, refl_arr_f, mu_arr_r, costheta_r) ;
        // left beam backward: left means it's the continuation of the original left beam, which is stored in the same array // Comment out to test bipolar model
        Stepper::euler_forward_one_nunubar( (*rho_array_ptr)[length + i], (*rho_array_store_ptr)[length -1 + i], (*rho_array_store_ptr)[length -1 - i], (*rho_another_array_store_ptr)[length - 1  - i], (*rho_another_array_store_ptr)[length-1+i], dt, spectrum, refl_arr_b, mu_arr, costheta_l);
        // right beam backward // Comment out to test bipolar model
        Stepper::euler_forward_one_nunubar( (*rho_another_array_ptr)[length + i], (*rho_another_array_store_ptr)[length -1 + i], (*rho_another_array_store_ptr)[length - 1  - i], (*rho_array_store_ptr)[length -1 - i],  (*rho_array_store_ptr)[length-1+i], dt, spectrum_r, refl_arr_b, mu_arr_r, costheta_r);
      }
   }

  void halo_euler_forward_one_bipolar(StateArray* rho_array_ptr, StateArray* rho_array_store_ptr, StateArray* rho_another_array_ptr, StateArray* rho_another_array_store_ptr, const double dt, const int totallength, const double spectrum[2], const double reflection, const double mu_arr[2], const double costheta[4]){ // loop through for iter iterations
    
    // spectrum[2] = {left beam, right beam}
     // mu_arr = {mu_left, mu_right} or {mu_1,mu_2} for the two beams
     // I define {mu_self, mu_the_other} when doing the calculations since the Hamiltonian takes in such parameters.
     // costheta[4] = { cos(2theta_left), cos(2theta_right) , cos(theta_right- theta_left), cos (theta_right+theta_left) }

      int length = totallength/2;
      // when calculating the right beam, the order of spectrum is reversed so I define the reversed spectrum
      double spectrum_r[2]; 
      double mu_arr_r[2]; 
      spectrum_r[0] = spectrum[1]; // spectrum_r is for the calculation of the right beam
      spectrum_r[1] = spectrum[0];
      mu_arr_r[0] = mu_arr[1];
      mu_arr_r[1] = mu_arr[0];
      // define the costheta's needed for each beam
      double costheta_l[3];
      double costheta_r[3];

      costheta_l[0] = costheta[0];
      costheta_l[1] = costheta[3];
      costheta_l[2] = costheta[2];

      costheta_r[0] = costheta[1];
      costheta_r[1] = costheta[3];
      costheta_r[2] = costheta[2];

      // build the reflection array
      double refl_arr_f[2];
      refl_arr_f[1] = reflection;
      refl_arr_f[0] = 1.0;
      double refl_arr_b[2];
      refl_arr_b[1] = 1.0;
      refl_arr_b[0] = reflection;

      // interaction_nunubar( state_type &h_store, const state_type &rho_counter, const state_type &rho_ya_counter, const state_type &rho_same_direction, const double spectrum[2], const double reflection[2], const double muf[2], const double costheta[3] )
      #pragma omp parallel for
      for(int i =0; i<length-1; i++){
        // the left beam forward
        Stepper::euler_forward_one_nunubar( (*rho_array_ptr)[i+1], (*rho_array_store_ptr)[i], (*rho_array_store_ptr)[totallength-2-i], (*rho_another_array_store_ptr)[totallength - 2  - i], (*rho_another_array_store_ptr)[i], dt, spectrum, refl_arr_f, mu_arr, costheta_l) ;
        // right beam forward
        Stepper::euler_forward_one_nunubar( (*rho_another_array_ptr)[i+1], (*rho_another_array_store_ptr)[i],  (*rho_another_array_store_ptr)[totallength-2-i], (*rho_array_store_ptr)[totallength - 2  - i], (*rho_array_store_ptr)[i], dt, spectrum_r, refl_arr_f, mu_arr_r, costheta_r) ;
      }
   }

}

#endif //HALO_PARALLEL_LOOPER_H
