//
// Created by Lei Ma on 9/6/17.
//

#ifndef HALO_PARALLEL_HAMILTONIAN_H
#define HALO_PARALLEL_HAMILTONIAN_H

#include "initializer.h"


namespace Hamiltonian
{


   void vacuum(state_type & hvac, double spectrum = 1.0){

      const double omegav_2 = spectrum * omegav * 0.5;

      // hvac[0] = -c2theta * omegav_2;
      // hvac[1] = s2theta * omegav_2;
      // hvac[2] = 0.0;
      
      hvac[2] = -c2theta * omegav_2;
      hvac[0] = s2theta * omegav_2;
      hvac[1] = 0.0;
   }

   void interaction( state_type& hnn, const state_type &rho_counter, const double reflection = 1.0, const double muf = 5.0, const double costheta = -1.0 ){

      // state_type hnn;

      for(int i=0;i<3;i++){
         hnn[i] = (1.0-costheta) * muf * reflection * rho_counter[i];
      }

      // return hnn;
   }

   // Neutrino interactions with antineutrinos
   void interaction_nunubar( state_type &h_store, const state_type &rho_counter, const state_type &rho_ya_counter, const state_type &rho_same_direction, const double spectrum[2], const double reflection[2], const double muf[2], const double costheta[3]){
      // muf[2] = {mu_self,mu_the_other}
      // spectrum = {spect_self, spect_the_other}
      // reflection[2] = { reflection on the other side as the self beam, reflection on the same side of the self beam }
      // The two dimensional arrays: first is for neutrino, second element is for antineutrino
      // reflection can be different for different flavors
      // costheta[3] = { cos(theta_self+theta_self) , cos(theta_other+theta_self), cos(theta_other - theta_self) }, 

      // For actual halo problem
      // for(int i=0;i<3;i++){
      //    h_store[i] = ( 1.0 - costheta[0] )* spectrum[0] * muf[0] * reflection[0] * rho_counter[i] + ( 1.0 - costheta[1] ) * spectrum[1] * muf[1] * reflection[0] * rho_ya_counter[i] + ( 1.0 - costheta[2] )* spectrum[1] * muf[1] * reflection[1] * rho_same_direction[i] ;
      // }

      // To test bipolar model, I change the rho_ya_counter reflection to reflection[1] which is set to 1
      for(int i=0;i<3;i++){
         h_store[i] = ( 1.0 - costheta[0] )* spectrum[0] * muf[0] * reflection[1] * rho_counter[i] + ( 1.0 - costheta[1] ) * spectrum[1] * muf[1] * reflection[1] * rho_ya_counter[i] + ( 1.0 - costheta[2] )* spectrum[1] * muf[1] * reflection[0] * rho_same_direction[i] ;
      }
      // return hnn;
   }


}

#endif //HALO_PARALLEL_HAMILTONIAN_H
