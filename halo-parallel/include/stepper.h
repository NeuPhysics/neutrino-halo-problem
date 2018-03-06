//
// Created by Lei Ma on 9/6/17.
//

#ifndef HALO_PARALLEL_STEPPER_H
#define HALO_PARALLEL_STEPPER_H

#include "initializer.h"
#include "helper.h"
#include <cmath>

namespace Stepper
{
   state_type vacuum_euler_forward( const state_type &rho_self ,  const double dt ) {

      

      // state_type hself;
      state_type hvac;

      Hamiltonian::vacuum(hvac);
      
      // for(int i =0;i<3;i++){
      //    hself[i] = h.interaction( rho_counter )[i];
      //    hvac[i] = h.vacuum()[i];
      // } 
     
      state_type rhs;

      for(int i=0;i<3;i++){// Calculate the RHS of \partial_z rho

         int idx1 = (i+1)%3;
         int idx2 = (i+2)%3;

        //  rhs[i] = 2 * ( (hvac[idx2]) * rho_self[idx1] - rho_self[idx2] * (hvac[idx1] ) );
         rhs[i] = 2 * ( (hvac[idx1]) * rho_self[idx2] - rho_self[idx1] * (hvac[idx2] ) );
      }

      for(int i=0;i<3;i++){

         rhs[i] = rhs[i] * dt + rho_self[i];

      }


      return rhs;
   }
   
   void euler_forward( state_type& rhs, const state_type &rho_self , const state_type &rho_counter, const double dt ) {


      state_type hself;
      state_type hvac;

      Hamiltonian::interaction( hself, rho_counter );
      Hamiltonian::vacuum(hvac);
      
      // for(int i =0;i<3;i++){
      //    hself[i] = h.interaction( rho_counter )[i];
      //    hvac[i] = h.vacuum()[i];
      // } 
     
      // state_type rhs;

      for(int i=0;i<3;i++){// Calculate the RHS of \partial_z rho

         int idx1 = (i+1)%3;
         int idx2 = (i+2)%3;

        //  rhs[i] = 2 * ( (hvac[idx2] + hself[idx2]) * rho_self[idx1] - rho_self[idx2] * (hvac[idx1] + hself[idx1]) );
         rhs[i] = 2 * ( (hvac[idx1] + hself[idx1]) * rho_self[idx2] - rho_self[idx1] * (hvac[idx2] + hself[idx2]) );
         // The following code is for {{1,0},{0,1}} density matrix
        //  rhs[i] = ( (hvac[idx1] + hself[idx1]) * rho_self[idx2] - rho_self[idx1] * (hvac[idx2] + hself[idx2]) );
      }

      for(int i=0;i<3;i++){

         rhs[i] = rhs[i] * dt + rho_self[i];

      }


      // return rhs;
   }
   
   void euler_forward_one( state_type &rho_self_future, const state_type &rho_self , const state_type &rho_counter, const double dt, const double reflection = 1.0, const double muf = 5.0, const double costheta = -1.0) {


      state_type hself;
      state_type hvac;

      Hamiltonian::interaction( hself, rho_counter, reflection, muf, costheta );
      Hamiltonian::vacuum(hvac);
      

      for(int i=0;i<3;i++){// Calculate the RHS of \partial_z rho
         /*
            It can be optimized by combining with the last for loop!
         */

         int idx1 = (i+1)%3;
         int idx2 = (i+2)%3;

        //  rho_self_future[i] = 2.0 * ( (hvac[idx2] + hself[idx2]) * rho_self[idx1] - rho_self[idx2] * (hvac[idx1] + hself[idx1]) );
         rho_self_future[i] = 2.0 * ( (hvac[idx1] + hself[idx1]) * rho_self[idx2] - rho_self[idx1] * (hvac[idx2] + hself[idx2]) );
      }

      for(int i=0;i<3;i++){

         rho_self_future[i] = rho_self_future[i] * dt + rho_self[i];

      }


   }
   
   void euler_forward_one_incline( double alpha, state_type &rho_self_future, const state_type &rho_self , const state_type &rho_counter, const double dt, const double reflection = 1.0, const double muf = 5.0, const double costheta = -1.0) {
     // This function also updates a htot and store it in the state assigned

     state_type hself;
     state_type hvac;
     state_type htot;
     double sumhrec = 0.0;
     double eleh = 0.0;

     Hamiltonian::interaction(hself, rho_counter, reflection, muf, costheta);
     Hamiltonian::vacuum(hvac);

     for (int i = 0; i < 3; i++)
     {
       eleh = hself[i] + hvac[i];
       htot[i] = eleh;
       sumhrec += eleh * eleh;
     }
     
     sumhrec = 1.0 / (std::sqrt(sumhrec));

     for (int i = 0; i < 3; i++)
     { // Calculate the RHS of \partial_z rho
       /*
            It can be optimized by combining with the last for loop!
         */

       int idx1 = (i + 1) % 3;
       int idx2 = (i + 2) % 3;

       // rho_self_future[i] = 2 * ( (hvac[idx2] + hself[idx2]) * rho_self[idx1] - rho_self[idx2] * (hvac[idx1] + hself[idx1]) );
       rho_self_future[i] = 2 * ((htot[idx2]) * rho_self[idx1] - rho_self[idx2] * (htot[idx1]));
     }
     double ele = 0.0;
     double inprod = 0.0;
     double sumrec = 0.0;

     for (int i = 0; i < 3; i++)
     {
       ele = rho_self_future[i] * dt + rho_self[i];
       rho_self_future[i] = ele;
       inprod += ele * htot[i];
     }

     int ipsign = Helper::sgnf(inprod);

     for (int i = 0; i < 3; i++)
     {
       ele = ipsign * alpha * htot[i] * sumhrec + (1 - alpha) * rho_self_future[i];
       rho_self_future[i] = ele;
       sumrec += ele * ele;
     }
     sumrec = 1.0 / (std::sqrt(sumrec));

     for (int i = 0; i < 3; i++)
     {
       rho_self_future[i] = rho_self_future[i] * sumrec;
     }
   }
   
   void evolution_op_one( state_type &rho_self_future, const state_type &rho_self , const state_type &rho_counter, const double dt, const double reflection = 1.0, const double muf = 5.0, const double costheta = -1.0) {

     state_type hself;
     state_type hvac;
     state_type htot;
     state_type cpurho;

     Hamiltonian::interaction(hself, rho_counter, reflection, muf, costheta);
     Hamiltonian::vacuum(hvac);

     double sumh = 0.0;
     double eleh = 0.0;
     for (int i = 0; i < 3; i++)
     {
       eleh = hself[i] + hvac[i];
       htot[i] = eleh;
       sumh += eleh * eleh;
      //  std::cout <<  "eleh: " << eleh << std::endl;
      //  std::cout <<  "hself element: " << hself[i] << std::endl;
      //  std::cout <<  "hvac element: " << hvac[i] << std::endl;
     }

     sumh = std::sqrt(sumh);

    //  std::cout <<  "sumh: " << sumh << std::endl;

     state_type u_evo;

     u_evo[0] = -htot[0] / sumh; // The actual vector form of vector h
     u_evo[1] = -htot[1] / sumh;
     u_evo[2] = -htot[2] / sumh;

    //  std::cout <<  "u_evo[0]: " << u_evo[0] << std::endl;
    //  std::cout <<  "u_evo[1]: " << u_evo[1] << std::endl;
    //  std::cout <<  "u_evo[2]: " << u_evo[2] << std::endl;


     double sumurho;
     sumurho = u_evo[0] * rho_self[0] + u_evo[1] * rho_self[1] + u_evo[2] * rho_self[2];
    //  std::cout <<  "sumurho: " << sumurho << std::endl;

     double cs_val = cos(sumh * dt) * sin(sumh * dt);
     double c2_val = cos(2.0 * sumh * dt) ;
     double s_sq_val = 2.0 * pow(sin(sumh * dt), 2.0);
    //  std::cout <<  "cs_val: " << cs_val << std::endl;
    //  std::cout <<  "c2_val: " << c2_val << std::endl;
    //  std::cout <<  "s_sq_val: " << s_sq_val << std::endl;

     cpurho[0] =  2.0 * cs_val * ( rho_self[1] * u_evo[2] - rho_self[2] * u_evo[1] ) ;
     cpurho[1] =  2.0 * cs_val * ( rho_self[2] * u_evo[0] - rho_self[0] * u_evo[2] ) ;
     cpurho[2] =  2.0 * cs_val * ( rho_self[0] * u_evo[1] - rho_self[1] * u_evo[0] ) ;
    //  std::cout <<  "cpurho[0]: " << cpurho[0] << std::endl;
    //  std::cout <<  "cpurho[1]: " << cpurho[1] << std::endl;
    //  std::cout <<  "cpurho[2]: " << cpurho[2] << std::endl;

    //  rho_self_future[0] = rho_self_future[2];
    //  rho_self_future[1] = rho_self_future[0];
    //  rho_self_future[2] = -rho_self_future[1];
       rho_self_future[0] = c2_val * rho_self[0] + s_sq_val * sumurho * u_evo[0] + cpurho[0];
       rho_self_future[1] = c2_val * rho_self[1] + s_sq_val * sumurho * u_evo[1] + cpurho[1];
       rho_self_future[2] = c2_val * rho_self[2] + s_sq_val * sumurho * u_evo[2] + cpurho[2];

    //  std::cout <<  "rho_self_future[0]: " << rho_self_future[0] << std::endl;
    //  std::cout <<  "rho_self_future[1]: " << rho_self_future[1] << std::endl;
    //  std::cout <<  "rho_self_future[2]: " << rho_self_future[2] << std::endl;
   }


   void euler_forward_one_nunubar( state_type &rho_self_future, const state_type &rho_self , const state_type &rho_counter, const state_type &rho_ya_counter, const state_type &rho_same_direction, const double dt, const double spectrum[2], const double reflection[2], const double muf[2], const double costheta[3]) {
      // spectrum[2] = { the corresponding self spectrum which are are updating, spectrum of the other beam }; The order can not be switched since the vacuum Hamiltonian is calculated based on spectrum[0].
      // interaction_nunubar( state_type &h_store, const state_type &rho_counter, const state_type &rho_ya_counter, const state_type &rho_same_direction, const double spectrum[2], const double reflection[2], const double muf[2], const double costheta[3] )

      state_type hself;
      state_type hvac;

      Hamiltonian::interaction_nunubar(hself, rho_counter, rho_ya_counter, rho_same_direction, spectrum, reflection, muf, costheta);
      double eta = spectrum[0]/abs(spectrum[0]);
      Hamiltonian::vacuum(hvac, eta);

      for(int i=0;i<3;i++){// Calculate the RHS of \partial_z rho

         int idx1 = (i+1)%3;
         int idx2 = (i+2)%3;

         rho_self_future[i] = 2 * ( (hvac[idx1] + hself[idx1]) * rho_self[idx2] - rho_self[idx1] * (hvac[idx2] + hself[idx2]) );
      }

      for(int i=0;i<3;i++){

         rho_self_future[i] = rho_self_future[i] * dt + rho_self[i];

      }


   }

}


#endif //HALO_PARALLEL_STEPPER_H
