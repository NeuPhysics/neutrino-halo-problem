//
// Created by Lei Ma on 9/8/17.
//

#ifndef HALO_PARALLEL_HELPER_H
#define HALO_PARALLEL_HELPER_H

#include "initializer.h"
#include <ctime>

#include <omp.h> // for openmp

namespace Helper
{
   
   double sgnf(double val) {
      return val/std::abs(val);
   } 

   // Reverse an array of state_type and store in a new array; Tested.
   void state_array_reverse( const state_type state [], state_type state_reversed [], const int length){
   
      for(int i=0; i< length; i++){
         
         for(int j=0; j < 3; j++) {
            state_reversed[i][j] = state[length -1 - i][j];
         }

      }

      // return state_reversed;

   }

   void state_array_copy( const state_type state [], state_type state_copied [], const int length) {


      for(int i=0; i< length; i++){
         
         for(int j=0; j < 3; j++) {
            state_copied[i][j] = state[i][j];
         }

      }

   }

   void sa_ptr_swap( StateArray* &stateptr, StateArray* &stateptr_swap) {

      StateArray *saptr_temp = stateptr;
      stateptr = stateptr_swap;
      stateptr_swap = saptr_temp;

   }  
   
   void state_avg( StateArray* statearrayavg, StateArray* statearray, const int length, double alpha = 0.5) {

      double norm_avg = 1.0;
      double ele = 0.0;
      double sumrecp;

      #pragma omp parallel for
      for(int i=0; i< length; i++) {
            sumrecp = 0.0; // Have to set it to 0 for each iteration in i
            for(int j=0;j < 3;j++) {
                  ele = (alpha) * (*statearray)[i][j] + (1 - alpha) * (*statearrayavg)[i][j];
                  (*statearrayavg)[i][j] = ele; 
                  sumrecp = sumrecp + ele*ele;
            }
           sumrecp = 1/( norm_avg*std::sqrt(sumrecp) ); 
            for(int j=0;j < 3;j++) {
                 (*statearrayavg)[i][j] = (*statearrayavg)[i][j] *sumrecp;
            }
      }
   }  


   double squared_difference( StateArray *current_state, StateArray *past_state, const int length, const int savestep = 1 ) {
      
         double sum = 0;
      
         for(int i = 0; i < length; i+=savestep) {
            sum += pow( (*current_state)[i][0] - (*past_state)[i][0], 2.0 );
            // sum += abs( (*current_state)[i][0] - (*past_state)[i][0] ); // Calculation of absolute value is much much slower than square.
         }
         
         return sum;
      
   }
      

}

// End of Helper namespace



class Timing {
public:
   Timing() :rawtime{time(0)} {} //Constructor
   string timestamp()
   {
      // time (&rawtime);
      tm *ctm = localtime(&rawtime);
      int year = ctm->tm_year + 1900;
      int month = ctm->tm_mon+1;
      int day = ctm->tm_mday;
      int hour = ctm->tm_hour;
      int min = ctm->tm_min;
      int sec = ctm->tm_sec;

      return to_string(year) + "-" + to_string(month)+ "-"+to_string(day)+ "-" + to_string(hour) + "-"+ to_string(min) + "-" + to_string(sec);
   }

   string time_pretty()
   {
      return ctime (&rawtime);
   }

   clock_t stopwatch_reset()
   {
      const clock_t stopwatch_begin = clock();
      return stopwatch_begin;
   }

   void stopwatch_info(const clock_t stopwatch_begin_time, const int NIterations, const double calc_iter=1000.0)
   {
      const clock_t stopwatch_end = clock();
      cout << "Total clock time: " << float( stopwatch_end - stopwatch_begin_time ) /  CLOCKS_PER_SEC << endl;
      cout << "Clock time for 1000 iterations: " << calc_iter * float( stopwatch_end - stopwatch_begin_time  ) /  CLOCKS_PER_SEC / NIterations << endl;
   }
   
   clock_t wall_stopwatch_reset()
   {
      const clock_t stopwatch_begin = omp_get_wtime();
      return stopwatch_begin;
   }

   void wall_stopwatch_info(const clock_t stopwatch_begin_time, const int NIterations, const double calc_iter=1000.0)
   {
      const clock_t stopwatch_end = omp_get_wtime();
      cout << "Total clock time: " << float( stopwatch_end - stopwatch_begin_time ) << endl;
      cout << "Clock time for 1000 iterations: " << calc_iter * float( stopwatch_end - stopwatch_begin_time  ) / NIterations << endl;
   }

private:
   time_t rawtime;

};



#endif
