//
// Created by Lei Ma on 9/6/17.
//

#ifndef HALO_PARALLEL_RECORDER_H
#define HALO_PARALLEL_RECORDER_H

#include <iostream>
#include <fstream>
#include "initializer.h"

using namespace std;

namespace Recorder
{
   void cout_state( StateArray *st_array, const int length, const int savestep = 1, const int idx = 2) {
      
      for(int i = 0; i < length-1; i+=savestep) {
         cout << (*st_array)[i][idx] << ", ";
      }

      cout << (*st_array)[length-1][idx] << endl;

   }

   void record_z( double *zarray, const int length, ofstream& dumpObj, const int savestep = 1 ) {
      
      // int length = sizeof(zarray)/8;

      for(int i=0;i < length-1; i += savestep ){
         dumpObj << zarray[i] << ", " ;
      }

      dumpObj << zarray[length-1] << endl;

   }

   void record_state( StateArray *st_array, const int length, ofstream& dumpObj, const int savestep = 1, const int idx = 2) {
      
      for(int i = 0; i < length-1; i+=savestep) {
         dumpObj << (*st_array)[i][idx] << ", ";
      }

      dumpObj << (*st_array)[length-1][idx] << endl;

   }
   
   void cout_single_state( StateArray *st_array, const int position) {
      
      cout << (*st_array)[position][0] << ", ";
      cout << (*st_array)[position][1] << ", ";
      cout << (*st_array)[position][2] << endl;

   }


}

#endif //HALO_PARALLEL_RECORDER_H
