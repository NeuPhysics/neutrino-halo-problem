//
// Created by Lei Ma on 9/6/21.
//

#ifndef HALO_PARALLEL_STATE_H
#define HALO_PARALLEL_STATE_H

#include "initializer.h"

class StateArray {
public:
   StateArray(const int s) :elem{new state_type[s]}, sz{s} {} // Constructor
   state_type& operator[](const int i)  // getter
   { 
      return elem[i];
   }
   int size() { return sz;}  // return size of the array of states
   void init(const double iv[3]) // initialize the array to 0's
   {
      for(int i=0;i<sz;i++) { // Might need long int when the array size gets huge!
         elem[i][0] = iv[0];
         elem[i][1] = iv[1];
         elem[i][2] = iv[2];
      }
   }
   
   ~StateArray()
   {
     delete[] elem;
   }
private:
   state_type* elem;
   int sz;
};


#endif //HALO_PARALLEL_STATE_H
