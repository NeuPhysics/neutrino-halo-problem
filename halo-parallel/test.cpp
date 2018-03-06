/*
This file is used to test the headers.
*/



#include "include/state.h"
#include "include/hamiltonian.h"
#include "include/stepper.h"
#include "include/looper.h"
#include "include/recorder.h"
#include "include/helper.h"
// #include "include/loader.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <typeinfo>

using namespace std;

int main() {

    const int SIZE = 10;
   StateArray *yasaptr = new StateArray(SIZE);
   StateArray *saptr = new StateArray(SIZE);
  

   double init_value[3] = {1.0,0.0,0.0};
   (*yasaptr).init(init_value);
   (*saptr).init(init_value);
   
   (*yasaptr)[SIZE-2][0] = 0.9;
   (*yasaptr)[SIZE-2][1] = 0.43;

   cout << "Testing Evolution Operator Stepper: " << endl; 
   for(int i=0 ; i<10*SIZE; i++) {    
   Looper::halo_evolution_op_one(saptr, yasaptr, 1.0/SIZE, SIZE, 0.1, 5.0);
   Recorder::cout_state( saptr, SIZE);
   Recorder::cout_state( yasaptr, SIZE);
        Helper::sa_ptr_swap(saptr, yasaptr);
   }
    // Stepper::evolution_op_one( (*saptr)[1], (*yasaptr)[0], (*yasaptr)[SIZE-2], 1.0/SIZE, 0.1, 5.0, -1.0 );


   cout << "BEFORE AVG" << endl;
   cout << (*saptr)[0][0] << (*saptr)[0][1] <<  (*saptr)[0][2]<< endl;
   cout << (*yasaptr)[0][0] << (*yasaptr)[0][1] << (*yasaptr)[0][2] << endl;

   Helper::state_avg(saptr, yasaptr, SIZE, 0.0);

   cout << "After AVG: alpha =0.0 " << endl;
   cout << (*saptr)[0][0] << (*saptr)[0][1] <<  (*saptr)[0][2]<< endl;
   cout << (*yasaptr)[0][0] << (*yasaptr)[0][1] << (*yasaptr)[0][2] << endl;

   Helper::state_avg(saptr, yasaptr, SIZE, 0.5);

   cout << "After AVG: alpha=0.5 " << endl;
   cout << (*saptr)[0][0] << (*saptr)[0][1] <<  (*saptr)[0][2]<< endl;
   cout << (*yasaptr)[0][0] << (*yasaptr)[0][1] << (*yasaptr)[0][2] << endl;

   (*saptr)[0][0] = 1.0;
   (*yasaptr)[0][0] = 10.0;
   cout << (*saptr)[0][0] << endl;
   cout << (*yasaptr)[0][0] << endl;
   
//   printf("Address of sa is %p\n", sa); 
//   printf("Address of yasa is %p\n", yasa); 
   printf("saptr is %p\n", (void*)saptr); 
   printf("yasaptr is %p\n", (void*)yasaptr); 
   Helper::sa_ptr_swap(saptr, yasaptr);
   printf("Address of sa is %p\n", (void*)saptr); 
   printf("Address of yasa is %p\n", (void*)yasaptr); 
//   printf("saptr is %p\n", saptr); 
//   printf("yasaptr is %p\n", yasaptr); 

   cout << (*saptr)[0][0] << endl;
   cout << (*yasaptr)[0][0] << endl;


   (*yasaptr)[0][0] = 2.0 ;

   cout << (*saptr)[0][0] << endl;
   cout << (*yasaptr)[0][0] << endl;
        

   // Looper::halo_euler_forward_one(saptr, yasaptr, 0.01, 1000);

//   Stepper::euler_forward_one( (*saptr)[1], (*yasaptr)[0], (*saptr)[0], 0.001);

    cout << "PROGRAM END" << endl;

    return 0;
}
