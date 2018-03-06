/*
halo_sim is meant to solve the halo problem using a simultaneous algrimth.
*/

#include "include/state.h"
#include "include/hamiltonian.h"
#include "include/stepper.h"
#include "include/looper.h"
#include "include/recorder.h"
#include "include/helper.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

// #include <omp.h> //For openmp

using namespace std;

int main(int argc, char *argv[]) { // three arguments: { # of total iteractions, # of steps within one calculations, # of threads  }

     /* 
        Check the number of parameters
     */
     if (argc < 4) {
        // Tell the user how to run the program
        cerr << "Usage: " << argv[0] << " takes 3 parameters:  { # of total iteractions, # of steps within one calculations, # of threads  }" << endl;
        // "Usage messages" are a conventional way of telling the user
        return 1;
    }

    /* 
        Set openmp number of threads
    */
  //  omp_set_dynamic(0);
  //  omp_set_num_threads( stoi(argv[3]) );


    /*
        Initializing some constants
    */
    const int Ntop = stoi(argv[1]); // Set how many times the overall iteraction is
    const int STEPS= stoi(argv[2]); // Set size of the array of z in z direction
    const int SIZE = 2*STEPS; // Set the size of the array for all states
    const double range0 = 0.0; // Set initial coordinate for z
    const double range1 = 1.0; // Set end of z range
    const double STEPSIZE = ( range1 - range0 )/SIZE; // Calculate step size of the calculation
    const int STEPSAVE = Ntop/50; // How many overall iteractions before saving data
    int saveflag = 0; // indicator for the save cycles used to indicate STEPSAVE
    const int STEPSAVEZ = 100; // How many steps before save a data point within each iteraction

    /* 
        data for saving
    */
    ofstream data("halo_parallel_" + std::to_string(SIZE)+"_"+ std::to_string(STEPSIZE) + ".csv" ); 


    /*
        Write basic information about the computation
    */
    cout << "PROGRAM START" << endl;
    cout << "Halo Problem Forward and Backward:" << endl << "Total number of iterations: " + to_string(Ntop) << endl << "Size of rhos: " + to_string(SIZE) << endl << "Range: " + to_string(range1) << endl << "Step size: " + to_string(STEPSIZE) << endl; 
    cout << "Save Steps: " << to_string(STEPSAVE) << endl;

    /* 
        Output warning if the step size and range and size of array are not consistant
    */
    if( (SIZE*STEPSIZE + range0) - range1 > 0.1 ){
        cout << "NOT REACHING THE END POINT!" << endl;
    }
    


    state_type rho_init = { 1.0, 0.0, 0.0 };

    /*
        The state arrays
    */
    StateArray sa_store(SIZE);
    StateArray sa(SIZE);
    StateArray* sa_store_ptr = &sa_store;
    StateArray* sa_ptr = &sa;

    sa.init();
    sa_store.init();

    /*
        Assign initial values to sa_store which is used to calculate sa
    */
    for(auto i = 0; i < 3; i++) {
        sa_store[0][i] = rho_init[i]; 
    }


    // Track clock time
    const clock_t begin_time = clock();


    for(int i=0; i<Ntop; i++ ){

        if(saveflag%STEPSAVE == 0) {
            Recorder::record_state( sa_ptr, SIZE, data , STEPSAVEZ );
        }

        Stepper::euler_forward_one( (*rho_array_ptr)[i+1], (*rho_array_store_ptr)[i], (*rho_array_ptr)[totallength - 2  - i], dt) ;

        Helper::sa_ptr_swap(sa_ptr, sa_store_ptr);


        ++saveflag;

    }


    /* 
        Calculate clock time for the iterations
    */
    cout << "Total clock time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    cout << "Clock time for 1000 iterations: " << 1000.0 * float( clock () - begin_time ) /  CLOCKS_PER_SEC / Ntop << endl;

    // Output data for benchmark
    // cout << to_string(Ntop) << ", " <<  to_string(SIZE) << ", " << argv[3] << ", " << float( clock () - begin_time ) /  CLOCKS_PER_SEC / Ntop << endl; // Ntop, SIZE, Threads, Time/Iteration



    cout << "PROGRAM END" << endl;

    return 0;
}
