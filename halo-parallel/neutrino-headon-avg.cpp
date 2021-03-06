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

#include <omp.h> //For openmp

using namespace std;

int main(int argc, char *argv[]) { // six arguments: { range, reflection coefficient, neutrino potential mu, # of total iteractions, # of steps within one calculations, # of threads }

     /* 
        Check the number of parameters
     */
     if (argc < 7) {
        // Tell the user how to run the program
        cerr << "Usage: " << argv[0] << " takes 6 parameters:  { range, reflection coefficient, neutrino potential mu, # of total iteractions, # of steps within one calculations, # of threads }" << endl;
        // "Usage messages" are a conventional way of telling the user
        return 1;
    }

    /* 
        Set openmp number of threads
    */
    omp_set_dynamic(0);
    omp_set_num_threads( stoi(argv[6]) );


    /*
        Initializing some constants
    */
    const double RCOEFF = stod(argv[2]); // Set reflection coefficient
    const double MIU = stod(argv[3]); // Set reflection coefficient
    const int Ntop = stoi(argv[4]); // Set how many times the overall iteraction is
    const int STEPS = stoi(argv[5]); // Set size of the array of z in z direction
    const int TH = stoi(argv[6]); // Set number of threads

    const int SIZE = 2*STEPS; // Set the size of the array for all states
    const double range0 = 0.0; // Set initial coordinate for z
    const double range1 = stod(argv[1]); // Set end of z range
    const double STEPSIZE = ( range1 - range0 )/STEPS; // Calculate step size of the calculation
    const int STEPSAVE = Ntop/10000; // How many overall iteractions before saving data
    int saveflag = 0; // indicator for the save cycles used to indicate STEPSAVE
    int termflag = 0; // indicator for the termination of the calculation
    const int STEPSAVEZ = SIZE/200; // How many steps before save a data point within each iteraction

    const double avg_alpha = 0.0; // within [0,1]

    /*
        Timing object
    */
    Timing tmg;


    /* 
        data for saving
    */
    ofstream data("neutrino-headon-avg-" + to_string(avg_alpha).substr(0,3) + "_MU_" + to_string(MIU) +"_REFL" + to_string(RCOEFF) + "_ITER" + to_string(Ntop)+"_STEPS"+ to_string(STEPS) + "_RANGE" + to_string(range1) + "_TH_" + to_string(TH) +"_t" + tmg.timestamp() + ".csv" ); 


    /*
        Write basic information about the computation
    */
    cout << "PROGRAM START： " << tmg.time_pretty() << endl;
    cout << "Halo Problem Forward and Backward:" << endl;
    cout << "Reflection coefficients: " << to_string(RCOEFF) << endl;
    cout << "Neutrino potential mu: " << to_string(MIU) << endl;
    cout << "Average coefficient alpha: " << to_string(avg_alpha) << endl;
    cout << "Vacuum Omega: " << to_string(omegav) << endl;
    cout << "Vacuum Mixing angle sin 2theta: " << to_string(s2theta) << endl;
    cout << "Total number of iterations: " + to_string(Ntop) << endl;
    cout << "Size of rhos: " + to_string(SIZE) << endl;
    cout << "Range: from " + to_string(range0) + " to " + to_string(range1) << endl;
    cout << "Step size: " + to_string(STEPSIZE) << endl;
    cout << "Threads: " << to_string(TH) << endl;
    cout << "Save Steps: " << to_string(STEPSAVE) << endl;

    /* 
        Output warning if the step size and range and size of array are not consistant
    */
    if( (STEPS*STEPSIZE + range0) - range1 > 0.1 ){
        cout << "NOT REACHING THE END POINT!" << endl;
    }
    


    state_type rho_init = { 0.0, 0.0, 1.0 };

    /*
        The state arrays
    */
    StateArray* sa_store_ptr = new StateArray(SIZE);// = &sa_store;
    StateArray* sa_ptr = new StateArray(SIZE);// = &sa;
    double* z_array = new double[SIZE];

    /*
        Initalize all arrays to 0
    */
    double init_value[3] = {0.0,0.0,1.0};
    (*sa_ptr).init(init_value);
    (*sa_store_ptr).init(init_value);

    for(auto i=0; i < SIZE; i++)
    {
        z_array[i] = i * STEPSIZE;
    }
    
    Recorder::record_z(z_array, SIZE, data, STEPSAVEZ);

    /*
        Assign initial values to sa_store which is used to calculate sa
    */
    for(auto i = 0; i < 3; i++) {
        (*sa_store_ptr)[0][i] = rho_init[i]; 
        (*sa_ptr)[0][i] = rho_init[i]; 
    }

    // cout << "Initial condition: " << (*sa_store_ptr)[0][0] << ", " << (*sa_store_ptr)[0][1]  << ", "<< (*sa_store_ptr)[0][2] << endl;



    // Track clock time in openmp
    const clock_t begin_time = tmg.wall_stopwatch_reset();


    /*
        Main loop to find equilibrium
    */
    for(int i=0; i<Ntop; i++ ){

        // Average the two for damping, this method is slow!!!!
        // Helper::state_avg(sa_store_ptr, sa_ptr, SIZE, avg_alpha); // assign the average to sa_store_ptr

        // Evolve
        Looper::halo_euler_forward_one_avg(sa_ptr, sa_store_ptr, STEPSIZE, SIZE, avg_alpha, RCOEFF, MIU);
        
        // Save data and check termination condition
        if(saveflag%STEPSAVE == 0) {
            Recorder::record_state( sa_ptr, SIZE-1, data , STEPSAVEZ );
            Recorder::record_state( sa_ptr, SIZE-1, data , STEPSAVEZ, 1 ); // Recording the first element of state
            Recorder::record_state( sa_ptr, SIZE-1, data , STEPSAVEZ, 0 ); // Recording the second element of state
            // Recorder::cout_state( sa_store_ptr, SIZE-1, STEPSAVEZ );
            // Recorder::cout_state( sa_ptr, SIZE-1, STEPSAVEZ );
            if( (Helper::squared_difference(sa_store_ptr, sa_ptr, SIZE-1, STEPSAVEZ) < 1.0e-10) && (saveflag > 2) ){
                termflag++;
                if( termflag > 50 ){
                    cout << "Termination Condition Match: Terminating... " << endl;
                    break;
                }
            }
        }



        // swap the pointer to reuse the storage
        Helper::sa_ptr_swap(sa_ptr, sa_store_ptr);

        // Update saveflag to keep track of which step we are calculating
        ++saveflag;

    }


    /* 
        Calculate clock time for the iterations
    */
    tmg.wall_stopwatch_info(begin_time, Ntop);

    // Output data for benchmark
    // cout << to_string(Ntop) << ", " <<  to_string(SIZE) << ", " << argv[3] << ", " << float( clock () - begin_time ) /  CLOCKS_PER_SEC / Ntop << endl; // Ntop, SIZE, Threads, Time/Iteration

    cout << "At z=0, forward beam state: ";
    Recorder::cout_single_state( sa_store_ptr, 0);
    cout << "At z=0, backward beam state: ";
    Recorder::cout_single_state( sa_store_ptr, SIZE-2);
    
    cout << "At z="+ to_string(range1) +", forward beam state: ";
    Recorder::cout_single_state( sa_store_ptr, SIZE/2-1);
    cout << "At z="+ to_string(range1) +", backward beam state: ";
    Recorder::cout_single_state( sa_store_ptr, SIZE/2-1);

    cout << "PROGRAM END" << endl;

    return 0;
}
