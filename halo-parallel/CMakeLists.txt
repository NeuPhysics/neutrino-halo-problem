cmake_minimum_required(VERSION 3.7)
project(halo_parallel)

set(CMAKE_CXX_STANDARD 11)

# set(SOURCE_FILES test.cpp include/state.h include/looper.h include/hamiltonian.h include/stepper.h include/recorder.h include/loader.h)
# add_executable(test-suite ${SOURCE_FILES})
#set(SOURCE_FILES neutrino-headon.cpp include/state.h include/looper.h include/hamiltonian.h include/stepper.h include/recorder.h include/loader.h)
#add_executable(neutrino-headon ${SOURCE_FILES})
set(SOURCE_FILES neutrino-headon-nunubar.cpp include/state.h include/looper.h include/hamiltonian.h include/stepper.h include/recorder.h include/loader.h)
add_executable(neutrino-headon-nunubar ${SOURCE_FILES})

#set(SOURCE_FILES neutrino-headon-avg.cpp include/state.h include/looper.h include/hamiltonian.h include/stepper.h include/recorder.h include/loader.h)
#add_executable(neutrino-headon-avg ${SOURCE_FILES})

# set(SOURCE_FILES neutrino-headon-incline.cpp include/state.h include/looper.h include/hamiltonian.h include/stepper.h include/recorder.h include/loader.h)
# add_executable(neutrino-headon-incline ${SOURCE_FILES})
# set(SOURCE_FILES neutrino-headon-evo.cpp include/state.h include/looper.h include/hamiltonian.h include/stepper.h include/recorder.h include/loader.h)
# add_executable(neutrino-headon-evo ${SOURCE_FILES})
# for parallel using openmp
SET(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc/7.2.0/bin/g++)
# SET(CMAKE_CXX_COMPILER /usr/bin/g++)
SET(GCC_COVERAGE_COMPILE_FLAGS "-g -fopenmp -Wall")
# SET(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )