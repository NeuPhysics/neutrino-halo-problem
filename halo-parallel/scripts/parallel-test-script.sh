echo "Ntop, Size, OMPThread, Time" > halo_sim_parallel_benchmark.record;

for NTOP in 100 1000 2000 3000
do
   for ITER in 1000 10000 20000 30000 40000
   do
      for OMPTH in 1 2 3 4 5 6 7 8 9 10 11
      do
         echo "$NTOP, $ITER, $OMPTH"
         nohup ./gxx/euler-parallel-openmp-supernova.o $NTOP $ITER $OMPTH >> halo_sim_parallel_benchmark.record;
      done
   done
done
