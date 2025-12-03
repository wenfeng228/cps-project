#!/bin/bash
#SBATCH --job-name=strong_scaling
#SBATCH -C intel
#SBATCH --output=strong_scaling.txt
#SBATCH --error=strong_scaling_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00

exe=../code/simulator_omp
chmod +x $exe

echo "model,graph,factor,time1,time2,total_time,final_step,final_l2,gamma,lambda,budget,max_step,tolerance,thread_count,num_nodes,num_edges" > result/strong_scaling.csv

for thread in 1 2 4 6 8 10 12 14 16; do
    echo "Running: Threads=$thread"
    export OMP_NUM_THREADS=$thread
    for i in {1..3}; do
        echo "Running: iteration=$i"
        srun -n 1 -c $thread $exe 0.1 0.25 100 100000 0 "../data/nodes_6.csv" "../data/edges_6.csv" "../output" 10 1000 >> result/strong_scaling.csv
    done
done

echo "wrote result/strong_scaling.csv"
