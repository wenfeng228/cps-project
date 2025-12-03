#!/bin/bash
#SBATCH --job-name=weak_scaling
#SBATCH -C intel
#SBATCH --output=weak_scaling.txt
#SBATCH --error=weak_scaling_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00

exe=../code/simulator_omp
chmod +x $exe

echo "model,graph,factor,time1,time2,total_time,final_step,final_l2,gamma,lambda,budget,max_step,tolerance,thread_count,num_nodes,num_edges" > result/weak_scaling.csv

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

for exp in {0..4}; do
    thread=$((2**exp))
    data=$((exp+3))
    export OMP_NUM_THREADS=$thread
    echo "Running: Threads=$thread, Data=*_${data}.csv"
    for i in {1..3}; do
        echo "Running: iteration=$i"
        srun -n 1 -c $thread $exe 0.1 0.25 100 100000 0 "../data/nodes_${data}.csv" "../data/edges_${data}.csv" "../output" 10 1000 >> result/weak_scaling.csv
    done
done

echo "done. wrote result/weak_scaling.csv"
