#!/bin/bash
#SBATCH --job-name=complexity_scaling
#SBATCH -C intel
#SBATCH --output=complexity_scaling.txt
#SBATCH --error=complexity_scaling_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00

exe=../code/simulator_omp
chmod +x $exe

echo "model,graph,factor,time1,time2,total_time,final_step,final_l2,gamma,lambda,budget,max_step,tolerance,thread_count,num_nodes,num_edges" > result/complexity_scaling.csv

for data in 0 1 2 3 4 5 6 7; do
    echo "Running: data=$data"
    export OMP_NUM_THREADS=4
    for i in {1..3}; do
        echo "Running: iteration=$i"
        srun -n 1 -c 4 $exe 0.1 0.25 100 100000 0 "../data/nodes_${data}.csv" "../data/edges_${data}.csv" "../output" 10 1000 >> result/complexity_scaling.csv
    done
done

echo "wrote result/complexity_scaling.csv"
