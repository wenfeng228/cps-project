#!/bin/bash
#SBATCH --job-name=sensitivity_lambda
#SBATCH -C intel
#SBATCH --output=sensitivity_lambda.txt
#SBATCH --error=sensitivity_lambda_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00

exe=../code/simulator_omp
chmod +x $exe

echo "model,graph,factor,time1,time2,total_time,final_step,final_l2,gamma,lambda,budget,max_step,tolerance,thread_count,num_nodes,num_edges" > result/sensitivity_lambda.csv

for lambda in 0.1 0.2 0.3 0.4; do
    echo "Running: lambda=$lambda"
    export OMP_NUM_THREADS=4
    srun -n 1 -c 4 $exe 0.1 $lambda 100 100000 1 "../data/nodes_6.csv" "../data/edges_6.csv" "../output" 10 1000 >> result/sensitivity_lambda.csv
    mv ../output/adjacency_list_uniform_trace_6.csv ../output/mod1_lambda${lambda}_6.csv
    mv ../output/adjacency_list_local_spillover_trace_6.csv ../output/mod2_lambda${lambda}_6.csv
    mv ../output/adjacency_list_global_spillover_trace_6.csv ../output/mod3_lambda${lambda}_6.csv
    mv ../output/compressed_sparse_row_uniform_trace_6.csv ../output/mod4_lambda${lambda}_6.csv
    mv ../output/compressed_sparse_row_local_spillover_trace_6.csv ../output/mod5_lambda${lambda}_6.csv
    mv ../output/compressed_sparse_row_global_spillover_trace_6.csv ../output/mod6_lambda${lambda}_6.csv
done

echo "wrote result/sensitivity_lambda.csv"
