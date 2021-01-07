# GlobOptBernstein
Constrained global optimization using parallelized Bernstein algorithm.

The paper is avaliable at https://arxiv.org/abs/2003.01758

## Requirements
To run the files in this repository, you will need the following:
- CUDA in the correct version compatible with your MATLAB version. 
For more information, please refer to https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html;jsessionid=165811ee071d12a1e7413420a560
- MATLAB Parallel Computing Toolbox

## Usage
- run run_benchmark_evaluation.m and run_benchmark_evaluation_P7.m to generate and save results for section V.B Benchmark Evaluation.
- run Benchmark_Evaluation/analyze_benchmark_results.m to generate plots of results for section V.B Benchmark Evaluation.
- run run_increasing_constraints_problems.m for section section V.C to generate and save results in section V.C Increasing Constraint Problems.
- run increasing_constraints/analyze_increasing_constraints_results.m to generate plots of results for section V.C Increasing Constraint Problems.
- run hardware_demo/analyze_PCBA_hardware_demo_results.m to generate plots of results for section VI.
- if you want the raw hardware data, please send an email to Shreyas Kousik

## Authors
Bohao Zhang

Shreyas Kousik

ROAHM LAB
University of Michigan

Visit our website: http://roahmlab.com
