This repository contains a snapshot of the codes written in my thesis.

1- POST.m file is the code for the power functions of the Point-Optimal Sign-Based tests for predictive regressions
under the assumption that the errors are independently distributed.

2- PCCPOS_HPC.m is the code for the power functions of the PCC-POS tests, intended to be calculated with a high performance computing cluster. 

 2(a)- Par.slurm is to be accompanied with PCCPOS_HPC.m.

 2(b)- wrapperPar.m is to be accompanied with Par.slurm and PCCPOS_HPC.m.

 2(c)- kendalltau.m is to be accompanied with PCCPOS_HPC.m.	

  E.g. SLURM command: sbatch -N 1 -p par6.q --exclusive -a 1-12 par.slurm

3- BonferroniCausalityTest.R is the code for the power function of the bound-type procedure for testing the null-hypothesis of Granger non-causality.

4- SimAnneal_POS.m is the code for the simulated annealing algorithm for the porjection technique of the  POS-based tests  

 4(a)- wrapperParAnnealing.m is the wrapper for SimAnneal_POS.m

5- ProjectionTechnique_GridSearch.m is the code for the projection technique using the grid search algorithm.

 5(a)- cartprod.m is to be used with ProjectionTechnique_GridSearch.m
 5(b)- ind2subVect.m is to be used with ProjectionTechnique_GridSearch.m	
