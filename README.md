# GW-to-identify-subnetworks

This is for the paper: Sun, Yahui, et al. "A physarum-inspired prize-collecting steiner tree approach to identify subnetworks for drug repositioning." BMC systems biology 10.5 (2016): 128.

This repository contains the MATLAB coding of the GW algorithm, which is used to identify subnetworks (solve Prize-Collecting Steiner Tree Problem).

There are seven MATLAB files in total: GW.m is the main code of GW algorithm; 
D_01_b.mat is a drug similarity network; Function_merge_cpn.m is a function in the GW growth phase to merge two components; 
Function_OutputDegree.m is a function to calculate the degrees of vertices; 
Function_StrongPrune2.m is a function of the strong pruning algorithm; 
Function_pdegree.m is a function used in the strong pruning process to calculate the prune degree; 
Function_solution is a function to calculate the quality of the identified subnetwork.

To run GW algorithm, put all the MATLAB files in the same document, and run GW.m.

For any issue, please feel free to email me, syhhit@gmail.com -Yahui
