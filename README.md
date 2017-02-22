# Modelling-Markovian-Systems
Simulation of an Ising Ring Chain (Generalized Hepp Coleman) Model which is utilized as a controllably Markovian system. The model consists of a ring of spins surrounding a central two-state system (S). The central two-level system S possesses homogeneous couplings due to the overlaps of symmetric spacial wave function of S with those of spins. Refer to [[1]](https://arxiv.org/abs/quant-ph/0509007) for further details.

![ising ring](https://cloud.githubusercontent.com/assets/20701981/23200158/cda59c16-f8f8-11e6-996f-4aafdcc409ec.jpg)

#Program Listing

1. RGroundStateCalc: This file calculates the ground state of the hamiltonian and is saved into a file "Rgstate".

2. RDiagEvolOperator: In this file first the energy eigenvalues of the evolution operator of the excited state are stored into a variable and then a diagonal matrix is formed using this variable. Similarly for the ground state evolution operator and hence both the unitary matrices which dictate the evolution of the system from the intial state "Rgstate" is obtained. For varying vaues of lambda different unitary evolution matrices are calculated and the states evolved according to them and Loschmidt Echo calculated. 

![3d actual](https://cloud.githubusercontent.com/assets/20701981/23200658/8abbff6e-f8fb-11e6-8fe0-c461d67818ed.png)

A distinct minima is observed when varying the interaction parameter. 
