# Quantum-computer code

This code is meant to be used on the cloud based Dwave Leap IDE to sample a quantum processor

Inputs are read into the main program through a text file with format r0x,r0y,r0z,v0x,v0y,v0z. Each line corresponds to a different satellite.

The inputs are passed to the varialbe_deltas function in parambulator to calculate the time, delta V, or propellant needed to travel from one satellite to another. The values are used as weights for the edges between our points (satellites)

Finally the nodes and edges, with weights, are added to the routes graph and passed to the Leap Hybrid Sampler to be run on the QPU.

As of June 17th 2020 the predefined traveling_salesperson function is used to find the optimal route.
