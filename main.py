#The problem is to optimize the delta V, time, or distance needed to travel
#between x satellites in given orbits based on their orbit parameters ri & vi.
#This code defines the funtions necessary to input initial conditions and 
#read results of the QPU.

#============================================

from dwave.system import LeapHybridSampler
import pandas as pd
from parambulator import propellant_calc
import networkx as nx
import dwave_networkx as dnx


Data_name = "Orbit_parameters2.txt" #text file with all orbit parameters for satellites 
Mo = 1000 #Total mass of space craft in kg
Isp = 1660 # Specific impulse of engine in s

#pars input data
params = pd.read_csv(Data_name, header=None)
params.columns = ['r0x','r0y','r0z','v0x','v0y','v0z']

#Calculate fuel requirements to move from one satellite to another
deltas = propellant_calc(params['r0x'],params['r0y'],params['r0z'],\
                      params['v0x'],params['v0y'],params['v0z'],Mo,Isp)   
    
#Turn matrix of Deltas into nodes and edges for reading
#Graph is directed so traveling from u to v does not have same weight
#as traveling from v to u                         
routes = nx.DiGraph()
for i in range(len(deltas)):
    for j in range(i+1,len(deltas)):
        routes.add_weighted_edges_from([(i,j,deltas[i][j]), (j,i,deltas[j][i])])

#choose sampler
sampler = LeapHybridSampler()

#Test routes on quantum chip
bqm = dnx.traveling_salesperson(routes,sampler)

#Verify path is valid (all locations visited once) and result is in the list of routes
if (set(bqm) == set(routes)):
  print('Best route: ' +str(bqm))
else:
  print('Valid path not found')
