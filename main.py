#The problem is to optimize the delta V, time, or distance needed to travel
#between x satellites in given orbits based on their orbit parameters ri & vi.
#This code defines the funtions necessary to input initial conditions and 
#read results of the QPU.

#============================================

from dwave.system import LeapHybridSampler
import pandas as pd
from parambulator import variable_deltas
import networkx as nx
import dwave_networkx as dnx


Data_name = "Orbit_parameters.txt" #text file with all orbit parameters for satellites 
Mp = 300 #mass of propellant in kg
Mo = 500 #Total mass of space craft in kg
Isp = 300 # Specific impulse of engine in s

#variable to return in deltas
optimize = 't'

#pars input data
params = pd.read_csv(Data_name, header=None)
params.columns = ['r0x','r0y','r0z','v0x','v0y','v0z']

#Calculate delta V's, fuel requirements, and delta t's to move from one orbit to another
deltas = variable_deltas(params['r0x'],params['r0y'],params['r0z'],\
                      params['v0x'],params['v0y'],params['v0z'],\
                          Mp,Mo,Isp,optimize)

#create maximum value of worst route for weighting each edge
total_max = 0
if optimize == 'm':
    total_max = Mp
else:
    for i in range(len(deltas)-1):
        total_max += max(deltas[i][i+1:len(deltas)])    
    
#Turn matrix of Deltas into nodes and edges for reading
#Graph is undirected so traveling from u to v has same weight
#as traveling from v to u                         
routes = nx.Graph()
for i in range(len(deltas)):
    for j in range(i+1,len(deltas)):
        routes.add_weighted_edges_from({(i,j,(deltas[i][j])/total_max), (j,i,deltas[j][i]/total_max)})

#choose sampler
sampler = LeapHybridSampler()

#Test routes on quantum chip
bqm = dnx.traveling_salesperson(routes,sampler)
