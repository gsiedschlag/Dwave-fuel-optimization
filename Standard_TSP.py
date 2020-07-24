import numpy as np
import time
from dwave.system import LeapHybridSampler
import networkx as nx
import dwave_networkx as dnx
import pandas as pd

def create_cities(N):
    """
    Creates an array of random points of size N.
    """
    cities = []
    for i in range(N):
        cities.append(np.random.rand(2) * 10)
    return np.array(cities)

def distance_between_points(point_A, point_B):
    return np.sqrt((point_A[0] - point_B[0])**2 + (point_A[1] - point_B[1])**2)

def get_distance_matrix(cities):
    number_of_cities = len(cities)
    matrix = np.zeros((number_of_cities, number_of_cities))
    for i in range(number_of_cities):
        for j in range(i, number_of_cities):
            matrix[i][j] = distance_between_points(cities[i], cities[j])
            matrix[j][i] = matrix[i][j]
    return matrix

#single run code
#cities = create_cities(int(input('Enter number of nodes to use: ')))

Data_name = "20_TSP_city_locations.txt" #text file with all city locations 

params = pd.read_csv(Data_name, header=None)
params.columns = ['x','y']
nodes = int(input('Enter number of nodes to use: '))
cities = params[:nodes][:]

cost_matrix = get_distance_matrix(cities)

routes = nx.Graph()
for i in range(len(cities)):
    for j in range(i+1,len(cities)):
        routes.add_weighted_edges_from({(i,j,cost_matrix[i][j]), (j,i,cost_matrix[j][i])})

#choose sampler
sampler = LeapHybridSampler()

start_time = time.time()

#Test routes on quantum chip
bqm = dnx.traveling_salesperson(routes,sampler)

#Verify path is valid (all locations visited once) and result is in the list of routes
if (set(bqm) == set(routes)):
  print('Solution took ', round(time.time() - start_time, 3), 's for cities = ', str(len(cities)))
  print('Best route: ' +str(bqm))
else:
  print('Valid path not found')
