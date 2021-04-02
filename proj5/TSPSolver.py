#!/usr/bin/python3

from which_pyqt import PYQT_VER
import random

if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools


class TSPSolver:
	def __init__(self, gui_view):
		self._scenario = None

	def setupWithScenario(self, scenario):
		self._scenario = scenario

	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def defaultRandomTour(self, time_allowance=60.0):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time() - start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation(ncities)
			route = []
			# Now build the route using the random permutation
			for i in range(ncities):
				route.append(cities[perm[i]])
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def greedy(self, time_allowance=60.0):
		random.seed(time.time())
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0  # How many tries it took to find a path
		bssf = None
		original_cities = [cities[city] for city in range(len(cities))] # Saves the cities for re-runs
		start_time = time.time()
		city_index = 0

		while not foundTour and time.time() - start_time < time_allowance:

			cities = [original_cities[city] for city in range(len(original_cities))]
			route = []
			# current_city = cities.pop(random.randint(0, len(cities) - 1))
			current_city = cities.pop(city_index)
			city_index += 1
			origin_city = current_city
			route.append(origin_city)

			while len(cities) > 0 and time.time() - start_time < time_allowance:

				min_distance = np.inf
				min_city = None

				# Find the minimum city
				index = -1
				for i in range(len(cities)):
					assert (type(cities[i]) == City)
					if current_city.costTo(cities[i]) < min_distance:
						min_distance = current_city.costTo(cities[i])
						min_city = cities[i]
						index = i
				cities.pop(index)  # Remove the min city from list of cities

				if min_city is not None:
					current_city = min_city

					if len(cities) == 0:
						# Check if there's a path back to the origin city
						if min_city.costTo(origin_city) < float('inf'):
							route.append(current_city)
							foundTour = True
							break
						else:
							foundTour = False
							count += 1
							break

					route.append(current_city)
				else:
					# There's no solution, rerun
					count += 1
					foundTour = False
					break

		bssf = TSPSolution(route)
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''

	def branchAndBound(self, time_allowance=60.0):
		pass

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''

	def fancy(self, time_allowance=60.0):
		pass
