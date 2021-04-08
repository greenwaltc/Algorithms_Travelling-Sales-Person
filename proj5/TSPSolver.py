#!/usr/bin/python3

from which_pyqt import PYQT_VER
from queue import PriorityQueue
import numpy as np
import copy
import sys
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
        original_cities = [cities[city] for city in range(len(cities))]  # Saves the cities for re-runs
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

        max_queue_size = 0
        num_states = 0
        num_pruned = 0
        num_solutions = 0
        random.seed(time.time())
        results = {}
        self.cities = self._scenario.getCities()
        self.ncities = len(self.cities)
        self.lowerBound = 0
        self.state_num = 0

        '''Initialize BSSF to the greedy algorithm solution'''
        self.bssf = self.greedy(time_allowance)['soln']

        '''Initialize state priority queue'''
        stateQueue = PriorityQueue()

        '''	Create the root of the state tree
			Reduce the cost matrix
			Set lower bound to cost of first reduction'''
        root = self.state(self.ncities)
        root.cost_matrix = [[-1 for i in range(self.ncities)] for k in range(self.ncities)]
        self.initializeState(None, root)
        root.city_num = 0  # Always start at first city
        root.path.append(root.city_num)
        self.lowerBound = root.cost

        stateQueue.put((root.cost / root.depth, root))  # Dividing by depth encourages digging deeper first
        start_time = time.time()
        '''Begin the algorithm'''
        while not stateQueue.empty() and time.time() - start_time < time_allowance:
            if stateQueue.qsize() > max_queue_size:
                max_queue_size = stateQueue.qsize()
            state = stateQueue.get()[1]
            if state.cost < self.bssf.cost:
                '''Make each child state'''
                for j in range(self.ncities):
                    if time.time() - start_time > time_allowance:
                        break  # Over on time
                    if state.cost_matrix[state.city_num][j] != math.inf:
                        # There is a path from this city to the next
                        '''Set up initial values for child'''
                        child = self.state(self.ncities)
                        num_states += 1
                        state.children.append(child)
                        child.parent = state
                        child.state_num = self.state_num
                        self.state_num += 1
                        child.city_num = j
                        child.depth = child.parent.depth + 1
                        child.cost_matrix = copy.deepcopy(
                            child.parent.cost_matrix)  # Don't want the parent's matrix values to change
                        child.path = copy.deepcopy(child.parent.path)
                        child.path.append(j)

                        '''Inf out appropriate row and column'''
                        row = child.parent.city_num
                        col = child.city_num
                        for k in range(self.ncities):
                            child.cost_matrix[row][k] = math.inf
                        for k in range(self.ncities):
                            child.cost_matrix[k][col] = math.inf

                        '''Prevent premature cycles'''
                        path_len = len(child.path)
                        index = path_len - 1
                        while index >= 0:
                            child.cost_matrix[child.city_num][child.path[index]] = math.inf
                            index -= 1

                        '''Calculate State Cost'''
                        cost_reduction = self.reduceMatrix(child)
                        cost_step = child.parent.cost_matrix[child.parent.city_num][child.city_num]
                        cost_prev_state = child.parent.cost
                        child.cost = cost_prev_state + cost_step + cost_reduction

                        '''If the state is a leaf node and 
						it's less than BSSF so far, update
						BSSF and continue to next state'''
                        if len(child.path) == self.ncities:
                            if child.cost < self.bssf.cost:
                                route = []
                                for i in range(self.ncities):
                                    route.append(self.cities[child.path[i]])
                                self.bssf = TSPSolution(route)
                            num_solutions += 1
                            continue

                        '''Add child state to the queue'''
                        if self.bssf.cost > child.cost > self.lowerBound:
                            stateQueue.put((child.cost / child.depth, child))
                        else:
                            if stateQueue.qsize() > max_queue_size:
                                max_queue_size = stateQueue.qsize()
                            num_pruned += 1
            else:  # Found the local optimal solution
                if stateQueue.qsize() > max_queue_size:
                    max_queue_size = stateQueue.qsize()
                num_pruned += stateQueue.qsize()
                break
        end_time = time.time()
        results['cost'] = self.bssf.cost  # if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = num_solutions
        results['soln'] = self.bssf
        results['max'] = max_queue_size
        results['total'] = num_states
        results['pruned'] = num_pruned
        return results

    def initializeState(self, parent, child):
        # This is the root of the state tree, or state one
        if parent == None:
            child.city_num = 0  # first state assumes always starting at first city
            child.state_num = self.state_num
            self.state_num += 1
            child.depth = 1
            '''Initialize first state cost matrix'''
            for i in range(self.ncities):
                for k in range(self.ncities):
                    child.cost_matrix[i][k] = self.cities[i].costTo(self.cities[k])
            '''Reduce the cost matrix'''
            child.cost = self.reduceMatrix(child)

        pass

    def reduceMatrix(self, state):
        total_cost = 0
        '''Reduce row-by-row'''
        for i in range(self.ncities):
            row_min = math.inf
            '''Find the minimum value in the row'''
            for j in range(self.ncities):
                if state.cost_matrix[i][j] < row_min:
                    row_min = state.cost_matrix[i][j]

            '''Subtract minimum value from each position in row'''
            if row_min != math.inf and row_min != 0:
                total_cost += row_min
                for j in range(self.ncities):
                    state.cost_matrix[i][j] -= row_min

        '''Reduce column by column'''
        for j in range(self.ncities):
            col_min = math.inf
            '''Find the minimum value in the column'''
            for i in range(self.ncities):
                if state.cost_matrix[i][j] < col_min:
                    col_min = state.cost_matrix[i][j]

            '''Subtract minimum value from each position in column'''
            if col_min != math.inf and col_min != 0:
                total_cost += col_min
                for i in range(self.ncities):
                    state.cost_matrix[i][j] -= col_min

        return total_cost

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

    class state:
        '''
			n:param the number of nodes in the graph
		'''

        def __init__(self, n):
            self.state_num = -1
            self.parent = None
            self.children = []
            self.cost = -1
            self.city_num = None
            self.depth = 0
            self.cost_matrix = [[None], [None]]
            self.path = []  # Stores the path so far in terms of city indexes

        def __gt__(self, other):
            if (self.cost > other.cost):
                return True
            else:
                return False
