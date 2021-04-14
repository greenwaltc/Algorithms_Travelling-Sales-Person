#!/usr/bin/python3

from which_pyqt import PYQT_VER
from queue import PriorityQueue
from collections import defaultdict
import heapq
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
        original_cities = [cities[city] for city in
                           range(len(cities))]  # Saves the cities for re-runs
        start_time = time.time()
        city_index = 0

        while not foundTour and time.time() - start_time < time_allowance:

            cities = [original_cities[city] for city in
                      range(len(original_cities))]
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

        '''Initialize BSSF to the greedy algorithm solution'''
        bssf = self.greedy(time_allowance)['soln']
        # Since the greedy algorithm is quick and provides a good solution, it
        # will provide really early pruning for B&B and improve its efficiency
        # overall. Less states are generated and thus less pruning is needed.

        '''Initialize state priority queue'''
        stateQueue = PriorityQueue()

        '''	Create the root of the state tree
		Reduce the cost matrix
		Set lower bound to cost of first reduction'''
        root = self.state()
        root.cost_matrix = [[-1 for i in range(self.ncities)] \
                            for k in range(self.ncities)]
        self.initializeState(None, root)
        root.city_num = 0  # Always start at first city
        root.path.append(root.city_num)
        lowerBound = root.cost

        stateQueue.put((root.cost, root))

        start_time = time.time()
        '''Begin the algorithm'''
        while not stateQueue.empty() \
                and time.time() - start_time < time_allowance:
            if stateQueue.qsize() > max_queue_size:
                max_queue_size = stateQueue.qsize()
            state = stateQueue.get()[1]
            if state.cost > bssf.cost:
                num_pruned += 1
                continue
            '''Make each child state'''
            for j in range(self.ncities):
                if time.time() - start_time > time_allowance:
                    break  # Over on time
                if state.cost_matrix[state.city_num][j] != math.inf:
                    # There is a path from this city to the next

                    '''Set up initial values for child'''
                    child = self.state()
                    self.initializeState(state, child, j)
                    num_states += 1

                    self.infRowCol(child)

                    '''Calculate State Cost'''
                    cost_reduction = self.reduceMatrix(child)
                    cost_step = child.parent.cost_matrix \
                            [child.parent.city_num][child.city_num]
                    cost_prev_state = child.parent.cost
                    child.cost = \
                        cost_prev_state + cost_step + cost_reduction

                    '''If the state is a leaf node and
                    it's less than BSSF so far, update
                    BSSF and continue to next state'''
                    if len(child.path) == self.ncities:
                        if child.cost < bssf.cost:
                            '''Make BSSF route'''
                            route = []
                            for i in range(self.ncities):
                                route.append(self.cities[child.path[i]])
                            bssf = TSPSolution(route)
                        num_solutions += 1
                        continue

                    '''Add child state to the queue'''
                    if bssf.cost > child.cost > lowerBound:
                        stateQueue.put(((child.cost / child.depth), child))
                        # Encourages digging deeper first
                    else:
                        num_pruned += 1

        end_time = time.time()
        results['cost'] = bssf.cost
        results['time'] = end_time - start_time
        results['count'] = num_solutions
        results['soln'] = bssf
        results['max'] = max_queue_size
        results['total'] = num_states
        results['pruned'] = num_pruned
        return results

    def initializeState(self, parent, child, j=0):
        if parent == None:
            # This is the root of the state tree, or state one
            root = child
            root.city_num = 0
            # first state assumes always starting at first city
            root.depth = 1
            '''Initialize first state cost matrix'''
            for i in range(self.ncities):
                for k in range(self.ncities):
                    root.cost_matrix[i][k] = self.cities[i].costTo(
                        self.cities[k])
            '''Reduce the cost matrix'''
            root.cost = self.reduceMatrix(root)
        else:
            # This is a child state
            child.parent = parent
            child.city_num = j
            child.depth = child.parent.depth + 1
            # Don't want the parent values to be overwritten so
            # make a deep copy
            child.cost_matrix = \
                copy.deepcopy(child.parent.cost_matrix)
            child.path = copy.deepcopy(child.parent.path)
            child.path.append(j)

    '''O(3n) ~= o(n)'''
    def infRowCol(self, state):
        '''Inf out appropriate row and column'''
        row = state.parent.city_num
        col = state.city_num
        for k in range(self.ncities):
            state.cost_matrix[row][k] = math.inf
        for k in range(self.ncities):
            state.cost_matrix[k][col] = math.inf

        '''Prevent premature cycles'''
        path_len = len(state.path)
        index = path_len - 1
        while index >= 0:
            row = state.city_num
            col = state.path[index]
            state.cost_matrix[row][col] = math.inf
            index -= 1

    def reduceMatrix(self, state):
        total_cost = 0
        '''Reduce row-by-row'''
        '''O(n^2)'''
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

    class state:
        def __init__(self):
            self.parent = None
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

    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

    def fancy(self, time_allowance=60.0):
        BSSF = self.greedy(time_allowance)['soln']
        results = {}
        run = 1
        max_runs = 8
        start_time = time.time()

        while run <= max_runs and time.time() - start_time < time_allowance:
            temp_sol = \
                self.fancy_helper(time_allowance, BSSF, start_time)
            if temp_sol.cost < BSSF.cost:
                BSSF = temp_sol
            run += 1  # This loop runs max_runs number of times

        end_time = time.time()
        results['cost'] = BSSF.cost
        results['time'] = end_time - start_time
        results['count'] = None
        results['soln'] = BSSF
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    '''
    Christofides Algorithm for solving the TSP. Runs in O(n^4)
    '''
    def fancy_helper(self, time_allowance=60.0, BSSF=None, start_time=time.time()):
        improvement_found = True
        BSSF_route = copy.deepcopy(BSSF.route)
        ncities = len(BSSF_route)
        results = {}

        while improvement_found and time.time() - start_time < time_allowance:
            improvement_found = False  # Start assumption that we have best sol
            # For each edge 1
            for i in range(0, ncities):

                # Gets us closer to the time allotment
                if time.time() - start_time > time_allowance:
                    break

                edge1 = (BSSF_route[i], BSSF_route[(i+1) % ncities])

                # All other possible edges
                for j in range(0, ncities):

                    if j == i or j == i - 1 or j == i + 1:
                        continue  # An edge can't be replaced with itself

                    # Gets us closer to the time allotment
                    if time.time() - start_time > time_allowance:
                        break

                    edge2 = (BSSF_route[j], BSSF_route[(j+1) % ncities])

                    # Begin swapping edges
                    # For each possible permutation of edges, try making a path
                    # and if the new path's cost is less than the previous bssf
                    # then update the bssf and note that we've found an
                    # improvement
                    for p in range(1, 23):

                        # There are 24 permutations possible. The first is
                        # always the original edges, so skip 0

                        # Gets us closer to the time allotment
                        if time.time() - start_time > time_allowance:
                            break

                        new_edge1, new_edge2 = \
                            self.permutation_edge(p, edge1, edge2)

                        # Replace the cities in the previous bssf path with
                        # the new cities
                        new_route = copy.deepcopy(BSSF_route)
                        new_route[i] = new_edge1[0]
                        new_route[(i+1) % ncities] = new_edge1[1]
                        new_route[j] = new_edge2[0]
                        new_route[(j+1) % ncities] = new_edge2[1]

                        # Make the solution and compare cost. If the new route
                        # is invalid the cost will be infinity.
                        solution = TSPSolution(new_route)
                        if solution.cost < BSSF.cost:
                            # Update the BSSF
                            improvement_found = True
                            BSSF = solution
                            print("Improvement found!")
                    if improvement_found:
                        break
                if improvement_found:
                    break
            if improvement_found:
                break
        return BSSF

    '''
    Reorders the original edge according to the permutation number
    
    Possible permutations:
    {A,B,C,D} {A,B,D,C} {A,C,B,D} {A,C,D,B} {A,D,B,C} {A,D,C,B} 
    {B,A,C,D} {B,A,D,C} {B,C,A,D} {B,C,D,A} {B,D,A,C} {B,D,C,A} 
    {C,A,B,D} {C,A,D,B} {C,B,A,D} {C,B,D,A} {C,D,A,B} {C,D,B,A} 
    {D,A,B,C} {D,A,C,B} {D,B,A,C} {D,B,C,A} {D,C,A,B} {D,C,B,A}
    
    where A,B,C,D are cities and the original edge is in the form:
        edge1: (A,B); edge2: (C,D)
    '''
    def permutation_edge(self, p, edge1, edge2):
        A = 0; B = 1; C = 0; D = 1
        if p == 0:
            return None, None # This if statement should never be reached
        if p == 1:
            return (edge1[A], edge1[B]), (edge2[D], edge2[C])
        if p == 2:
            return (edge1[A], edge2[C]), (edge1[B], edge2[D])
        if p == 3:
            return (edge1[A], edge2[C]), (edge2[D], edge1[B])
        if p == 4:
            return (edge1[A], edge2[D]), (edge1[B], edge2[C])
        if p == 5:
            return (edge1[A], edge2[D]), (edge2[C], edge1[B])
        if p == 6:
            return (edge1[B], edge1[A]), (edge2[C], edge2[D])
        if p == 7:
            return (edge1[B], edge1[A]), (edge2[D], edge2[C])
        if p == 8:
            return (edge1[B], edge2[C]), (edge1[A], edge2[D])
        if p == 9:
            return (edge1[B], edge2[C]), (edge2[D], edge1[A])
        if p == 10:
            return (edge1[B], edge2[D]), (edge1[A], edge2[C])
        if p == 11:
            return (edge1[B], edge2[D]), (edge2[C], edge1[A])
        if p == 12:
            return (edge2[C], edge1[A]), (edge1[B], edge2[D])
        if p == 13:
            return (edge2[C], edge1[A]), (edge2[D], edge1[B])
        if p == 14:
            return (edge2[C], edge1[B]), (edge1[A], edge2[D])
        if p == 15:
            return (edge2[C], edge1[B]), (edge2[D], edge1[A])
        if p == 16:
            return (edge2[C], edge2[D]), (edge1[A], edge1[B])
        if p == 17:
            return (edge2[C], edge2[D]), (edge1[B], edge1[A])
        if p == 18:
            return (edge2[D], edge1[A]), (edge1[B], edge2[C])
        if p == 19:
            return (edge2[D], edge1[A]), (edge2[C], edge1[B])
        if p == 20:
            return (edge2[D], edge1[B]), (edge1[A], edge2[C])
        if p == 21:
            return (edge2[D], edge1[B]), (edge2[C], edge1[A])
        if p == 22:
            return (edge2[D], edge2[C]), (edge1[A], edge1[B])
        if p == 23:
            return (edge2[D], edge2[C]), (edge1[B], edge1[A])
