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
        '''
            n:param the number of nodes in the graph
        '''

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

    '''
    Christofides Algorithm for solving the TSP. Runs in O(n^4)
    '''
    def fancy(self, time_allowance=60.0):
        random.seed(time.time())
        results = {}
        self.cities = self._scenario.getCities()
        self.ncities = len(self.cities)

        # Create the minimum spanning tree
        mst = self.MST(0)
        # mst has a format of {{from_city: {<to_cities>}}}
        # Or, a set of dictionaries to sets

        # Initialize all the cities into Cities
        even_odd_cities = [self.City(self.cities[i])
                           for i in range(self.ncities)]

        # Find the odd edges
        for group in mst:
            # flip group
            to_set = mst[group]
            for to in to_set:
                groupEven = even_odd_cities[group].isEven
                even_odd_cities[group].isEven = False if groupEven else True
                toEven = even_odd_cities[to].isEven
                even_odd_cities[to].isEven = False if toEven else True
        odd_cities = [even_odd_cities[i] for i in range(self.ncities)
                      if even_odd_cities[i].isEven is False]

        pass

    '''
    An implementation of Prim's algorithm to find the minimum spanning tree of
    the graph of cities.
    Algorithm steps found here: https://bradfieldcs.com/algos/graphs/prims-spanning-tree-algorithm/
    
    :returns: forest of cities/edges F
    '''
    def MST(self, starting_vertex):
        mst = defaultdict(set)
        edgeExists = self._scenario.getEdgeExists()
        visited = set([starting_vertex])
        possible_cities = [j for j in range(self.ncities)
                           if edgeExists[starting_vertex][j]]
        edges = [
            self.Edge(self.cities[starting_vertex],
                      self.cities[possible_cities[to]])
            for to in range(len(possible_cities))
        ]

        heapq.heapify(edges)

        while edges:
            edge = heapq.heappop(edges)
            frm = edge.origin_city._index
            to = edge.destination_city._index

            if to not in visited:
                visited.add(to)
                mst[frm].add(to)
                possible_cities = [j for j in range(self.ncities)
                                   if edgeExists[to][j]]
                for to_next in possible_cities:
                    if to_next not in visited:
                        heapq.heappush(edges, self.Edge(self.cities[to],
                                                        self.cities[to_next]))
        return mst

    class Edge:
        def __init__(self, origin, destination):
            assert(type(origin) == City and
                   type(destination) == City)
            self.origin_city = origin
            self.destination_city = destination
            self.cost = origin.costTo(destination)

        def __gt__(self, other):
            if (self.cost > other.cost):
                return True
            else:
                return False

    class City:
        def __init__(self, city):
            self.city = city
            self.isEven = True



