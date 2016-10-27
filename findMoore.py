import numpy as np
import itertools as it
import igraph as ig

import argparse
import sys
import pickle
import datetime

from math import sqrt, floor

class Initialiser:
	def __init__(self, new_val, accelerate, continue_file, save_option):
		self.new_val = new_val
		self.accelerate = accelerate
		self.continue_file = continue_file
		self.save_option = save_option

	def initialise(self):
		global wb
		# set up saver
		global saver
		saver = Saver(self.save_option)
		
		# from user input, decide whether to start a new search or continue an old one
		if self.new_val:
			# if starting a new search then 
			# instantiate a new WorkingBasket with proper valence
			wb = WorkingBasket(self.new_val)
			# have the WorkingBaket run it's new_basket method
			wb.new_basket()
			# Optional: run an accelerator to get the known-present edges in without running through the huristic checks
			if self.accelerate:
				accelerate() # todo

		elif self.continue_file:
			# if continuing an old one then
			# unpickle the matrix
			wb = pickle.load(open(self.continue_file, "rb"))
			# override user's save option (helping out forgetful users)
			#  in the event of a continued search, the user will likely want to save progress.
			saver.save_option = True # you're welcome self.
			# Optional: check to make sure the loaded subgraph is well formatted, otherwise trust it # todo

		else:
			sys.exit("[!] You broke it. See --help for more details.")
		# create a manager and pass off control to the manager
		global manager
		manager = Manager()
		manager.search()

class WorkingBasket:
	def __init__(self, valence):
		self.val = valence # valence of target Moore graph
		self.v = self.val * (self.val - 1) # number of vertices in the basket
		self.e = self.v * self.val / 2 # number of edges in the target basket
		self.fruit = tuple(tuple(i*(self.val-1)+j for j in range(0,self.val-1)) for i in range(0,self.val))

	def new_basket(self):
		self.matrix = np.zeros(shape=(self.v, self.v), dtype=np.int32) # working basket zeroed out
		# rule out intra-fruit edges and self-loops
		intra_fruit = it.chain.from_iterable([it.combinations(self.fruit[f],2) for f in range(self.val)])
		[self.label(edge, -(self.e + 1)) for edge in intra_fruit]
		[self.label(edge, -(self.e + 1)) for edge in [(i,i) for i in range(self.v)]]
		self.current_edge_label = 0 # label number of the most recent edge placed

	def label(self, edge, label):
		# changes the label of an edge in self.matrix
		self.matrix[edge[0],edge[1]] = label
		self.matrix[edge[1],edge[0]] = label

	def flat_matrix(self):
		# return a uint8 np array that 'flattens' self.matrix
		# results in a much smaller matrix (1/4th the size) for use in parallelisation of
		#  heuristics. Is expensive to compute so use wisely.
		return np.vectorize(self._flatten, otypes=[np.bool])(self.matrix)

	def degree_of(self, vertex):
		# computes the degree of vertex in the current working basket
		return (self.matrix[vertex] > 0).sum()

	def pot_degree_of(self, vertex):
		# computes the number of potential edges at vertex
		return (self.matrix[vertex] == 0).sum()

	def _flatten(self, value):
		# !critical component of feasibility_check_2() -- do not alter
		if value < 0:
			return False
		else:
			return True

class Manager:
	def __init__(self):
		global wb
		global saver
		global heuristic_conductor
		heuristic_conductor = HeuristicConductor()

	def search(self):
		while True:
			print 'edge_label is ' + str(wb.current_edge_label)
			# EdgePicker tries to pick an available edge to label WorkingBasket.current_edge_label
			try:
				new_edge = choose_new_edge('deep')
			except EdgesExhausted: 
				# if there are no more available edges to pick, then we're finished and we need to check
				#  whether we've found a Moore graph or whether one isn't possible
				endgame()
			else:
				# edges are available to place. But *first* we need to check whether the current subgraph is feasible
				# !Checking this first is critical to the proper operation of the algorithm. Do not place a new edge before 
				#  running the heuristics first!

				# HeuristicChecker runs heuristic checks to see if current basket can be ruled out as a subgraph
				if not heuristic_conductor.passes():
					# if the most recently placed edge can be ruled out then:
					# edge placement number WorkingBasket.current_edge_label is marked and reverted
					wb.matrix[wb.matrix == wb.current_edge_label] = -(wb.current_edge_label -1)
					wb.matrix[wb.matrix == -wb.current_edge_label] = 0
					# edge placement number WorkingBasket.current_edge_label is decremented
					wb.current_edge_label -= 1
					# progress is saved / WorkingBaseket is pickled
					saver.save()

				else:
					# if the most recently placed edge looks good according to the heuristics then:
					# WorkingBasket.current_edge_label is incremented
					wb.current_edge_label += 1
					# the new_edge is labelled with wb.current_edge_label
					wb.label(new_edge, wb.current_edge_label)
					# mark newly unavailable edges with -WorkingBasket.current_edge_label
					label_non_edges(new_edge, -wb.current_edge_label)
					# progress is saved / WorkingBasket is pickled
					if (wb.current_edge_label % 1000) == 0:
						saver.save()

class HeuristicConductor:
	def __init__(self):
		global wb

	def passes(self):
		# return True iff all huristic checks say wb may still be a feasible subgraph
		return feasibility_check_1() and feasibility_check_2()
	pass #to do

def choose_new_edge(mode):
	global wb

	if not any((any((x == 0 for x in wb.matrix[y])) for y in range(wb.v))):
		raise EdgesExhausted('Edges exhasted')
	else:
		# select the new edge based on mode (deep, wide, or random). Use deep by default. Return a tuple.
		if mode == 'deep':
			# return a vertex that has the least number of potential edge ends
			pot_degree_list = [wb.pot_degree_of(v) for v in range(wb.v)]
			v0 = pot_degree_list.index(min([x for x in pot_degree_list if x !=0]))
			v1 = [v for v in wb.matrix[v0]].index(0)
			return (v0,v1)
		else:
			sys.exit("[!] choose_new_edge(mode) has invalid mode set.")

class EdgesExhausted(Exception):
	# this exception is raised by choose_new_edge if there are no avaialable edges from which to choose
	pass

class Saver:
	def __init__(self, save_option):
		self.save_option = save_option
		global wb

	def save(self):
		if self.save_option:
			filename = 'moore_search-val-'+ str(wb.val) + '-' + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")+".p"
			pickle.dump(wb, open(filename, "wb"))

class Verifier:
	def __init__(self, valence, candidate_graph):
		self.candidate_graph = candidate_graph # expected to be an instantiation of the CandidateGraph class
		self.val = valence # valence of the target Moore graph
		self.v = 1 + self.val + self.val * (self.val - 1) # number of vertices in the target Moore graph
		self.e = self.v * self.val / 2 # number of edges in the target Moore graph
		self.fruit = tuple(tuple(i*(self.val-1)+j for j in range(0,self.val-1)) for i in range(0,self.val))

	def verify(self):
		# diameter 2
		self.diameter_check = self.candidate_graph.graph.diameter() == 2
		# girth 5
		self.girth_check = self.candidate_graph.graph.girth() == 5
		# degree self.val
		self.degree_check = all([d == self.val for d in self.candidate_graph.graph.degree()])
		# order self.v
		self.order_check = self.candidate_graph.graph.vcount() == self.v
		# edge count
		self.edge_check = self.candidate_graph.graph.ecount() == self.e
		# return Moore graph staus
		return self.diameter_check and self.girth_check and self.degree_check and self.order_check and self.edge_check

class CandidateGraph:
	def __init__(self, basket):
		self.basket = basket # the WorkingBasket's self.matrix
		self.basket_v = len(basket[:,1]) # number of vertices in the basket
		self.val = int((sqrt(4 * self.basket_v + 1) + 1) / 2) # expected valence of the target Moore graph
		self.fruit = tuple(tuple(i*(self.val-1)+j for j in range(0,self.val-1)) for i in range(0,self.val))

		self.v = 1 + self.val + self.val * (self.val - 1) # number of vertices in the complete graph
		self.build_adj_mat() # construct the adjacency matrix of the hopeful Moore graph
		self.graph = ig.Graph.Adjacency(self.adj_mat.tolist()) # igraph graph constructed from self.adj_mat
		self.graph.to_undirected()

	def build_adj_mat(self):
		# use self.basket to create an adjacency matrix for the complete candidate graph 
		self.adj_mat = np.vectorize(self._project, otypes=[np.int8])(self.basket) # gives the adj matrix of the basket
		# add (self.val + 1) more rows/cols to adj_mat, all zeroed out
		b = np.zeros((self.v,self.v))
		b[:-(self.val + 1),:-(self.val + 1)] = self.adj_mat
		self.adj_mat = np.array(b, dtype=np.int32)
	
		parent_verts = tuple(i for i in range(self.basket_v,self.basket_v + self.val))
		# connect each parent vertex to each vertex in it's corresponding fruit
		parent_edges = it.chain.from_iterable(tuple(tuple((parent_verts[p],v) for v in self.fruit[p]) for p in range(self.val)))
		[self._add_edge(edge) for edge in parent_edges]
		# connect the root vertex to each parent vertex
		root_edges = tuple((self.v - 1, parent_verts[p]) for p in range(self.val))
		[self._add_edge(edge) for edge in root_edges]

	def _project(self, value):
		if value > 0:
			return 1
		else:
			return 0

	def _add_edge(self, edge):
		# adds an edge to self.adj_mat
		self.adj_mat[edge[0],edge[1]] = 1
		self.adj_mat[edge[1],edge[0]] = 1

def endgame():
	global wb
	candidate_graph = CandidateGraph(wb.matrix)
	verifier = Verifier(wb.val, candidate_graph)
	if verifier.verify():
		# we found a Moore graph!
		# save the adjacencyMatrix
		filename = 'moore_graph-val-' + str(wb.val) + "-" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")+".p"
		pickle.dump(list(candidate_graph.graph.get_adjacency()) , open(filename, "wb"))
		sys.exit("[!] A Moore graph was found and saved.")
	else: # not a moore graph.
		sys.exit("[!] No Moore graph was found.")

def label_non_edges(new_edge, label):
	# Note, label should be a negative integer.
	# Calling this function should always be imediately preceeded by calling wb.label(new_edge, wb.current_edge_label)
	# This function assumes new_edge has just been added to WorkingBasket
	# The goal of this function is to find all of the edges that are now ruled out because of the placement of new_edge (and that
	#  had *not* already been ruled out prior to the placement of new_edge) and label those new non-edges with -wb.current_edge_label
	# This is one area of code that may benefit from parallelisation
	global wb
	## Make negative the edges that would cause a vertex to have degree > wb.val-1
	#    note that the addition of new_edge will only have raised the degree of the vertices new_edge[0] and new_edge[1], so we need
	#     only check those vertices
	for edge_end in range(2):
		if wb.degree_of(new_edge[edge_end]) == wb.val -1: # then no more edges can be added to vertex new_edge[0]
			# so we find all potential edges at new_edge[0] and label them -wb.current_edge_label
			[wb.label((new_edge[edge_end],v), label) for v in np.where(wb.matrix[new_edge[edge_end]] == 0)[0]]

	## Make negative all other edges from v0 to fruit_of(v1), and all other edges from v1 to fruit_of(v0) !!
	for edge_end in range(2):
		[wb.label((new_edge[edge_end],v),-wb.current_edge_label)for v in fruit_of(new_edge[(edge_end + 1)%2]) if wb.matrix[new_edge[edge_end]][v] == 0]

	## Make negative the edges that would create triangles
	#    begin with the vertex new_edge[0] and look at every vertex already incident to edge[0]
	#    that collection of vertices is np.where(wb.matrix[new_edge[0]] > 0)[0]
	# if wb.matrix[vertex][new_edge[1]] is zero the mark that edge negative (because adding it would cause a triangle to form)
	for edge_end in range(2):	
		[wb.label((neighbor,new_edge[(edge_end + 1)%2]),label) for neighbor in np.where(wb.matrix[new_edge[edge_end]] > 0)[0] if wb.matrix[neighbor][new_edge[(edge_end + 1)%2]] == 0]
	## Make negative the edges that would create rectangles (three cases)
	# case 1
	#	 suppose new_edge = (v0, v1)
	#    find the neighbors(v0). They are nlabel
	#    find the neighbors(v1). They are np.where(wb.matrix[new_edge[1]] > 0)[0]
	#	 we want to rule out all edges (i,j) where i is in neighbors(v0) and j is in neighbors(v1), because those edges would cause
	#     a rectangle to form
	[wb.label((i,j),label) for (i,j) in ((i,j) for i in np.where(wb.matrix[new_edge[0]] > 0)[0] for j in np.where(wb.matrix[new_edge[1]] > 0)[0]) if wb.matrix[i][j] == 0]

	# cases 2 & 3
	#    suppose new_edge = (v0, v1)
	#    find the neighbors of the neighbors of v0
	for edge_end in range(2):	
		neighbors_of_v0 = set()
		neighbors_of_v0.update([n for n in np.where(wb.matrix[new_edge[edge_end]] > 0)[0] if n != new_edge[edge_end]]) # exclude new_edge[1] for efficiency
		array_of_nn_of_v0 =  [np.where(wb.matrix[neighbor] > 0)[0] for neighbor in neighbors_of_v0] # array of neighbors of neighbors of v0
		nn_of_v0 = set() # neighbors of neighbors of v0
		[nn_of_v0.update(x) for x in array_of_nn_of_v0]
		nn_of_v0.discard(new_edge[edge_end]) #neighbors of neighbors of v0 excluding v0
		# rule out all edges (i,j) where i is in neighbors(v0) and j is in neighbors(v1) (if they haven't been ruled out already), 
		#  because those edges would cause a rectangle to form
		[wb.label((x,new_edge[(edge_end + 1)%2]),label) for x in nn_of_v0 if wb.matrix[x][new_edge[(edge_end + 1)%2]] == 0]

def feasibility_check_1():
	# return True iff for every vertex v, and every fruit f not containing v, there is either already an edge from v to f or else there
	#  is a potential edge (labeled '0' in WorkingBasket) from v to f
	global wb
	return all([all([any(wb.matrix[f,v] >= 0) for f in wb.fruit if v not in f]) for v in range(wb.v)])

def feasibility_check_2():
	# return True iff for every v0 and every v1 not in the fruit of v0, there exists a path of length at most two (through either 
	#  already placed edges and/or potential edges labeled 0).
	global wb
	neighbors = wb.flat_matrix() # a boolean adjacency matrix -- True wherever an already placed or a potential edge exists, False elsewhere
	neighbors_of_neighbors = neighbors.dot(neighbors) # i,j is True if there exists a path of length 2 from i to j
	intra_fruit_matrix = np.zeros(shape=(wb.v, wb.v), dtype=bool)
	[[[label_true(intra_fruit_matrix,e[0],e[1]) for e in [(i,j) for i in fruit for j in fruit]] for v in fruit] for fruit in wb.fruit]
	return (neighbors + neighbors_of_neighbors + intra_fruit_matrix).all()

def accelerate():
	# for all v in fruit(0), place an edge from v to the vth edge in every other fruit
	for v in wb.fruit[0]:
		for f in range(1,wb.val):
			wb.current_edge_label += 1
			wb.label((v,wb.fruit[f][v]), wb.current_edge_label)
			label_non_edges((v,wb.fruit[f][v]), -wb.current_edge_label)
	print 'accelerate() complete'

def label_true(matrix, i, j):
	matrix[i][j] = matrix[j][i] = True

def fruit_of(vertex):
	# returns the fruit of the input vertex
	global wb
	return wb.fruit[int(floor(vertex/(wb.val - 1)))]

def town_of(v0):
	# returns a set of all current and potential neighbors of vertex, as well as all current and potential neighbors of
	#  neighbors of vertex. This is used to determine whether another vertex is within (or *can be* within) distance two of vertex
	global wb
	neighbors_of_v0 = set() # potential neighbors of v0
	neighbors_of_v0.update([n for n in np.where(wb.matrix[v0] >= 0)[0]])
	nn_of_v0 = set() # potential neighbors of neighbors of v0
	[nn_of_v0.update(x) for x in [np.where(wb.matrix[neighbor] >= 0)[0] for neighbor in neighbors_of_v0]]
	town_of_v0 = neighbors_of_v0 | nn_of_v0 # set of all vertices that are potential or current neighbors or neighbors of neighbors of v0
	return town_of_v0

def test_setup(): # todo remove this function
	global wb
	global new_edge
	wb.current_edge_label = 99
	wb.label((0,1),4)
	wb.label((0,2),4)
	wb.label((2,3),4)
	new_edge = (0,2)

def test_setup_2(): # todo remove this function
	global wb
	global candidate_graph
	global verifier
	candidate_graph = CandidateGraph(wb.matrix)
	verifier = Verifier(wb.val, candidate_graph)	


if __name__ == "__main__":
	# if script is executed at the CLI parse the CLI arguments
	usage = '''usage: %(prog)s [-n <valence>] [-a] [-s][-c <pickle_file.p>]'''
	parser = argparse.ArgumentParser(usage=usage)
	parser.add_argument("-n", "--new", type=int, help="Begin a new search for a Moore of valence NEW_VAL", action="store", dest="new_val")
	parser.add_argument("-a", "--accelerate", help="Place known edges in the basket without checking heuristics when ACCELERATE option is set. Speeds up search without any risk of failure.", action="store_true", dest="accelerate")
	parser.add_argument("-s", "--save", help="Saves progress of search so may be continued later. Ideal for large valence graphs.", action="store_true", dest="save_option")
	parser.add_argument("-c", "--continue", type=str, help="Continue the search saved in CONTINUE_FILE", action="store", dest="continue_file")
	args = parser.parse_args()

	# verify consistency of options selection
	if not (args.new_val or args.continue_file):
		sys.exit("[!] Please choose either -n or -c to run this script. See --help for more details.")
	if args.accelerate and not args.new_val:
		sys.exit("[!] To use the --accelerate feauture you must specify a valnce using -n. See --help for more details.")
	if (args.new_val or args.accelerate) and args.continue_file:
		sys.exit("[!] The -c (or --continue) feature cannot be used with -n or -a. See --help for more details.")

	# hand off control to the initialiser
	initialiser = Initialiser(args.new_val, args.accelerate, args.continue_file, args.save_option)
	initialiser.initialise()
