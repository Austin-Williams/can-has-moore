import numpy as np
import itertools as it
import igraph as ig

class WorkingBasket:
	def __init__(self, valence):
		self.val = valence # valence of target Moore graph
		self.v = self.val * (self.val - 1) # number of vertices in the basket
		self.e = self.v * self.val / 2 # number of edges in the target basket
		self.fruit = tuple(tuple(i*(self.val-1)+j for j in range(0,self.val-1)) for i in range(0,self.val))

	def new_basket(self):
		self.matrix = np.zeros(shape=(self.v, self.v), dtype=np.int32) # working basket zeroed out
		# rule out inter-fruit edges and self-loops
		inter_fruit = it.chain.from_iterable([it.combinations(self.fruit[f],2) for f in range(self.val)])
		[self.label(edge, -(self.e + 1)) for edge in inter_fruit]
		[self.label(edge, -(self.e + 1)) for edge in [(i,i) for i in range(self.v)]]
		self.current_edge = 0 # label number of the most recent edge placed

	def label(self, edge, label):
		# changes the label of an edge in self.matrix
		self.matrix[edge[0],edge[1]] = label
		self.matrix[edge[1],edge[0]] = label

	def construct_full_graph(self):
		# outputs the adjacency matrix of the graph (inluding root and parent vertices)
		#  of the graph currently represented by self.matrix

		# convert 
		adj_mat = np.vectorize(self._flatten, otypes=[np.int8])(self.matrix)
		pass

	def flat_matrix(self):
		# return a uint8 np array that 'flattens' self.matrix
		# results in a much smaller matrix (1/4th the size) for use in parallelisation of
		#  heuristics
		return np.vectorize(self._flatten, otypes=[np.int8])(self.matrix)

	def _flatten(self, value):
		if value == 0:
			return 0
		else:
			return value / abs(value)

class Verifier:
	def __init__(self, valence, basket):
		self.val = valence # valence of the target Moore graph
		self.v = 1 + self.val + self.val * (self.val - 1) # number of vertices in the target Moore graph
		self.e = self.v * self.val / 2 # number of edges in the target Moore graph
		self.basket = basket # the WorkingBasket's self.matrix
		self.basket_v = self.val * (self.val - 1) # number of vertices we expect to see in the basket
		self.fruit = tuple(tuple(i*(self.val-1)+j for j in range(0,self.val-1)) for i in range(0,self.val))
		self.build_adj_mat() # construct the adjacency matrix of the hopeful Moore graph
		self.graph = ig.Graph.Adjacency(self.adj_mat.tolist()) # igraph graph constructed from self.adj_mat

	def verify(self):
		# diamter 2
		self.diameter = self.graph.diameter()
		self.diameter_check = self.diameter == 2
		# girth 5
		self.girth = self.graph.girth()
		self.girth_check = self.girth == 5
		# degree self.val
		self.degree_check = all([d == self.val for d in self.graph.degree()])
		# order self.v
		self.order_check = self.graph.vcount == self.v
		# return Moore graph staus
		return self.diameter_check and self.girth_check and self.degree_check and self.order_check

	def build_adj_mat(self):
		# use self.basket to create an adjacency matrix for the full candidate graph 
		self.adj_mat = np.vectorize(self._project, otypes=[np.int8])(self.basket) # gives the adj matrix of the basket
		self.adj_mat.resize((self.v, self.v)) # adds (self.val + 1) more rows/cols to adj_mat, all zeroed out
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

class Initialiser:
	# from user input, decide whether to start a new search or continue an old one

		# if starting a new search then 
			# instantiate a new WorkingBasket with proper valence
			# have the WorkingBaket run it's new_basket method
			# Optional: run an accelerator to get the known-present edges in without running through the huristic checks

		# if continuing an old one then
			# unpickle the matrix
			# compute the valence
			# instantiate a new WorkingBasket with proper valence
			# compute the follwoing from the unpickled file and store the values in the WorkingBasket
				# self.matrix
				# self.current_edge
			# Optional: check to make sure the loaded subgraph is feasible, otherwise trust it

	# indicate that all went well and we're ready to have the conductor start working

	pass

class Conductor:
	# HuristicChecker runs heuristic checks to see if current basket can be ruled out as a subgraph
		# if it can be ruled out then:
			# edge placement number WorkingBasket.current_edge is marked and reverted
			# edge placement number WorkingBasket.current_edge is decremented
			# progress is saved / WorkingBaseket is pickled

		# if it cannot be ruled out, then:
			# WorkingBasket.current_edge is incremented
			# EdgePicker tries and picks an available edge to label WorkingBasket.current_edge
				# If there are no more available edges to pick (raise Exception), then we're finished and we need to check
				#  whether we've found a Moore graph or whether one isn't possible.
					# In this event, it's worthwhile to save / pickle the progress before evaluating.
			# Helpers mark newly unavailable edges with -WorkingBasket.current_edge
			# progress is saved / WorkingBaseket is pickled
	# repeat 
	pass

if __name__ == "__main__":
	# for quick / easy testing, will remove later
	b = WorkingBasket(3)
	b.new_basket()
	v = Verifier(3, b.matrix)
