#!/opt/local/bin/python3.5
"""
	Filter points
	Objective of this script is to optimize for the highest point count given restriction that 
		there is constraint on minimal spatial distance between any tow points 


	Algorithm overview:
		Operates on the fact that there is no better solution to brake connected graph while minimizing 
			number of nodes deleted than deletion of the nodes with higher degree first

		Points are represented using a undirectional graph in which each edge represents a spatial distance 
			lower than allowed by constraint. In order to construct such graph structure weighted Laplacian 
			symmetrical matrix was constructed and distance between nodes were used as weights.


		In order to optimize for performance only paths among nodes possessing edges were calculated and
			each path was processed separately.
			

	Author:  Hrvoje Novosel<hrvojedotnovose@gmail.com>
	Created: Tue Aug 16 22:36:03 CEST 2016
	Tested and developed using Python 3.5.1 (Darwin Kernel Version 15.3.0 x86_64)

"""

import os
import sys
import numpy as np
from multiprocessing import Pool
from lib.geo import geoclass as geo

# filter distance in meters (radius)
MIN_DISTANCE = 1500


def read_points(fname):
	"""
		Import file structure:
			Source structure: csv file, having variable names in the first row and
			three columns in the following manner: pt_id, pt_lat, pt_lon
			CSV Deliminetor: comma
	"""
	ret = []
	
	with open(fname, 'r') as fp:
		fp.readline()
		for ln in fp:
			# ok, don't have enough time for generic solution
			fields = ln.rstrip().split(',')
			# appending should go lat, long, ....
			ret.append( [ float(fields[2]), float(fields[1]), int(fields[0]) ] )

	print("Record count of loaded dataset: {}".format(len(ret)))
	return ret


def save_points(fname, pts_lst, selected):
	"""
		Quick save of the results
	"""
	with open(fname, 'w') as fp:
		fp.write('y;x;id\n')
		for i in selected:
			fp.write(';'.join([str(x) for x in pts_lst[i]]))
			fp.write('\n')


def calculate_distances(vals):
	"""
		Calculates haversine distance between two points
	"""
	pt_a = vals[1]
	# has every other pts!
	pts  = vals[2]
	ret_lst = []

	for i in range(len(pts)):
		dist = geo.haversine_distance( \
				pt_a[0], pt_a[1], \
				pts[i][0], pts[i][1]
			)
		if dist < MIN_DISTANCE:
			ret_lst.append(( (i * 2 + 1), dist))

	return (vals[0], ret_lst)


def create_adj_matrix(pts_list):
	"""
		Calculates laplacian adjacency matrix, placing distance only to edges in violation of min distance constraint
	"""
	adj_mat = np.zeros([len(pts_lst), len(pts_lst)], dtype=np.int)

	# lets go parallel
	ivals = range(0, len(pts_list), 2)
	jvals = [pts_list[x] for x in range(1, len(pts_list), 2)]
	
	pool = Pool(10)
	prog = 0
	for val in pool.map(calculate_distances, [(x, pts_list[x], jvals) for x in ivals]):
		# print(val)
		if len(val[1]) > 0:
			i = val[0]
			for dist in val[1]:
				j = dist[0]
				adj_mat[i][j] = dist[1]
				adj_mat[j][i] = dist[1]
		# useful progress message
		prog_precent = prog/len(ivals) * 100
		if prog_precent % 10 == 0:
			sys.stdout.write('\rCreating distance matrix: {}%'.format(prog_precent))
			sys.stdout.flush()
		prog += 1

	sys.stdout.write('\rCreating distance matrix: 100.0%\n')
	return adj_mat


def create_path(node, adj_mat, path):
	"""
		Construct path for all nodes returning number of elements in the path?
	"""
	# recursive build up of path
	# how graph is undirectional order doesn't matter
	for nodes in np.nonzero(adj_mat[node]):
		for n in nodes:
			if not n in path:
				path.append(n)
				path = create_path(n, adj_mat, path)

	# print('node: {} -> path: {}'.format(node, path))
	return path


def pop_edge_node(node_lst, adj_mat):
	"""
		Recursively pops nodes out of list until all edges are removed
	"""
	for n1 in node_lst:
		for n2 in node_lst:
			if adj_mat[n1][n2] > 0:
				node_lst.pop(0)
				return pop_edge_node(node_lst, adj_mat)

	return node_lst


if __name__ == '__main__':

	if len(sys.argv) != 4:
		# just a quick debug line printout
		if len(sys.argv) == 2:
			print('./filter.py sample_pts.csv filtered_pts.csv 1500')
		else:
			print('Usage:')
			print('python filter.py input_ds output_ds constraining_distance (m)')
			print('eg:')
			print('python filter.py input.csv output.csv 1500')
		sys.exit()

	# globals
	FNAME = sys.argv[1]
	SNAME = sys.argv[2]

	# filter distance in meters (radius)
	try:
		MIN_DISTANCE = float(sys.argv[3])
	except:
		print('Incorrect distance value')
		sys.exit()


	pts_lst = read_points(FNAME)
	adj_mat = create_adj_matrix(pts_lst)


	# split matrix in two parts, ones without any point in proximity and ones for the graph
	selected   = [] # list of selected object ids
	unselected = [] # list of pointers to adj_matrix
	for i in range(len(adj_mat)):
		if np.sum(adj_mat[i]) == 0:
			selected.append(i)
		else:
			unselected.append(i)


	print('Number no constraints infringement selected nodes: {}'.format(len(selected)))

	# better solution sorting schema would be on degree 
	unselected = sorted( unselected, key=lambda x: np.count_nonzero(adj_mat[x]) )


	processed = [False for _ in range(len(unselected))]
	while False in processed:

		node = unselected[processed.index(False)]
		path = create_path(node, adj_mat, [node])

		# sort found path according to the node degree
		degree_sorted = sorted(path, key=lambda x: unselected.index(x) * -1)

		# print (degree_sorted)
		
		# start poping nodes until no edges has been left
		degree_sorted.pop(0)
		degree_sorted = pop_edge_node(degree_sorted, adj_mat)
		for nod in degree_sorted:
			selected.append(nod)

		# print (degree_sorted)


		# set the flags
		for p in path:
			processed[unselected.index(p)] = True
		# break
		
	# and finally process the results
	print('Total number of selected points {}, number of discarded points {}'.format(len(selected), len(pts_lst) - len(selected)))
	save_points(SNAME, pts_lst, selected)
	print('Points have been saved to file ' + SNAME)
	print('Done :)')
