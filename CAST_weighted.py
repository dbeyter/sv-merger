import os
import sys
import operator


# CAST (Cluster Affinity Search Technique) algorithm.

# active_elements: bool list with true of size of the num of vertices. 
# true becomes false when vertex is added to a clique. 
# distance_matrix: list of lists for all pairwise distance between elements (N by N)
# distance_thresh: distance threshold.

def find_cliques(dm_map, sv_pns, distance_matrix, distance_thresh):

	# if len(distance_matrix) == 1:
	# 	return([set([0])])
	# else:
	# 	print "DM: " + str(distance_matrix)

	list_of_cliques = []

	# while not all elements added to a clique
	active_elements_map = [True for el in dm_map]
	clique_memberships = [False for el in dm_map]

	vertex_pn_scoremap = [ [0]*len(set(sv_pns)) for k in range(len(dm_map)) ]
	vertex_pn_edgenummap = [ [0]*len(set(sv_pns)) for k in range(len(dm_map)) ]
	vertex_total_uniqpns = [0]*len(dm_map)
	vertex_totalscore = [0]*len(dm_map)

	cli_r = 0
	while sum(active_elements_map) > 0:
		ADDING_REQUIRED = True
		#print( "Current active els: " + str(sum(active_elements_map)))

		active_indices = set([ind for ind in range(len(active_elements_map)) if active_elements_map[ind]])

		v = find_maximal_degree_active_vertex(active_indices, distance_matrix, dm_map, sv_pns, distance_thresh)

		seen_cliques_set = set()
		cli_r = cli_r ^ (1 << v)
		seen_cliques_set.add(cli_r)
		clique_elements_sofar = set([v])

		# init
		v_pn = sv_pns[v]
		for el in range(len(vertex_pn_scoremap)):
			if el != v:
				vertex_pn_scoremap[el][v_pn] = distance_matrix[dm_map[el]][dm_map[v]]
				vertex_totalscore[el] = vertex_pn_scoremap[el][v_pn]
				vertex_pn_edgenummap[el][v_pn] = 1
				vertex_total_uniqpns[el] = 1
		#print( " START: " + str(clique_elements_sofar))

		close_not_in = True
		far_not_out = True
		
		
		not_in_clique = active_indices.difference(clique_elements_sofar)


		while close_not_in or far_not_out:

			if ADDING_REQUIRED:

				min_dist = float('Inf')
				closest_el = -1
				for el in not_in_clique:
					el_dist = vertex_totalscore[el]
					#print "==NOT IN CLIQUE== clique distance of " + str(el) + " is: " + str(el_dist)
					if el_dist <= distance_thresh and el_dist < min_dist:
						min_dist = el_dist
						closest_el = el
				
				#print( "Closest el is: " + str(closest_el))
				#print("Min distance is: " + str(min_dist))
				if closest_el != -1:
					clique_elements_sofar.add(closest_el)
					#clique_memberships[closest_el] = True
					not_in_clique.remove(closest_el)


					# Loop prevention
					cli_r = cli_r ^ (1 << closest_el)
					if cli_r not in seen_cliques_set:
						seen_cliques_set.add(cli_r)
					else:
						#print "!!!!!!!!!!!!!!!!!!!THE SAME CLIQUE HAS BEEN OBSERVED BEFORE AT THE ADDITION!!!!!!!!!!!!!!!!!!!"
						ADDING_REQUIRED = False



					clique_el = closest_el
					clique_el_pn = sv_pns[closest_el]	
					
					#update scores
					for el in active_indices:
						# print("el is: " + str(el))
						# print("closest_el is: " + str(closest_el))
						# print("dm_map[el]: " + str(dm_map[el]))
						# print("dm_map[closest_el]: " + str(dm_map[closest_el]))
						# print("clique_el_pn: " + str(clique_el_pn)) 
						if el != clique_el:
							vertex_pn_scoremap[el][clique_el_pn] , vertex_totalscore[el], vertex_pn_edgenummap[el][clique_el_pn], vertex_total_uniqpns[el] = \
								element_changed_in_cluster_score_updates(True, el, sv_pns[el], clique_el, clique_el_pn, distance_matrix[dm_map[el]][dm_map[clique_el]], \
								vertex_pn_scoremap[el][clique_el_pn], vertex_pn_edgenummap[el][clique_el_pn], vertex_totalscore[el], vertex_total_uniqpns[el])
						
							#update cluster pn count
							#vertex_pn_scoremap[el][clique_el_pn] = new_vertex_pn_score
							#vertex_totalscore[el] = new_vertex_cluster_score
							

							# See if correct.
							#vertex_pn_edgenummap[el][clique_el_pn] += 1
							#if vertex_pn_edgenummap[el][clique_el_pn] == 1:
							#	vertex_total_uniqpns[el] += 1

					#update pn_weights
					#print("Selected pn: " + str(clique_el_pn))
					#pn_weights[clique_el_pn] += 1
					#cluster_pn_num = cluster_pn_num + old_cluster_pn_num_change

					close_not_in = ADDING_REQUIRED
					far_not_out = True


				else:
					close_not_in = False

				#print("(added?) Clique elements so far: " + str(clique_elements_sofar))
				#print("weights so far: " + str(pn_weights))

			## check if there is a distant gene in C.
			max_dist = float('-Inf')
			farthest_el = -1
			for el in clique_elements_sofar:
				el_dist = vertex_totalscore[el]
				#print("==IN CLIQUE== clique distance of " + str(el) + " is: " + str(el_dist))
				if el_dist > distance_thresh and el_dist > max_dist:
					max_dist = el_dist
					farthest_el = el

			#print("Farthest el is: " + str(farthest_el))
			#print("Max distance is: " + str(max_dist))
			if farthest_el != -1:
				clique_elements_sofar.remove(farthest_el)
				#clique_memberships[farthest_el] = False
				not_in_clique.add(farthest_el)

				# Loop prevention
				cli_r = cli_r ^ (1 << farthest_el)
				if cli_r not in seen_cliques_set:
					seen_cliques_set.add(cli_r)

				clique_el = farthest_el
				clique_el_pn = sv_pns[farthest_el]
				# if pn_weights[clique_el_pn] > 1:
				# 	old_cluster_pn_num_change = 0
				# else:
				# 	old_cluster_pn_num_change = -1

				#update scores
				for el in active_indices:
					if el != clique_el:

						vertex_pn_scoremap[el][clique_el_pn] , vertex_totalscore[el], vertex_pn_edgenummap[el][clique_el_pn], vertex_total_uniqpns[el] = \
							element_changed_in_cluster_score_updates(False, el, sv_pns[el], clique_el, clique_el_pn, distance_matrix[dm_map[el]][dm_map[clique_el]], \
							vertex_pn_scoremap[el][clique_el_pn], vertex_pn_edgenummap[el][clique_el_pn], vertex_totalscore[el], vertex_total_uniqpns[el])
							#update cluster pn count
							#vertex_pn_scoremap[el][clique_el_pn] = new_vertex_pn_score
							#vertex_totalscore[el] = new_vertex_cluster_score
							#vertex_pn_edgenummap[el][clique_el_pn] -= 1
							#if vertex_pn_edgenummap[el][clique_el_pn] == 0:
							#	vertex_total_uniqpns[el] -= 1

				#update pn_weights
				# Cannot not exist so just remove the weight by 1
				#pn_weights[clique_el_pn] -= 1
				
				#cluster_pn_num = cluster_pn_num + old_cluster_pn_num_change


				far_not_out = True
				close_not_in = ADDING_REQUIRED

			else:
				far_not_out = False
			#print("(removed?) Clique elements so far: " + str(clique_elements_sofar))
			#print("weights so far: " + str(pn_weights))

			#print("close_not_in: " + str(close_not_in))
			#print("far_not_out: " + str(far_not_out))
				
				
		# add the clique_elements_sofar to the list of cliques because it's done.
		list_of_cliques.append(clique_elements_sofar)
		
		#print("Final clique is: " + str(clique_elements_sofar))
		
		# any element that is added is no longer active.
		for el in clique_elements_sofar:
			active_elements_map[el] = False

		# Clean the remaining els for their total scores, and pns encountered in clique.
		#print("new init begin")
		for el in not_in_clique:
			for el_c in clique_elements_sofar:
				pn_c = sv_pns[el_c]
				vertex_pn_scoremap[el][pn_c] = 0
				vertex_pn_edgenummap[el][pn_c] = 0
			vertex_total_uniqpns[el] = 0
			vertex_totalscore[el] = 0
		#print("new init end")
		# vertex_pn_scoremap = [ [0]*len(set(sv_pns)) for k in range(len(dm_map)) ]
		# vertex_pn_edgenummap = [ [0]*len(set(sv_pns)) for k in range(len(dm_map)) ]
		# vertex_total_uniqpns = [0]*len(dm_map)
		# vertex_totalscore = [0]*len(dm_map)



	return(list_of_cliques)


def find_maximal_degree_vertex(active_elements_map, distance_matrix, dm_map, sv_pns, distance_thresh):

	element_degrees = [0 for el in active_elements_map]
	for i in range(len(active_elements_map)):
		if active_elements_map[i]:
			near_pns = [ sv_pns[j] for j in range(len(dm_map)) if (distance_matrix[dm_map[i]][dm_map[j]]<=distance_thresh and active_elements_map[j])   ]
			set_near_pns = set(near_pns)
			#print("i: " + str(i) + " near_vertices_len: " + str(len(near_pns)) + " uniq pn count: " + str(len(set_near_pns)))
			element_degrees[i] = len(set_near_pns)

	# done.
	# Get the max degree vertex.
	max_index, max_value = max(enumerate(element_degrees), key=operator.itemgetter(1))
	return(max_index)



def find_maximal_degree_active_vertex(active_indices, distance_matrix, dm_map, sv_pns, distance_thresh):

	#element_degrees = [0 for el in range(len(sv_pns))]
	#print("max_deg begin")
	max_el = -1
	max_degrees = float('-Inf')
	for i in active_indices:

			near_pns = [ sv_pns[j] for j in active_indices if (distance_matrix[dm_map[i]][dm_map[j]]<=distance_thresh)   ]
			set_near_pns = set(near_pns)
			#print("i: " + str(i) + " near_vertices_len: " + str(len(near_pns)) + " uniq pn count: " + str(len(set_near_pns)))
			if len(set_near_pns) > max_degrees:
				max_el = i
				max_degrees = len(set_near_pns)

	# done.
	#print("max_deg end")
	# Get the max degree vertex.
	#max_index, max_value = max(enumerate(element_degrees), key=operator.itemgetter(1))
	return(max_el)


#element_changed_in_cluster_score_updates(addition, vertex, vertex_pn, clique_el, clique_el_pn, d_vertex_clique_el, vertex_pn_score, vertex_to_pn_edgenum, vertex_cluster_score, vertex_to_uniqpn_num)
def element_changed_in_cluster_score_updates(addition, vertex, vertex_pn, clique_el, clique_el_pn, d_vertex_clique_el, vertex_pn_score, vertex_to_pn_edgenum, vertex_cluster_score, vertex_to_uniqpn_num):
	
	#print("vertex: " + str(vertex))
	#print("clique_el: " + str(clique_el))  
	
	total_vertex_score = vertex_cluster_score * vertex_to_uniqpn_num
	#print("vertex_cluster_score: " + str(vertex_cluster_score)) 
	#print("vertex_to_uniqpn_num: " + str(vertex_to_uniqpn_num)) 
	#print("total_vertex_score (mult): " + str(total_vertex_score)) 

	total_vertex_score -= vertex_pn_score
	#print("vertex_pn_score: " + str(vertex_pn_score)) 
	#print("left total_vertex_score: " + str(total_vertex_score))

	total_vertex_pn_score = vertex_pn_score * vertex_to_pn_edgenum
	#print("vertex_pn_score: " + str(vertex_pn_score)) 
	#print("vertex_to_pn_edgenum: " + str(vertex_to_pn_edgenum)) 
	#print("total_vertex_pn_score (mult): " + str(total_vertex_pn_score)) 

	# adding a vertex to the cluster.
	if addition:
		new_vertex_to_pn_edgenum = vertex_to_pn_edgenum + 1
		if new_vertex_to_pn_edgenum == 1:
			new_vertex_to_uniqpn_num = vertex_to_uniqpn_num + 1
		else:
			new_vertex_to_uniqpn_num = vertex_to_uniqpn_num
		
		#print("new_vertex_to_pn_edgenum: " + str(new_vertex_to_pn_edgenum)) 
		new_vertex_pn_score = (total_vertex_pn_score + d_vertex_clique_el)/float(new_vertex_to_pn_edgenum)

		#print("d_vertex_clique_el: " + str(d_vertex_clique_el) )
		#print("new_vertex_pn_score: " + str(new_vertex_pn_score) )

	# removing a vertex from the cluster.
	else:
		new_vertex_to_pn_edgenum = vertex_to_pn_edgenum - 1
		if new_vertex_to_pn_edgenum == 0:
			new_vertex_pn_score = 0
			new_vertex_to_uniqpn_num = vertex_to_uniqpn_num - 1
		else:
			new_vertex_pn_score = (total_vertex_pn_score - d_vertex_clique_el)/float(new_vertex_to_pn_edgenum)
			new_vertex_to_uniqpn_num = vertex_to_uniqpn_num 


	new_vertex_cluster_score = (total_vertex_score + new_vertex_pn_score)/float(new_vertex_to_uniqpn_num)
	#print("new_vertex_to_uniqpn_num: " + str(new_vertex_to_uniqpn_num))
	#print("total_vertex_score + new_vertex_pn_score: " + str(total_vertex_score + new_vertex_pn_score) )
	#print("new_vertex_cluster_score: " + str(new_vertex_cluster_score) )


	#update what needs updating.
	return(new_vertex_pn_score, new_vertex_cluster_score, new_vertex_to_pn_edgenum, new_vertex_to_uniqpn_num)


