

# On the TRF conbined csv.
CHR_COL = 0
CLUSTER_BEGIN_COL = 1
CLUSTER_END_COL = 2
SV_IDX_COL = 3
SAMPLE_COL = 4
METHOD_COL = 5
EVENT_COL = 6
CLUSTERID_COL = 7


def read_pieced_data(part_cluster_fn):
	sv_cluster_dict = {}
	
	part_cluster_f = open(part_cluster_fn, 'r')
	for l in part_cluster_f:
		if l[0] != '#':
			els = l.strip().split('\t')
			
			index = els[CLUSTERID_COL]
			
			chrom = els[CHR_COL]
			cluster_begin = int(els[CLUSTER_BEGIN_COL])
			cluster_end = int(els[CLUSTER_END_COL])
			sv_idx = els[SV_IDX_COL]
			event = els[EVENT_COL]

			sv = [chrom, cluster_begin, cluster_end, sv_idx, event]
			if index not in sv_cluster_dict:
				sv_cluster_dict[index] = [sv]
			else:
				sv_cluster_dict[index].append(sv)
			
	part_cluster_f.close()
	return(sv_cluster_dict)




def get_sv_sample(SVs_to_cluster):
	sv_pns = [SVs_to_cluster[i][SAMPLE_COL] for i in range(len(SVs_to_cluster))]
	pn_numdict = {}
	ind = 0
	
	sv_pn_nums = []
	for i in range(len(sv_pns)):
		pn = sv_pns[i]
		if pn not in pn_numdict:
			pn_numdict[pn] = ind
			ind += 1
		sv_pn_nums.append(pn_numdict[pn])


	return(sv_pn_nums, pn_numdict)



def find_interval_overlap(x1, x2, y1, y2):

	overlaps = ((x1 <= y2) and (y1 <= x2))
	x_overlap_rate = 0
	y_overlap_rate = 0

	if overlaps:
		x_len = x2 - x1 + 1
		y_len = y2 - y1 + 1
		begin_overlap = max(x1, y1)
		end_overlap = min(x2, y2)
		overlap_amount = end_overlap - begin_overlap + 1
		x_overlap_rate = float(overlap_amount)/x_len
		y_overlap_rate = float(overlap_amount)/y_len

	return( (overlaps, x_overlap_rate, y_overlap_rate) )

def get_interval_distance(sv1_begin, sv1_end, sv2_begin, sv2_end, sv_type):
	(overlaps, x_overlap_rate, y_overlap_rate) = find_interval_overlap(sv1_begin, sv1_end, sv2_begin, sv2_end)
	if overlaps:
		min_overlap = min(x_overlap_rate, y_overlap_rate)
		max_overlap_left = 1.0 - min_overlap
		return(max_overlap_left)
	else:
		return(1.0)

def form_simplified_distance_matrix(SVs_to_cluster, sv_type):
	distance_matrix = [ [-1 for i in range(len(SVs_to_cluster))] for k in range(len(SVs_to_cluster))]

	for i in range(0, len(SVs_to_cluster)):
		for j in range(i, len(SVs_to_cluster)):
			if i == j:
				distance_matrix[i][j] = 0
			else:	
				distance_matrix[i][j] = get_interval_distance(SVs_to_cluster[i][0], SVs_to_cluster[i][1], SVs_to_cluster[j][0], SVs_to_cluster[j][1], sv_type)
				distance_matrix[j][i] = distance_matrix[i][j]
	return(distance_matrix)

def gather_minimal_distance_matrix(SVs_to_cluster, sv_type):
	unique_sv_dict = {}
	index = 0
	uniq_sv_map = []
	unique_sv_coords = []
	for i in range(len(SVs_to_cluster)):
		sv = SVs_to_cluster[i]
		key = (sv[CLUSTER_BEGIN_COL], sv[CLUSTER_END_COL])
		if key not in unique_sv_dict:
			unique_sv_dict[key] = index
			unique_sv_coords.append(key)
			index += 1			
		uniq_sv_map.append(unique_sv_dict[key])

	# uniq_sv_map is ready, sv locations are stored. Can generate the distance matrix
	minimal_distance_matrix = form_simplified_distance_matrix(unique_sv_coords, sv_type)
	return(minimal_distance_matrix, uniq_sv_map)


