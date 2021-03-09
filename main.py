import os
import sys
import measure_SV_distance as SVdist
import CAST_weighted as CC_w



##### Merge the data using the CAST algorithm
def cast_cluster_pieced_data(part_cluster_fn, event, overlap_percentage, merged_out_fn):
	overlap_distance_thresh = (1 - (overlap_percentage/float(100)) )

	sv_cluster_dict = SVdist.read_pieced_data(part_cluster_fn)

	keys = sorted(sv_cluster_dict.keys())
	
	for k in keys:
		#print(k)
		clique_ids = 0
		minimal_distance_matrix, dm_map = SVdist.gather_minimal_distance_matrix(sv_cluster_dict[k], event)
		sv_pn_nums, pn_numdict = SVdist.get_sv_sample(sv_cluster_dict[k])
		list_of_cliques = CC_w.find_cliques(dm_map, sv_pn_nums, minimal_distance_matrix, overlap_distance_thresh)
		#print(sv_cluster_dict[k])
		#print(list_of_cliques)
		for cli in list_of_cliques:
			for sv_el in cli:
				sv_id = sv_cluster_dict[k][sv_el][3]
				print('\t'.join([sv_id, k+'.'+str(clique_ids)]))
			clique_ids += 1		
		
if __name__ == '__main__':
	part_cluster_fn = sys.argv[1]
	event = sys.argv[2] 
	overlap_percentage = float(sys.argv[3])
	merged_out_fn = sys.argv[4]
	cast_cluster_pieced_data(part_cluster_fn, event, overlap_percentage, merged_out_fn)