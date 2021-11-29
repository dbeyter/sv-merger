import os
import sys
import measure_SV_distance as SVdist
import CAST_weighted as CC_w

import precluster_svs as pcl
import find_trregion_overlaps as tro




### Separate the file into SVs within TRR and outside TRRs.
def separate_sv_inout_trr(part_cluster_fn):

	part_cluster_f = open(part_cluster_fn, 'r')

	intrr_fn = part_cluster_fn + '.intrr'
	outtrr_fn = part_cluster_fn + '.outtrr'
	
	part_cluster_intrr_f = open(intrr_fn, 'w')
	part_cluster_outtrr_f = open(outtrr_fn, 'w')
	for l in part_cluster_f:
		if l.split('\t')[pcl.TRR_IDX] != '.':
			part_cluster_intrr_f.write(l)
		else:
			part_cluster_outtrr_f.write(l)
	
	part_cluster_outtrr_f.close()
	part_cluster_intrr_f.close()
	return(intrr_fn, outtrr_fn)



##### Merge the data using the CAST algorithm
def cast_cluster_pieced_data(part_cluster_fn, sv_merge_out_fn, event, overlap_percentage):
	overlap_distance_thresh = (1 - (overlap_percentage/float(100)) )

	sv_cluster_dict = SVdist.read_pieced_data(part_cluster_fn)

	keys = sorted(sv_cluster_dict.keys())
	

	sv_merge_out_f = open(sv_merge_out_fn, 'w')
	for k in keys:
		clique_ids = 0
		minimal_distance_matrix, dm_map = SVdist.gather_minimal_distance_matrix(sv_cluster_dict[k], event)
		sv_pn_nums, pn_numdict = SVdist.get_sv_sample(sv_cluster_dict[k])
		list_of_cliques = CC_w.find_cliques(dm_map, sv_pn_nums, minimal_distance_matrix, overlap_distance_thresh)
		
		for cli in list_of_cliques:
			for sv_el in cli:
				sv_id = sv_cluster_dict[k][sv_el][3]
				#print('\t'.join([sv_id, k+'.'+str(clique_ids)]))
				sv_merge_out_f.write('\t'.join([sv_id, k+'.'+str(clique_ids)]) + '\n')
			clique_ids += 1

	sv_merge_out_f.close()

if __name__ == '__main__':


	operation = sys.argv[1]

	if operation == "MERGE":
		svs_fn = sys.argv[2]
		trf_regions_fn = sys.argv[3]
		trr_overlap_out_fn = svs_fn + '.trr_overlap'
		relaxation_bp = 5
		
		# Find TRR vs non-TRR SVs.
		tro.find_trregions(svs_fn, trf_regions_fn, trr_overlap_out_fn, relaxation_bp)

		precluster_out_fn = svs_fn + '.precluster'
		overlap_percentage = 50
		use_trf_coordinates = 1
		
		# Pre-cluster SVs.
		pcl.combine_sv_info_with_trf_intersections(trr_overlap_out_fn, precluster_out_fn, overlap_percentage, use_trf_coordinates)

		# Separate in vs. out trr SVs. 
		intrr_fn, outtrr_fn = separate_sv_inout_trr(precluster_out_fn)

		event = sys.argv[4]
		intrr_overlap_percentage = 85
		intrr_sv_merge_out_fn = svs_fn + '.intrr.merged.csv'
		cast_cluster_pieced_data(intrr_fn, intrr_sv_merge_out_fn, event, intrr_overlap_percentage)
		print("intrr done.")

		outtrr_sv_merge_out_fn = svs_fn + '.outtrr.merged.csv'
		cast_cluster_pieced_data(outtrr_fn, outtrr_sv_merge_out_fn, event, overlap_percentage)
		print("outtrr done.")


	elif operation == "FIND_TRR_OVERLAPS":
		svs_fn = sys.argv[2]
		trf_regions_fn = sys.argv[3]
		trr_overlap_out_fn = sys.argv[4]
		relaxation_bp = int(sys.argv[5])
		
		# Find TRR vs non-TRR SVs.
		tro.find_trregions(svs_fn, trf_regions_fn, trr_overlap_out_fn, relaxation_bp)

	elif operation == "PRE_CLUSTER":
		trr_overlap_out_fn = sys.argv[2]
		precluster_out_fn = sys.argv[3]
		overlap_percentage = int(sys.argv[4])
		use_trf_coordinates = bool(int(sys.argv[5]))
		
		# Pre-cluster SVs.
		pcl.combine_sv_info_with_trf_intersections(trr_overlap_out_fn, precluster_out_fn, overlap_percentage, use_trf_coordinates)

	elif operation == "FIND_CLIQUES":
		preclustered_svs_fn = sys.argv[2]
		intrr_sv_merge_out_fn = sys.argv[3]
		outtrr_sv_merge_out_fn = sys.argv[4]
		
		event = sys.argv[5]
		intrr_overlap_percentage = int(sys.argv[6])
		outtrr_overlap_percentage = int(sys.argv[7])

		# Separate in vs. out trr SVs. 
		intrr_fn, outtrr_fn = separate_sv_inout_trr(preclustered_svs_fn)
		
		cast_cluster_pieced_data(intrr_fn, intrr_sv_merge_out_fn, event, intrr_overlap_percentage)
		cast_cluster_pieced_data(outtrr_fn, outtrr_sv_merge_out_fn, event, outtrr_overlap_percentage)
