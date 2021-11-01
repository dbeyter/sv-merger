import os
import sys

#import functools

from collections import namedtuple as nt


# On the TRF combined csv.
CHR = 0
BEGIN = 1
END = 2
SV_IDX = 3
SVLEN = 7

TRR_BEGIN = 8
TRR_END = 9
TRR_IDX = 10




# TODO: Switch them to the usage within the data structure to sort the SVs.
CHR_COL = 0
CLUSTER_BEGIN_COL = 1
CLUSTER_END_COL = 2
SV_IDX_COL = 3
SVLEN_COL = 4

AUX_BEGIN_COL = 5
AUX_END_COL = 6
TRR_IDX_COL = 7
CLUSTERID_COL = 8





def combine_sv_info_with_trf_intersections(trf_intersections_fn, out_fn, overlap_p, use_trf_coordinates=True):
	overlap_p = overlap_p/float(100)

	trf_coords, sv_trf_indices = parse_trf_containments(trf_intersections_fn)

	out_f = open(out_fn, 'w')

	all_vars = [ [ [], 0, 0] ]
	svs_f = open(trf_intersections_fn, 'r')

	counter = 0
	for l in svs_f:
		if l[0] != '#':
			if counter % 50000 == 0:
				print(str(counter) + " svs read.")
			counter += 1

			#chrom, sv_begin, sv_end, sv_idx, event, sv_len = l.strip().split('\t')
			els = l.strip().split('\t')
			chrom = els[CHR]
			sv_begin = int(els[BEGIN])
			sv_end = int(els[END])
			sv_idx = els[SV_IDX]
			sv_len = int(els[SVLEN])
			

			if sv_idx in sv_trf_indices:
				trf_begin_coord = int(trf_coords[sv_trf_indices[sv_idx]][0])
				trf_end_coord = trf_begin_coord + sv_len
				trf_index = sv_trf_indices[sv_idx]
			else:
				trf_begin_coord = sv_begin
				trf_end_coord = sv_end
				trf_index = "."

			if use_trf_coordinates:
				sv_els = [ chrom, trf_begin_coord, trf_end_coord, sv_idx, sv_len, sv_begin, sv_end, trf_index, "0"]
			else:
				sv_els = [ chrom, sv_begin, sv_end, sv_idx, sv_len, trf_begin_coord, trf_end_coord, trf_index, "0"]

			all_vars[0][0].append(sv_els)
	
	svs_f.close()

	print('All svs are read. sorting.')
	done_sv_clusters_list = separate_clusters(all_vars, overlap_p)
	# all_vars.sort(key=lambda x: (x[CHR_COL], x[CLUSTER_BEGIN_COL], x[CLUSTER_END_COL]))
	print('Initial Sort ended. writing.')

	print('ClusterID and pos sorting starting.')
	flat_sv_clusters_list = []
	for cluster in done_sv_clusters_list:
		for sv in cluster[0]:
			flat_sv_clusters_list.append(sv)
	#print(flat_sv_clusters_list[0])
	flat_sv_clusters_list.sort(key=lambda x: map(int, x[CLUSTERID_COL].split('.')) + [x[CLUSTER_BEGIN_COL]]  )
	
	# for sv in flat_sv_clusters_list:
	# 	write_line = '\t'.join(map(str, sv))
	# 	out_f.write(write_line + '\n')


	# Arrange the SVs
	arranging_svs = {}
	for i in range(len(flat_sv_clusters_list)):
		sv = flat_sv_clusters_list[i]
		clusterID = sv[CLUSTERID_COL]
		#print clusterID
		if clusterID not in arranging_svs:
			arranging_svs[clusterID] = [int(sv[CLUSTER_BEGIN_COL]), [(sv[CLUSTER_BEGIN_COL],i)] ] 
		else:
			if sv[CLUSTER_BEGIN_COL] < arranging_svs[clusterID][0]:
				arranging_svs[clusterID][0] = sv[CLUSTER_BEGIN_COL]
			arranging_svs[clusterID][1].append((sv[CLUSTER_BEGIN_COL],i))
	
	list_of_cluster_to_print = []
	for cID in arranging_svs:
		arranging_svs[cID][1].sort()
		list_of_cluster_to_print.append(arranging_svs[cID])
	list_of_cluster_to_print.sort()
	print('sorting ended.')

	#counter = 0

	extra_sv_data = {}
	svs_f = open(trf_intersections_fn, 'r')
	for l in svs_f:
		if l[0] != '#':
			els = l.strip().split('\t')
			extra_sv_data[els[3]] = "\t".join(els[4:7])
	# done.

	for i in range(len(list_of_cluster_to_print)):
		begin_coord = list_of_cluster_to_print[i][0]
		sv_indices = list_of_cluster_to_print[i][1]
		for si in range(len(sv_indices)):
			#print svs[s]
			coord, sv_idx = sv_indices[si]
			#write_line = '\t'.join(map(str, flat_sv_clusters_list[sv_idx]))
			sv_data = flat_sv_clusters_list[sv_idx]
			write_line = '\t'.join(map(str, sv_data[0:4] + [extra_sv_data[sv_data[3]]] + sv_data[4:]))

			out_f.write(write_line + '\n')

	print('Write ended.')

	
	out_f.close()


# SV clusters are presented as [], 0, 0. The first 0 means that it has gone through POS separating, and the next 0 means that it has gone through SVLEN separating.
# If a cluster has gone through both separations, then we can end it. 
def separate_clusters(list_of_clusters, overlap_p):

	done_sv_clusters_list = []
	while len(list_of_clusters) > 0:
		todo_by_POS(done_sv_clusters_list, list_of_clusters)
		todo_by_SVLEN(done_sv_clusters_list, list_of_clusters, overlap_p)

	return(done_sv_clusters_list)



#input all vars as a list of lists. Let's say all SVs has a clusterID.
def todo_by_POS(done_sv_clusters_list, undone_list_of_clusters):

	print(len(undone_list_of_clusters))
	if len(undone_list_of_clusters) > 0:
		# Check if the sv list to pop has been already operated on todoPOS
		if undone_list_of_clusters[-1][1] != 1:
			G = undone_list_of_clusters.pop()

			G[0].sort(key=lambda x: (x[CHR_COL], x[CLUSTER_BEGIN_COL], x[CLUSTER_END_COL]))
			pos_separated_G = []

			SVs_to_cluster = [[], 1, 0]
			last_endpos = G[0][0][CLUSTER_END_COL]
			G[0].append(['K'] +[float('Inf') for x in range(11)])
			cluster_id = 0
			for i in range(len(G[0])):
				cur_SV = G[0][i]
				#print cur_SV[CLUSTER_BEGIN_COL], last_endpos
				if cur_SV[CLUSTER_BEGIN_COL] < last_endpos:

					cur_SV[CLUSTERID_COL] = cur_SV[CLUSTERID_COL] + "." + str(cluster_id)
					SVs_to_cluster[0].append(cur_SV)
					#print 'Appending to cluster'
				else:
					#print("LEN and Last el of the SVs_to_cluster appended: " + str(len(SVs_to_cluster[0])) + " and " + str(SVs_to_cluster[0][-1]) )
					pos_separated_G.append((SVs_to_cluster))
					#print "LAST POS sep group before ending: " + str(pos_separated_G[-1])
					# if len(pos_separated_G)> 3:
					# 	print "secondtoLAST POS sep group before ending: " + str(pos_separated_G[-2])
					cluster_id += 1
					if cur_SV[CLUSTERID_COL]  != float('Inf'):
						cur_SV[CLUSTERID_COL] = cur_SV[CLUSTERID_COL] + "." + str(cluster_id)
					SVs_to_cluster = [[cur_SV], 1, 0]

				if last_endpos < cur_SV[CLUSTER_END_COL]:
					last_endpos = cur_SV[CLUSTER_END_COL]

			#print "POS sep group number: " + str(len(pos_separated_G))
			# for i in range(1,5):
			# 	#print "pos_separated_G: " + str(pos_separated_G[0-(i)])
			# 	print "pos_separated_G: " + str(len(pos_separated_G[0-(i)][0]))		

			if len(pos_separated_G) > 1:
				# append for further processing.
				for sv_cluster in pos_separated_G:
					undone_list_of_clusters.append(sv_cluster)
			elif len(pos_separated_G) == 1:
				if G[2] == 1:
					done_sv_clusters_list.append(pos_separated_G[0])
				else:
					undone_list_of_clusters.append(pos_separated_G[0])
			else:
				print("pos_separated_G should not be zero. Abort")
				raise(Exception)
	# 	else:
	# 		print "Sending to LEN"
	# else:
	# 	print " Undone list of clusters is done."
	return


def svlen_overlap(svlen1, svlen2):
	p = min(svlen1, svlen2)/float(max(svlen1, svlen2))
	return(p)

def todo_by_SVLEN(done_sv_clusters_list, undone_list_of_clusters, overlap_p):

	# Check if the sv list to pop has been already operated on todoPOS
	if len(undone_list_of_clusters) > 0:
		if undone_list_of_clusters[-1][2] != 1:	

			G = undone_list_of_clusters.pop()
			G[0].sort(key=lambda x: (x[CHR_COL], x[SVLEN_COL]))
			svlen_separated_G = []

			SVs_to_cluster = [[], 0, 1]
			prev_SV = G[0][0]
			G[0].append(['K'] +[float('Inf') for x in range(11)])
			cluster_id = 0
			for i in range(len(G[0])):
				cur_SV = G[0][i]
				svoverlap = svlen_overlap(cur_SV[SVLEN_COL], prev_SV[SVLEN_COL])
				# print "prev_SV: " + str(prev_SV)
				# print "cur_SV: " + str(cur_SV)
				# print "overlap: " + str(overlap_p)
				if svoverlap >= overlap_p:

					cur_SV[CLUSTERID_COL] = cur_SV[CLUSTERID_COL] +"." + str(cluster_id)
					#print cur_SV
					SVs_to_cluster[0].append(cur_SV)
				else:
					#prev_SV[CLUSTERID_COL] = prev_SV[CLUSTERID_COL] +"." + str(cluster_id)
					svlen_separated_G.append(SVs_to_cluster)
					cluster_id += 1

					if cur_SV[CLUSTERID_COL]  != float('Inf'):
						cur_SV[CLUSTERID_COL] = cur_SV[CLUSTERID_COL] +"." + str(cluster_id)
					SVs_to_cluster = [[cur_SV], 0, 1]
				prev_SV = cur_SV

			#print "LEN sep group number: " + str(len(svlen_separated_G))
			#print svlen_separated_G
			#raise(Exception)
			if len(svlen_separated_G) > 1:
			# append for further processing.
				for sv_cluster in svlen_separated_G:
					undone_list_of_clusters.append(sv_cluster)
			elif len(svlen_separated_G) == 1:
				if G[1] == 1:
					done_sv_clusters_list.append(svlen_separated_G[0])
				else:
					undone_list_of_clusters.append(svlen_separated_G[0])
			else:
				print("pos_separated_G should not be zero. Abort")
				raise(Exception)
	# 	else:
	# 		print "Sending to POS"
	# else:
	# 	print " Undone list of clusters is done."

	return



def parse_trf_containments(trf_intersections_fn):
	# Read the TRF intersection fn
	trf_intersections_f = open(trf_intersections_fn, 'r')
	trf_coords = {}
	sv_trf_indices ={}
	for l in trf_intersections_f:
		if l[0] != '#':
			els = l.strip().split('\t')
			sv_idx = els[SV_IDX]
			trf_idx = els[TRR_IDX]
			if trf_idx != ".":
				trf_begin = int(els[TRR_BEGIN])
				trf_end = int(els[TRR_END])
				
				if trf_idx not in trf_coords:
					trf_coords[trf_idx] = (trf_begin, trf_end)
				sv_trf_indices[sv_idx] = trf_idx
	trf_intersections_f.close()
	return(trf_coords, sv_trf_indices)

def main():

	sv_trr_intersections_fn = sys.argv[1]
	out_fn = sys.argv[2]
	overlap_percentage = int(sys.argv[3])
	use_trf_coordinates = bool(int(sys.argv[4]))
	combine_sv_info_with_trf_intersections(sv_trr_intersections_fn, out_fn, overlap_percentage, use_trf_coordinates)


if __name__ == '__main__':
	main()