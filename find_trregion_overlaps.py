import os
import sys

#import bisect
from intervaltree import Interval, IntervalTree



def read_regions(trf_regions_fn):

	tr_regs = []

	trf_regions_f = open(trf_regions_fn, 'r')
	for l in trf_regions_f:
		begin, end, _, idx = map(int, l.strip().split('\t')[1:])
		tr_regs.append((begin, end, idx))

	trf_regions_f.close()
	return(tr_regs)


def find_trregions(svs_fn, trf_regions_fn, out_fn, relaxation_bp=5):

	# get the tandem repeat regions.
	tr_regs = read_regions(trf_regions_fn)

	# add said regions into an interval tree.
	t = IntervalTree()
	for r in tr_regs:
		#t[r[0]:r[1]] = (r[2], r[3])
		t[r[0]:r[1]] = r[2]


	out_f = open(out_fn, 'w')
	svs_f = open(svs_fn, 'r')
	
	counter = 0
	for l in svs_f:
		if l[0] != '#':
			els = l.strip().split('\t')
			sv_begin = int(els[1]) 
			sv_end = int(els[2])
			sv_idx = els[3]

			
			# First find intersection regions, among them: 
			# find those that engulf the SV with some relaxation
			sv_tr_regions = []
			for i in sorted(t[sv_begin:sv_end]):
				if (sv_begin >= (i.begin - relaxation_bp) ) and (sv_end <=  (i.end + relaxation_bp)):
					#sv_tr_regions.append((i.begin, i.end, i.data[0], i.data[1]) )
					sv_tr_regions.append((i.begin, i.end, i.data) )


			# Find the largest such region
			if (len(sv_tr_regions) == 0):
				#tr_reg_data = ['.', '.', '.', '.']
				tr_reg_data = ['.', '.', '.']

			else:
				largest_sv_tr_region_idx = sorted(sv_tr_regions, key=lambda x : (x[1]-x[0]), reverse=True)[0][-1]
				tr_reg_data = tr_regs[largest_sv_tr_region_idx]

			write_line = l.strip() + '\t' + '\t'.join(map(str, tr_reg_data))
			out_f.write(write_line + '\n')

			counter += 1
			if counter % 100000 == 0: 
				print counter


	svs_f.close()
	out_f.close()

def main():

	svs_fn = sys.argv[1]
	trf_regions_fn = sys.argv[2]
	out_fn = sys.argv[3]
	relaxation_bp = int(sys.argv[4])
	find_trregions(svs_fn, trf_regions_fn, out_fn, relaxation_bp)


if __name__ == '__main__':
	main()

