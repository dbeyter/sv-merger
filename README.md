# sv-merger

# Requirements:
Must have python 2 or 3.

Additional python packages: intervaltree

Please install via: `pip install intervaltree`

# Sample execution
arg 1: Tab separated file containing the structural variants (SV) to be merged.

arg 2: File containing tandem repeat coordinates given under folder trf_coords.

arg 3: SV type, e.g. DEL, INS.
```
python main.py MERGE ./test_data/toy_SV_data.csv ./trf_coords/chr21.trf.sorted.gor DEL
```

Note: Merges SVs within and outside separately. Uses an overlap threshold of 50% for SVs outside of tandem repeats, and a threshold of 85% for SVs within the tandem repeats.

# In-depth execution
### Finding SVs within and outside tandem repeat regions (TRR)s.
arg 0: Name of the step to execute.

arg 1: Tab separated file containing the structural variants (SV) to be merged.

arg 2: File containing tandem repeat coordinates given under folder trf_coords.

arg 3: Output file name for this step.

arg 4: Relaxation parameter while finding SVs within/outside tandem repeats. E.g. When `(SV.begin >= TRR.begin - relaxation) and (SV.end <= TRR.end + relaxation)`, where SV and TRR denotes a structural variant and tandem repeat region, respectively, the SV is accepted as within the given TRR.

```
python main.py FIND_TRR_OVERLAPS ./test_data/toy_SV_data.csv ./trf_coords/chr21.trf.sorted.gor ./test_data/toy_SV_data.csv.trr_overlap  5
```
### Pre-clustering SVs.
arg 0: Name of the step to execute.

arg 1: File name containing SV and TRR overlaps, i.e. the output file name from the previous step.

arg 2: Output file name for this step.

arg 3: Merging overlap percentage. Note: If different merging overlap parameters will be used for SVs within and outside TRRs, use the smaller percentage in this step.

arg 4: Boolean flag for using the tandem repeat coordinates in SV pre-clustering and merging. Using 1 will carry the SVs within TRRs to the start site of their respective TRRs, using 0 will use the original SV sites.
 
```
python main.py PRE_CLUSTER ./test_data/toy_SV_data.csv.trr_overlap  ./test_data/toy_SV_data.csv.precluster 50 1
```

### Merging SVs.
arg 0: Name of the step to execute.

arg 1: File name containing pre-clustered SVs, i.e. output file name from the previous step.

arg 2: Output file name for merged SVs within TRRs.

arg 3: Output file name for merged SVs outside TRRs.

arg 4: SV type, e.g. DEL, INS.

arg 5: Overlap percentage for SVs within TRRs.

arg 6: Overlap percentage for SVs outside TRRs.

```
python main.py FIND_CLIQUES ./test_data/toy_SV_data.csv.precluster ./test_data/toy_SV_data.csv.intrr.merged.csv ./test_data/toy_SV_data.csv.outtrr.merged.csv DEL 85 50
```

# Columns for the input file containing the SVs

0: chromosome

1: begin site

2: end site (use begin site + insertion length for insertions)

3: SV id (unique identifier for the SVs)

4: Sample id (unique identifier for the sample the given SV is found in)

5: Method/algorithm finding the SV

6: SV type (e.g. DEL, INS)

7: SV length

# Output 

The columns are as follows:

1: SV ids

2: The final clique id for given SV id in the 1st column.

For each clique id, a representative SV can be chosen if there are more than 1 SV per clique id. One approach would be to pick an SV with the most frequent begin, end, or begin-and-end coordinate among the SVs within the same clique id.



