# sv-merger

# Sample execution
python main.py ./test_data/toy_SV_data_tomerge_simplified_ids.csv DEL 50

# Arguments

arg 1: Tab separated file containing the structural variants (SV) to be merged.

arg 2: SV type, e.g. DEL, INS.

arg 3: Minimum percentage of overlap to draw an edge between two SVs in clique formation during SV merging.

# Columns for the input file containing the SVs

0: chromosome

1: begin site

2: end site (use begin site + insertion length for insertions)

3: SV id (unique identifier for the SVs)

4: Sample id (unique identifier for the sample the given SV is found in)

5: Method/algorithm finding the SV

6: SV type (e.g. DEL, INS)

7*: Pre-clustering id. (unique identifier for a preclustering of SVs)

*A pre-clustering is needed to reduce time complexity in merging of SVs. Two SVs that are guaranteed not to be merged together should belong to different "pre-clusters". 

# Output 

The columns are as follows:

1: SV ids

2: The final clique id for given SV id in the 1st column.

For each clique id, a representative SV can be chosen if there are more than 1 SV per clique id. One approach would be to pick an SV with the most frequent begin, end, or begin-and-end coordinate among the SVs within the same clique id.
