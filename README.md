# sv-merger

# Sample execution
python main.py ./test_data/toy_SV_data_tomerge_simplified_ids.csv DEL 50

The input CSV file contains the structural variants (SV) to be merged.
The columns are as follows: 

1: chromosome

2: begin site

3: end site (use begin site + insertion length for insertions)

4: SV id (unique identifier for the SVs)

5: Sample id (unique identifier for the sample the given SV is found in)

6: Method/algorithm finding the SV

7: SV type (e.g. DEL, INS)

8*: Pre-clustering id. (unique identifier for a preclustering of SVs)

* A pre-clustering is needed to reduce time complexity in merging the SVs. Two SVs that has no chance to exist in a set of merged SVs should belong to different "pre-clusters". 

# Output 

The columns are as follows:

1: SV ids

2: The final clique id for given SV id in the 1st column.

For each clique id, a representative SV can be chosen if there are more than 1 SV per clique id. One approach would be to pick an SV with the most frequent begin, end, or begin-and-end coordinate among the SVs within the same clique id.
