#!/usr/bin/bash

#creat directory
mkdir out/promoter/aln_promoter/pr_identity
mkdir out/cds/aln_cds/cds_identity
#promoter
cat out/cds/clusters_cds_new/list.cds | parallel --gnu -j 30 './src/scripts/similarity_seq_ratio.py  -f out/promoter/aln_promoter/aln_promoter_trimal/{} -o out/promoter/aln_promoter/pr_identity/{}_id.txt' ::: &> out/logs/identity_pr_log.txt

#cds
cat out/cds/clusters_cds_new/list.cds | parallel --gnu -j 30 './src/scripts/similarity_seq_ratio.py  -f out/cds/aln_cds/aln_cds_trimal/{} -o out/cds/aln_cds/cds_identity/{}_id.txt' ::: &> out/logs/identity_cds_log.txt
