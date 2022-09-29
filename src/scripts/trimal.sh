#!/usr/bin/bash

mkdir out/cds/aln_local_cds/aln_local_cds_trimal
mkdir out/promoter/aln_local_pr/aln_local_promoter_trimal
#mkdir out/aln_promoter_hm/aln_promoter_hm_trimal

#trimal version trimAl 1.2rev59

#cds
cat out/cds/clusters_cds_new/list.cds | parallel --gnu -j 30 \
'~contrera/soft/trimal1.2/source/trimal -in  out/cds/aln_local_cds/{} -out  out/cds/aln_local_cds/aln_local_cds_trimal/{} \
-automated1' ::: &> out/logs/trimal_cds_local_log.txt


#promoter
cat out/promoter/promoter_cluster/list.pr | parallel --gnu -j 30 \
'~contrera/soft/trimal1.2/source/trimal -in out/promoter/aln_local_pr/{} -out out/promoter/aln_local_pr/aln_local_promoter_trimal/{} \
-automated1' ::: &> out/logs/trimal_pr_local_log.txt

#promoter hard mask
#cat out/promoter_cluster_hm/list.pr.hm | parallel --gnu -j 30 \
#'~contrera/soft/trimal1.2/source/trimal -in out/aln_promoter_hm/{} \
#-out  out/aln_promoter_hm/aln_promoter_hm_trimal/{} -automated1' ::: &>  out/logs/trimal_prhm_log.txt
