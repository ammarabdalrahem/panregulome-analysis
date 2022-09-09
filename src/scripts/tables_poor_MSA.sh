for i in $(cat out/promoter/promoter_cluster/list.pr ); do echo -n $i ;echo -n ",";cat out/promoter/promoter_cluster/$i |grep ">"|wc -l ;done > tables/total_seq_cluster_pr_.csv

