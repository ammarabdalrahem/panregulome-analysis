#to create tables/total_seq_cluster_pr.csv for calculation poor alignment 
for file in out/uncom_data/promoter_cluster/*; do echo -n "$(basename "$file"),"; grep ">" "$file" | wc -l; done > tables/total_seq_cluster_pr.csv

#to create tables/total_seq_cluster_cds.csv for calculation poor alignment 
for file in out/uncom_data/cds_cluster/*; do echo -n "$(basename "$file"),"; grep ">" "$file" | wc -l; done > tables/total_seq_cluster_cds.csv
