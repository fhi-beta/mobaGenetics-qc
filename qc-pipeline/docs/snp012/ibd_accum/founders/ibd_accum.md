# Kinship filtering in snp012
Relatedness filtering, {length(ids)} individuals.
## Relatedness
![](ibd_accum/kinship_plot.png)
## Relatedness
- Pairwise kinship
![](ibd_accum/kinship_density.png)
- Cummulative positive kinship
![](ibd_accum/cumulated_kinship_density.png)
Percolation of the relatedness graph using a Kinship threshold of {kinship_threshold}: {length(excluded_ids)} excluded, {length(unrelated_ids)} remaining.
## Relatedness after relatedness filtering
- Pairwise kinship
![](ibd_accum/kinship_density_unrelated.png)
- Cummulative positive kinship
![](ibd_accum/cumulated_kinship_density_unrelated.png)
Removal of samples with accumulated kinship using threshold of {accumulated_kinship_threshold}: {length(excluded_cumulative_kinship)} excluded, {length(retained_ids)} remaining.
