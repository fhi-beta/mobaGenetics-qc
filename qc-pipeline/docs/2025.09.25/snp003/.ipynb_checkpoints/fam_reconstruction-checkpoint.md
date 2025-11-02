# Fam file reconstruction in snp003
- Number of samples in the genotyping data: 12716.
## Samples not in Medical Birth Regsitry
63 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 138 |
| Parent-offspring| 8641 |
| Full siblings| 277 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 14 |
| Male | 2 |
| Female | 4112 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 4053 |
| Female | 0 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 3 |
| Male | 2301 |
| Female | 2231 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
63 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
12515 individuals in total. Breakdown excluding multiple same-sex parents:
 -  4247 children
 -  4006 mothers
 -  3973 fathers
 -  4221 mother-child pairs
 -  4187 father-child pairs
 -  4161 mother-father-child trios
 -  289 unrelated

Multiple same-sex parents (at the individual level):
 -  0 children with more than one mother detected
 -  2 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

4223 mother-child relationships expected.
- 4221 (99.95%) recovered by genetic relationships.
- 2 (0.05%) not recovered by genetic relationships.


4187 father-child relationships expected.
- 4187 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


4221 mother-child relationships detected.
- 4221 (100%) matched to registry.
- 0 (0%) not matched to registry.


4189 father-child relationships detected.
- 4187 (99.95%) matched to registry.
- 2 (0.05%) not matched to registry.


###  Samples
12716 samples in total. Breakdown excluding multiple same-sex parents:
 -  4302 children
 -  4006 mothers
 -  3974 fathers
 -  4275 mother-child pairs
 -  4241 father-child pairs
 -  4214 mother-father-child trios
 -  434 unrelated

Multiple same-sex parents (at the sample level):
 -  87 children with more than one mother detected
 -  30 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

4364 mother-child relationships expected.
- 4362 (99.95%) recovered by genetic relationships.
- 2 (0.05%) not recovered by genetic relationships.


4268 father-child relationships expected.
- 4268 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


4362 mother-child relationships detected.
- 4362 (100%) matched to registry.
- 0 (0%) not matched to registry.


4271 father-child relationships detected.
- 4268 (99.93%) matched to registry.
- 3 (0.07%) not matched to registry.


## Exclusion
- Number of samples excluded: 18
