# Fam file reconstruction in snp016a
- Number of samples in the genotyping data: 24762.
## Samples not in Medical Birth Regsitry
96 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 36 |
| Parent-offspring| 3194 |
| Full siblings| 327 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 9482 |
| Male | 20 |
| Female | 3002 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 4 |
| Male | 3903 |
| Female | 1 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 3212 |
| Male | 4091 |
| Female | 1047 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
96 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
24666 individuals in total. Breakdown excluding multiple same-sex parents:
 -  2949 children
 -  2695 mothers
 -  393 fathers
 -  2774 mother-child pairs
 -  406 father-child pairs
 -  231 mother-father-child trios
 -  18629 unrelated

Multiple same-sex parents (at the individual level):
 -  7 children with more than one mother detected
 -  0 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

2791 mother-child relationships expected.
- 2764 (99.03%) recovered by genetic relationships.
- 27 (0.97%) not recovered by genetic relationships.


370 father-child relationships expected.
- 366 (98.92%) recovered by genetic relationships.
- 4 (1.08%) not recovered by genetic relationships.


2781 mother-child relationships detected.
- 2764 (99.39%) matched to registry.
- 17 (0.61%) not matched to registry.


406 father-child relationships detected.
- 366 (90.15%) matched to registry.
- 40 (9.85%) not matched to registry.


###  Samples
24762 samples in total. Breakdown excluding multiple same-sex parents:
 -  2949 children
 -  2695 mothers
 -  394 fathers
 -  2774 mother-child pairs
 -  407 father-child pairs
 -  232 mother-father-child trios
 -  18724 unrelated

Multiple same-sex parents (at the sample level):
 -  7 children with more than one mother detected
 -  0 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

2791 mother-child relationships expected.
- 2764 (99.03%) recovered by genetic relationships.
- 27 (0.97%) not recovered by genetic relationships.


370 father-child relationships expected.
- 366 (98.92%) recovered by genetic relationships.
- 4 (1.08%) not recovered by genetic relationships.


2781 mother-child relationships detected.
- 2764 (99.39%) matched to registry.
- 17 (0.61%) not matched to registry.


407 father-child relationships detected.
- 366 (89.93%) matched to registry.
- 41 (10.07%) not matched to registry.


## Exclusion
- Number of samples excluded: 77
