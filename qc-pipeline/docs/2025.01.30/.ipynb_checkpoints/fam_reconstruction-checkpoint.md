# Psam file reconstruction
- Number of samples in the genotyping data: 229606.
## Samples not in Medical Birth Regsitry
795 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 3758 |
| Parent-offspring| 139070 |
| Full siblings| 25971 |
| 2nd degree| 1 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 4 |
| Male | 0 |
| Female | 86815 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 59781 |
| Female | 55 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 3 |
| Male | 42421 |
| Female | 40527 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
795 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
225657 individuals in total. Breakdown excluding multiple same-sex parents:
 -  79254 children
 -  66963 mothers
 -  46225 fathers
 -  77912 mother-child pairs
 -  55119 father-child pairs
 -  53777 mother-father-child trios
 -  33363 unrelated

Multiple same-sex parents (at the individual level):
 -  237 children with more than one mother detected
 -  74 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

77948 mother-child relationships expected.
- 77733 (99.72%) recovered by genetic relationships.
- 215 (0.28%) not recovered by genetic relationships.


52983 father-child relationships expected.
- 52756 (99.57%) recovered by genetic relationships.
- 227 (0.43%) not recovered by genetic relationships.


78149 mother-child relationships detected.
- 77733 (99.47%) matched to registry.
- 416 (0.53%) not matched to registry.


55193 father-child relationships detected.
- 52756 (95.58%) matched to registry.
- 2437 (4.42%) not matched to registry.


###  Samples
229606 samples in total. Breakdown excluding multiple same-sex parents:
 -  80130 children
 -  67006 mothers
 -  46282 fathers
 -  78743 mother-child pairs
 -  55953 father-child pairs
 -  54566 mother-father-child trios
 -  36377 unrelated

Multiple same-sex parents (at the sample level):
 -  2183 children with more than one mother detected
 -  2096 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

80612 mother-child relationships expected.
- 80395 (99.73%) recovered by genetic relationships.
- 217 (0.27%) not recovered by genetic relationships.


55674 father-child relationships expected.
- 55436 (99.57%) recovered by genetic relationships.
- 238 (0.43%) not recovered by genetic relationships.


80926 mother-child relationships detected.
- 80395 (99.34%) matched to registry.
- 531 (0.66%) not matched to registry.


58049 father-child relationships detected.
- 55436 (95.5%) matched to registry.
- 2613 (4.5%) not matched to registry.


## Exclusion
- Number of samples excluded: 629
