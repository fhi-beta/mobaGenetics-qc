# Psam file reconstruction
- Number of samples in the genotyping data: 229618.
## Samples not in Medical Birth Regsitry
795 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 3725 |
| Parent-offspring| 138835 |
| Full siblings| 26009 |
| 2nd degree| 143 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](mod8_psam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 25188 |
| Male | 87 |
| Female | 61789 |

![](mod8_psam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 15 |
| Male | 59858 |
| Female | 69 |

![](mod8_psam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 9927 |
| Male | 42554 |
| Female | 30742 |

![](mod8_psam_reconstruction/children_sex_plot.png)
## Parental relationships
795 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
225667 individuals in total. Breakdown excluding multiple same-sex parents:
 -  79184 children
 -  66872 mothers
 -  46185 fathers
 -  77788 mother-child pairs
 -  55043 father-child pairs
 -  53647 mother-father-child trios
 -  33573 unrelated

Multiple same-sex parents (at the individual level):
 -  234 children with more than one mother detected
 -  74 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

77948 mother-child relationships expected.
- 77608 (99.56%) recovered by genetic relationships.
- 340 (0.44%) not recovered by genetic relationships.


52990 father-child relationships expected.
- 52691 (99.44%) recovered by genetic relationships.
- 299 (0.56%) not recovered by genetic relationships.


78022 mother-child relationships detected.
- 77608 (99.47%) matched to registry.
- 414 (0.53%) not matched to registry.


55117 father-child relationships detected.
- 52691 (95.6%) matched to registry.
- 2426 (4.4%) not matched to registry.


###  Samples
229618 samples in total. Breakdown excluding multiple same-sex parents:
 -  80058 children
 -  66914 mothers
 -  46247 fathers
 -  78612 mother-child pairs
 -  55877 father-child pairs
 -  54431 mother-father-child trios
 -  36587 unrelated

Multiple same-sex parents (at the sample level):
 -  2170 children with more than one mother detected
 -  2078 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

80612 mother-child relationships expected.
- 80258 (99.56%) recovered by genetic relationships.
- 354 (0.44%) not recovered by genetic relationships.


55682 father-child relationships expected.
- 55354 (99.41%) recovered by genetic relationships.
- 328 (0.59%) not recovered by genetic relationships.


80782 mother-child relationships detected.
- 80258 (99.35%) matched to registry.
- 524 (0.65%) not matched to registry.


57955 father-child relationships detected.
- 55354 (95.51%) matched to registry.
- 2601 (4.49%) not matched to registry.


## Exclusion
- Number of samples excluded: 919
