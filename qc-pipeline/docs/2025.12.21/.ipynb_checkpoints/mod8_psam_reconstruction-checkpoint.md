# Psam file reconstruction
- Number of samples in the genotyping data: 229606.
## Samples not in Medical Birth Regsitry
795 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 3758 |
| Parent-offspring| 139070 |
| Full siblings| 25973 |
| 2nd degree| 3 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](mod8_psam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 0 |
| Female | 86819 |

![](mod8_psam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 59781 |
| Female | 55 |

![](mod8_psam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 42421 |
| Female | 40530 |

![](mod8_psam_reconstruction/children_sex_plot.png)
## Parental relationships
795 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
225657 individuals in total. Breakdown excluding multiple same-sex parents:
 -  79255 children
 -  66962 mothers
 -  46226 fathers
 -  77911 mother-child pairs
 -  55120 father-child pairs
 -  53776 mother-father-child trios
 -  33362 unrelated

Multiple same-sex parents (at the individual level):
 -  237 children with more than one mother detected
 -  74 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

77948 mother-child relationships expected.
- 77732 (99.72%) recovered by genetic relationships.
- 216 (0.28%) not recovered by genetic relationships.


52983 father-child relationships expected.
- 52757 (99.57%) recovered by genetic relationships.
- 226 (0.43%) not recovered by genetic relationships.


78148 mother-child relationships detected.
- 77732 (99.47%) matched to registry.
- 416 (0.53%) not matched to registry.


55194 father-child relationships detected.
- 52757 (95.58%) matched to registry.
- 2437 (4.42%) not matched to registry.


###  Samples
229606 samples in total. Breakdown excluding multiple same-sex parents:
 -  80099 children
 -  67002 mothers
 -  46272 fathers
 -  78745 mother-child pairs
 -  55948 father-child pairs
 -  54594 mother-father-child trios
 -  36379 unrelated

Multiple same-sex parents (at the sample level):
 -  2184 children with more than one mother detected
 -  2098 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

80612 mother-child relationships expected.
- 80394 (99.73%) recovered by genetic relationships.
- 218 (0.27%) not recovered by genetic relationships.


55674 father-child relationships expected.
- 55437 (99.57%) recovered by genetic relationships.
- 237 (0.43%) not recovered by genetic relationships.


80929 mother-child relationships detected.
- 80394 (99.34%) matched to registry.
- 535 (0.66%) not matched to registry.


58046 father-child relationships detected.
- 55437 (95.51%) matched to registry.
- 2609 (4.49%) not matched to registry.


## Exclusion
- Number of samples excluded: 630
