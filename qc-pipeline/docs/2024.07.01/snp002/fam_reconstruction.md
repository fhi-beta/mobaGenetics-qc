# Fam file reconstruction in snp002
## Samples not in Medical Birth Regsitry
5 samples with missing birth year, will be assumed to be parent.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 10 |
| Parent-offspring| 262 |
| Full siblings| 3 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 2 |
| Female | 468 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 462 |
| Female | 0 |

![](fam_reconstruction/father_sex_plot.png)
## Parental relationship
182 mother-child relationships expected.
- 181 (99.45%) recovered by genetic relationships.
- 1 (0.55%) not recovered by genetic relationships.


110 father-child relationships expected.
- 110 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


291 parent-offspring relationships detected
- 291 (100%) match to registry.
- 0 (0%) do not match to registry.


## Exclusion
- Number of samples excluded: 3
