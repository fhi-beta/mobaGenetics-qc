# Fam file reconstruction in snp017d
## Samples not in Medical Birth Regsitry
2290 samples with missing birth year, will be assumed to be parent.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 9 |
| Parent-offspring| 460 |
| Full siblings| 26 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 22 |
| Male | 1 |
| Female | 2148 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 1110 |
| Female | 2 |

![](fam_reconstruction/father_sex_plot.png)
## Parental relationship
262 mother-child relationships expected.
- 260 (99.24%) recovered by genetic relationships.
- 2 (0.76%) not recovered by genetic relationships.
200 father-child relationships expected.
- 200 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.
461 parent-offspring relationships detected
- 460 (99.78%) match to registry.
- 1 (0.22%) do not match to registry.
## Exclusion
- Number of samples excluded: 5
