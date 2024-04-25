# Fam file reconstruction in snp015a
## Samples not in Medical Birth Regsitry
5595 samples with missing birth year, will be assumed to be parent.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 2 |
| Parent-offspring| 4930 |
| Full siblings| 169 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 24 |
| Male | 1 |
| Female | 3882 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 3905 |
| Female | 1 |

![](fam_reconstruction/father_sex_plot.png)
## Parental relationship
2587 mother-child relationships expected.
- 2575 (99.54%) recovered by genetic relationships.
- 12 (0.46%) not recovered by genetic relationships.
2382 father-child relationships expected.
- 2365 (99.29%) recovered by genetic relationships.
- 17 (0.71%) not recovered by genetic relationships.
4949 parent-offspring relationships detected
- 4940 (99.82%) match to registry.
- 9 (0.18%) do not match to registry.
## Exclusion
- Number of samples excluded: 42
