# Fam file reconstruction in snp009
## Samples not in Medical Birth Regsitry
40 samples with missing birth year, will be assumed to be parent.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 5 |
| Parent-offspring| 10321 |
| Full siblings| 340 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 5676 |
| Male | 7 |
| Female | 9 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 5 |
| Male | 5377 |
| Female | 0 |

![](fam_reconstruction/father_sex_plot.png)
## Parental relationship
5443 mother-child relationships expected.
- 5422 (99.61%) recovered by genetic relationships.
- 21 (0.39%) not recovered by genetic relationships.


5202 father-child relationships expected.
- 5175 (99.48%) recovered by genetic relationships.
- 27 (0.52%) not recovered by genetic relationships.


10614 parent-offspring relationships detected
- 10597 (99.84%) match to registry.
- 17 (0.16%) do not match to registry.


## Exclusion
- Number of samples excluded: 63
