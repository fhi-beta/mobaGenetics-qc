

<style>
.column-left{
  float: left;
  width: 50%;
}
.column-right{
  float: right;
  width: 50%;
}
</style>
**Module 0 - Data conversion**
------------------------------

------------------------------------------------------------------------

Data conversion from GenomeStudio to PLINK format was carried out prior to data delivery.

**Module 1 - Data preparation**
-------------------------------

------------------------------------------------------------------------

Module 1 prepared data for QC procedures. In the exported data from GenomeStudio, samples were coded using project specific retrieval IDs and non-informative family IDs. Also, no sex information was available.

Information about declared sample sex and pedigree structure provided by MoBa were used to update the dataset.

Markers with poor cluster separation, low 10% GC score and high AA theta deviation were removed. Clustering metrics were provided by the SNP table exported from GenomeStudio.

The array contains some duplicated/triplicated markers. Duplicates/triplicates were removed to avoid potential problems downstream.

Illumina provides a conversion list for converting marker IDs used on the array to rsIDs. In the provided list some markers are named '.'. To avoid duplicate/non-informative IDs, the original chip ID was used in situations were no rsID was available.

Illumina technical markers (assigned to chromosome 0) were excluded and markers with poor clustering were eliminated using metrics available from the variant table exported from GenomeStudio. Problematic markers reported by other consortia were subsequently excluded.

NOTE: All markers in the PAR region was correctly assigned to chr 24, no update necessary.
NOTE: The lab in Rotterdam removed 7710 markers showing poor performance prior to delivery (markers not included in the dataset)

**Updating sample IDs:**
Samples in: 17949
Samples updated: 17949
Samples not updated: 0

**Updating parental information:**
Samples in: 17949
Samples assigned one or more parents: 5984
Samples not updated: 11965

**Updating sex:**
Samples in: 17949
Samples where sex was updated: 17949
Samples not updated for sex: 0

**Remove markers by cluster separation &lt; 0.4:**
Markers in: 692367
Markers removed: 18154
Markers remaining: 674213

**Remove markers by 10% GC score:**
Markers in: 674213
Markers removed: 19915
Markers remaining: 654298

**Remove markerst by AA theta dev:**
Markers in: 654298
Markers removed: 4123
Markers remaining: 650175

**Remove duplicated markers:**
Markers in: 650175
Markers removed: 480
Markers remaining: 649695

**Update SNPs to rsID:**
Markers in: 649695
Markers updated: 338017
Markers not updated: 311678

**Removed technical markers (chr0):**
Markers in: 649695
Markers removed: 1910
Markers remaining: 647785

**Module 2 - Identify core samples and infere pedigree**
--------------------------------------------------------

------------------------------------------------------------------------

### Predefined (static) input files

The input file *A* with pre-solved problematic pedigrees contained **7** resolved families with **20** index individuals. The input file *B* with pre-identified problematic samples (often accidental duplicates of other samples) contained **65** individual IDs.

### Upstream (dynamic) input files

Genetic files with **17742** individuals reached this module. The original .fam file listed **5815** fathers (V3 column) and **5832** mothers (V4 column) who had genotypes (i.e., were listed in V2 column); also **66** fathers and **58** mothers without genotypes. The final .fam file will have them reset as missing (respective numbers **0** and **0**).

### Thresholds and procedures for relationship and sex inference

The thresholds used to identify paren-offspring relationship were *Z*1&gt; 0.8; twin or dublicated samples - *P**I*\_*H**A**T*≥ 0.8; full-siblings relationship - *Z*1≥ 0.35, *Z*1≤ 0.65, *P**I*\_*H**A**T*≥ 0.35, *P**I*\_*H**A**T*≤ 0.65. The Y chromosome genotype count threshold used to separate males from females was *Y**C*&gt; 92. The X chromosome *F* threshold used to separate females from males was *F*&lt; 0.648. Genetic sex was inferred based on both criteria. When criteria disagreed (**4** cases), samples were flagged as not suitable for analyses *(phenotypeOK=FALSE)* and genetic sex was infered from the X chromosome data.

### Modifications to .fam file

Being found in the input file *B*, **20** samples in the .fam file got assigned their true family IDs, genetic parents and genetic sex. These samples are suitable for analyses and are not flagged as problematic. Being found in the input file *A*, **64** samples in the .fam file got assigned dummy family IDs (e.g. "prblm001"), got founder's status (i.e., parental IDs were set to "0") and their declared sex was set to their genetic sex. They were flagged as not suitable for future analyses *(phenotypeOK=FALSE)*. In sex inference for all these updated samples, there were **0** cases where Y-chromosome and X-chromosome data did not agree (likely Klinefelter). The declared and inferred sex did not match in **35** samples.

The remaining .fam file contained **11549** declared parent-offspring relationships (**5756** paternal, **5793** maternal). Genetic inferrence of the same data detected **11549** parent-offspring relationships and **0** pairs of dublicated (twin) samples. If there is a difference between declared and inferred numbers, the auto-generated .pdf report should be manually inspected to detect new sample-identity problems. In sex inference for all these samples, there were **4** cases where Y-chromosome and X-chromosome data did not agree (likely Klinefelter), and **4** cases where declared and inferred sex did not agree. These samples were flagged as not suitable for analyses *(phenotypeOK=FALSE)*.

### Summary

In total, the number of samples flagged as not suitable for analyses *(phenotypeOK=FALSE)* is **68**. We do not trust the identity of these samples. The remaining **17674** samples were flagged as OK *(phenotypeOK=TRUE)*.

The updated .fam file contained **11562** declared parent-offspring relationships (**5762** paternal, **5800** maternal). The genetic (Xchr) sex was assigned to all the samples.

### Problematic genotyping arrays

The samples from the input file *A* where enriched in these genotyping arrays:

|      | chip            |  prob\_count|
|------|:----------------|------------:|
| 2    | 201641480150    |            8|
| 12   | 201680470069    |            8|
| 24   | 201689730161    |            7|
| 1    | 201570320148    |            1|
| 3    | 201641510161    |            1|
| 4    | 201641770011    |            1|
| (the | table is trimme |           d)|

The samples from the input file *B* where enriched in these genotyping arrays:

|      | chip            |  prob\_count|
|------|:----------------|------------:|
| 11   | 201692730188    |            3|
| 1    | 201629310100    |            2|
| 2    | 201629310157    |            2|
| 5    | 201680470065    |            2|
| 3    | 201641780148    |            1|
| 4    | 201641790166    |            1|
| (the | table is trimme |           d)|

In module 2 an ethnically homogeneous core set of samples were identifed for use in module 3 (marker cleaning). Marker cleaning requires an ethnically homogeneous sample set in order to facilitate marker cleaning sensitive to ethnic outliers. Initially markers with MAF &gt; 10% were temporarily removed. Samples with missingness rate &gt; 5% were permanently removed. Markers with missingness &gt; 2% were temporarily removed. The resulting dataset was then used to assess and update the sex of samples where reported and genetic sex did not match.

Further, strand-ambiguous markers, non-autosomal markers, and markers with HWE p &lt; 1e-4 were temporarily removed before IBD estimation.

the dataset was subsequently split into a parent dataset and offspring dataset, the latter including only one child in case of siblings (selected at random). The resulting datasets were cleaned separately in the modules described in the following modules and merged in module 7.

**Markers and samples in**
Markers in start: 647785
Samples in start: 17949

**Temporary removal of markers with MAF &lt; 10%**
Markers in: 647785
Markers removed: 415378
Markers remaining: 232407

**Permanent removal of samples with missingness &gt; 5%**
Samples in: 17949
Samples removed: 207
Samples remaining 17742

**Temporary removal of markers with missingness &gt; 2%**
Markers removed: 232407
Markers removed: 3879
Markers remaining: 228528

**Temporary removal of non-autosomal markers**
Markers in: 228528

Markers removed: 8521
Markers remaining: 220007

**Temporary removal of markers with HWE p &lt; 1e-4**
Markers in: 220007
Markers removed: 842
Markers remaining: 219165

**Temporary removal of strand ambiguous markers**
Markers removed: 219165
Markers removed: 1205
Markers remaining: 217960

**Temporary remove markers in high LD**
Markers in: 217960
Markers removed: 6338
Markers remaining: 211622

**Prune set of markers using --indep-pairwise 200 100 0.1**
Markers in: 211622
Markers removed: 171386
Markers remaining: 40236

JONAS STUFF...

**Removal of pedigree inconsistent samples**
Samples in: 17742
Samples for removal: 64
Samples removed: 64
Samples remaining: 17678

**PCA after merge with HapMap**
Markers after HapMap merge (used for PCA): 20308

**Sample selection post PCA**
Samples in: 17678
Samples removed after PCA: 785
Samples remaining after PCA: 16893

**Split dataset into founders and offspring**
Samples in: 16893
Founders: 11316
Offspring: 5577

Founders
--------

**IBD estimation**
Samples in: 11316
Markers in: 647785

**Remove samples with excess accumulated PIHAT:**
Samples in: 11316
Samples for removal: 14
Samples removed: 14
Samples remaining: 11302

**Remove one in a pair of samples with PIHAT &gt; 0.1:**
Samples in: 11302
Samples for removal: 478
Samples removed: 461
Samples remaining: 10841

Offspring
---------

**IBD estimation**
Samples in: 5577
Markers in: 647785

**Remove samples with excess accumulated PIHAT:**
Samples in: 5577
Samples for removal: 11
Samples removed: 11
Samples remaining: 5566

**Remove one in a pair of samples with PIHAT &gt; 0.1:**
Samples in: 5566
Samples for removal: 90
Samples removed: 87
Samples remaining: 5479

**Module 3 - Identify good markers**
------------------------------------

------------------------------------------------------------------------

Founders
--------

**Number of markers and samples at start of cleaning:**
Samples in start: 10841
Markers in start: 647785

**Remove markers with missingness &gt; 10%:**
Markers in: 647785
Markers removed: 288
Markers remaining: 647497

**Remove individuals with missingsness &gt; 5%:**
Samples in: 10841
Samples removed: 1
Samples remaining: 10840

**Remove markers with missingness &gt; 5%:**
Markers in: 647497
Markers removed: 1417
Markers remaining: 646080

**Remove individuals with missingess &gt; 3%:**
Samples in: 10840
Samples removed: 19
Samples remaining: 10821

**Remove markers with missingness &gt; 2%:**
Markers in: 646080
Markers removed: 8322
Markers remaining: 637758

**Remove individuals with missingness &gt; 2%:**
Samples in: 10821
Samples removed: 18
Samples remaining: 10803

**Remove autosomal markers with HWE p &lt; 1e-7:**
Markers in: 637758
Markers removed: 1659
Markers remaining: 636099

**Remove samples with HET excess &gt; 4SD using common autosomal markers (MAF &gt; 0.01)**
Samples in: 10803
Samples removed: 1
Samples remaining: 10802

**Remove autosomal markers with HWE p &lt; 1e-6**
Markers in: 636099
Markers removed: 208
Markers remaining: 635891

**Remove samples with HET excess &gt; 4SD using rare autosomal markers (MAF &gt; 0.01)**
Samples in: 10802
Samples removed: 61
Samples remaining: 10741

**Remove markers with missingness &gt; 2%:**
Markers in: 635891
Markers removed: 3
Markers remaining: 635888

**Temporarily remove samples failing sex check (F: 0.2, 0.8):**
Samples in: 10741
Samples for removal: 5
Samples removed: 5
Samples out: 10736

**Markers into sex clean:**
X markers in: 15099
Y markers in: 712
PAR markers in: 564
MT markers in: 126

**Remove chrX and PAR markers with HWE p &lt; 1e-6 (only female):**
Markers (X + PAR) in: 15663
Markers removed: 29
Markers remaining: 15634

**Remove chrX marker if any male has at least one heterozygote genotype:**
Markers removed: 694
Markers remaining 635165

**Markers after sex clean:**
Autosomes markers out: 619387
X markers out: 14404
Y markers out: 712
PAR markers out: 536
MT markers out: 126
TOTAL: 635165

Offspring
---------

**Number of markers and samples at start of cleaning:**
Samples in start: 5479
Markers in start: 647785

**Remove markers with missingness &gt; 10%:**
Markers in: 647785
Markers removed: 246
Markers remaining: 647539

**Remove individuals with missingsness &gt; 5%:**
Samples in: 5479
Samples removed: 2
Samples remaining: 5477

**Remove markers with missingness &gt; 5%:**
Markers in: 647539
Markers removed: 1223
Markers remaining: 646316

**Remove individuals with missingess &gt; 3%:**
Samples in: 5477
Samples removed: 19
Samples remaining: 5458

**Remove markers with missingness &gt; 2%:**
Markers in: 646316
Markers removed: 7548
Markers remaining: 638768

**Remove individuals with missingness &gt; 2%:**
Samples in: 5458
Samples removed: 9
Samples remaining: 5449

**Remove autosomal markers with HWE p &lt; 1e-7:**
Markers in: 638768
Markers removed: 974
Markers remaining: 637794

**Remove samples with HET excess &gt; 4SD using common autosomal markers (MAF &gt; 0.01)**
Samples in: 5449
Samples removed: 4
Samples remaining: 5445

**Remove autosomal markers with HWE p &lt; 1e-6**
Markers in: 637794
Markers removed: 123
Markers remaining: 637671

**Remove samples with HET excess &gt; 4SD using rare autosomal markers (MAF &gt; 0.01)**
Samples in: 5445
Samples removed: 34
Samples remaining: 5411

**Remove markers with missingness &gt; 2%:**
Markers in: 637671
Markers removed: 5
Markers remaining: 637666

**Temporarily remove samples failing sex check (F: 0.2, 0.8):**
Samples in: 5411
Samples for removal: 2
Samples removed: 2
Samples out: 5409

**Markers into sex clean:**
X markers in: 15251
Y markers in: 712
PAR markers in: 564
MT markers in: 126

**Remove chrX and PAR markers with HWE p &lt; 1e-6 (only female):**
Markers (X + PAR) in: 15815
Markers removed: 25
Markers remaining: 15790

**Remove chrX marker if any male has at least one heterozygote genotype:**
Markers removed: 477
Markers remaining 637164

**Markers after sex clean:**
Autosomes markers out: 621013
X markers out: 14774
Y markers out: 712
PAR markers out: 539
MT markers out: 126
TOTAL: 637164

**Module 4 - Individuals for analyses**
---------------------------------------

------------------------------------------------------------------------

### Founders

**Markers and samples at beginning of module:**
Markers start: 635165
Samples start: 11316

**Remove markers not surviving QC in both parents and offspring:**
Markers in: 635165
Markers removed: 427
Markers remaining: 634738

**Remove samples with missingness rate &gt; 2%:**
Samples in: 11316
Samples removed: 42
Samples remaining: 11274

**Remove samples with HET excess &gt; 4SD using common autosomal markers (MAF &gt; 0.01):**
Samples in: 11274
Samples removed: 1
Samples remaining: 11273

**Remove samples with HET excess &gt; 4SD using rare autosomal markers (MAF &lt; 0.01):**
Samples in: 11273
Samples removed: 65
Samples remaining: 11208

**Remove samples with excess accumulated PIHAT:**
Samples in: removed 11208
Samples removed: 11
Samples remaining: 11197

**Remove one in a pair of samples with PI\_HAT &gt; 0.1:**
Samples in: 11197
Samples removed: 457
Samples remaining: 10740

### Offspring

**Markers and samples at beginning of module:**
Markers start: 637164
Samples start: 5577

**Remove markers not surviving QC in both parents and offspring:**
Markers in: 637164
Markers removed: 2426
Markers remaining: 634738

**Remove samples with missingness rate &gt; 2%:**
Samples in: 5577
Samples removed: 33
Samples remaining: 5544

**Remove samples with HET excess &gt; 4SD using common autosomal markers (MAF &gt; 0.01):**
Samples in: 5544
Samples removed: 4
Samples remaining: 5540

**Remove samples with HET excess &gt; 4SD using rare autosomal markers (MAF &lt; 0.01):**
Samples in: 5540
Samples removed: 38
Samples remaining: 5502

**Remove samples with excess accumulated PIHAT:**
Samples in: removed 5502
Samples removed: 7
Samples remaining: 5495

**Remove one in a pair of samples with PI\_HAT &gt; 0.1:**
Samples in: 5495
Samples removed: 86
Samples remaining: 5409

**Module 5 - Preparation for phasing and imputation**
-----------------------------------------------------

------------------------------------------------------------------------

**Samples and markers into module:**
Samples in: 17742
Markers in: 647785

**Remove markers not passing QC for both offspring and founders:**
Markers in: 647785
Markers shared: 634738
Markers removed: 13047
Markers remaining: 634738

**Remove markers above chr 23:**
Markers in: 634738
Markers removed: 1374
Markers remaining: 633364

**Set mendelian errors to missing:**
Mendelian errors zeroed: 186590

**HRC harmonizing:**
Markers in: 633364
Marker chromosomes changed: 0
Marker positions changed: 0
Marker strand flips: 37109
Marker allele flips: 575573
Markers excluded (not in HRC): 65089
Markers after exclusion: 568275

**Number of markers per chromosome sent to phasing:**
<table class="table" style="width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
Chromosome
</th>
<th style="text-align:right;">
N
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:right;width: 20em; ">
45268
</td>
</tr>
<tr>
<td style="text-align:left;">
2
</td>
<td style="text-align:right;width: 20em; ">
46161
</td>
</tr>
<tr>
<td style="text-align:left;">
3
</td>
<td style="text-align:right;width: 20em; ">
38636
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:right;width: 20em; ">
35311
</td>
</tr>
<tr>
<td style="text-align:left;">
5
</td>
<td style="text-align:right;width: 20em; ">
33424
</td>
</tr>
<tr>
<td style="text-align:left;">
6
</td>
<td style="text-align:right;width: 20em; ">
39844
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:right;width: 20em; ">
31182
</td>
</tr>
<tr>
<td style="text-align:left;">
8
</td>
<td style="text-align:right;width: 20em; ">
29286
</td>
</tr>
<tr>
<td style="text-align:left;">
9
</td>
<td style="text-align:right;width: 20em; ">
24176
</td>
</tr>
<tr>
<td style="text-align:left;">
10
</td>
<td style="text-align:right;width: 20em; ">
28150
</td>
</tr>
<tr>
<td style="text-align:left;">
11
</td>
<td style="text-align:right;width: 20em; ">
28293
</td>
</tr>
<tr>
<td style="text-align:left;">
12
</td>
<td style="text-align:right;width: 20em; ">
27051
</td>
</tr>
<tr>
<td style="text-align:left;">
13
</td>
<td style="text-align:right;width: 20em; ">
19613
</td>
</tr>
<tr>
<td style="text-align:left;">
14
</td>
<td style="text-align:right;width: 20em; ">
18290
</td>
</tr>
<tr>
<td style="text-align:left;">
15
</td>
<td style="text-align:right;width: 20em; ">
17219
</td>
</tr>
<tr>
<td style="text-align:left;">
16
</td>
<td style="text-align:right;width: 20em; ">
18758
</td>
</tr>
<tr>
<td style="text-align:left;">
17
</td>
<td style="text-align:right;width: 20em; ">
17044
</td>
</tr>
<tr>
<td style="text-align:left;">
18
</td>
<td style="text-align:right;width: 20em; ">
16219
</td>
</tr>
<tr>
<td style="text-align:left;">
19
</td>
<td style="text-align:right;width: 20em; ">
13131
</td>
</tr>
<tr>
<td style="text-align:left;">
20
</td>
<td style="text-align:right;width: 20em; ">
13854
</td>
</tr>
<tr>
<td style="text-align:left;">
21
</td>
<td style="text-align:right;width: 20em; ">
7705
</td>
</tr>
<tr>
<td style="text-align:left;">
22
</td>
<td style="text-align:right;width: 20em; ">
8191
</td>
</tr>
<tr>
<td style="text-align:left;">
X
</td>
<td style="text-align:right;width: 20em; ">
11469
</td>
</tr>
</tbody>
</table>
