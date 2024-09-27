[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## MoBa Genotypes QC
This reposirtory contains code and resources used for the quality control of genotyping data from the [Norwegian Mother, Father and Child Cohort Study (MoBa)](https://www.fhi.no/en/studies/moba/).

### Current release
The current release is available for alpha testing for expert users. We expect a new release to be publicly available for beta testing in the first semester of 2025. For more information on the release process, please contact the release coordinator [Dana Kristjansson](https://www.fhi.no/en/ab/departments-and-centres/genetics-and-bioinformatics/dana-hashim).

👉 [Release notes and documentation](qc-pipeline/docs/release_notes.md): Documentation on the current and previous releases.

### General information on the genotypes

#### Batch effects
Samples in MoBa have been genoptyped in different batches on different genotyping platforms. Variants are not excluded based on batch effects, we recommend to account for batch effects in your own analyses. A table listing the batch for each sample is available with the data release. Information on batches is available at the [mobagen repository](https://github.com/folkehelseinstituttet/mobagen/wiki/MoBaGenetics1.5).

#### Content
We conducted only minimal filtering to allow you to tailor the release to your project: the release includes genotypes with imputation quality > 0.3 and all samples including duplicates and outliers in the PCA. Accompanying tables are available with variant and sample information to help you filter the data for your specific application. Genotypes are available for autosomal and sex chromosomes, please note that variant coverage varies between batches, especially for the sex chromosomes. The genotypes are provided as phased dosages in the plink2 format, a single file containing genotypes of parents and children. For details on the format please refer to the [plink documentation](https://www.cog-genomics.org/plink/2.0/formats#pgen).

👉 For more information on the data processing, please go to the [QC pipeline documentation](qc-pipeline/readme.md). 

#### Build and imputation
Currently genotypes are only available on build 37 imputed against the [Haplotype Reference Consortium Release 1.1](https://ega-archive.org/datasets/EGAD00001002729) (11,227 samples).


### Other pipelines and releases
The current pipeline and its associated releases have been built on previous efforts:
- [MoBa PsychGen](github.com/psychgen/MoBaPsychGen-QC-pipeline): A quality-controlled set of genotypes created by the [PsychGen Centre for Genetic Epidemiology and Mental Health](https://www.fhi.no/en/me/the-psychgen-centre-for-genetic-epidemiology-and-mental-health).
- [MoBaGenetics v.1.0](https://github.com/folkehelseinstituttet/mobagen/wiki/MoBaGenetics1.0): A quality-controlled set of genotypes obtained after a first round of genotyping created by the [Genetics and Bioinformatics department of the Norwegian Institute of Public Health](https://www.fhi.no/om/organisasjon/genetikk-og-bioinformatikk/).

### Citation
When using code or data produced using this pipeline, please refer to this repository and cite the reference work by the PsychGen team [doi: 10.1101/2022.06.23.496289](https://doi.org/10.1101/2022.06.23.496289).

### Content
- [GWAS](gwas/readme.md): this folder contains code and resources used to run a test GWAS.
- [Phenotypes](phenotypes/readme.md): this folder contains code and resources for the handling of phenotypes used in the pipeline.
- [QC Pipeline](qc-pipeline/readme.md): this folder contains code, resources, and documentation for the QC pipeline
- [Quarantine](quarantine/readme.md): this folder contains deprecated executables and resources.

### Troubleshooting and support
If despite our efforts at enforcing good practices in our work, you encounter issues or have suggestions of improvement, please use our [issue tracker](https://github.com/fhi-beta/mobaGenetics-qc/issues).
For technical questions, feedback, and assistance regarding the release please contact the [Genetics and Bioinformatics department of the Norwegian Institute of Public Health](https://www.fhi.no/om/organisasjon/genetikk-og-bioinformatikk/). 
To reach out to the community working on genetics in MoBa, please join our [Slack](https://join.slack.com/t/mobagen/shared_invite/zt-2r90hzo82-HNllFHDugSxJeknpJ9jT0w).

### Credits
This repository is maintained by the [Genetics and Bioinformatics department of the Norwegian Institute of Public Health](https://www.fhi.no/om/organisasjon/genetikk-og-bioinformatikk/), and notably Øystein Kapperud, Dana Kristjansson, Even Birkeland, and Marc Vaudel. The department is deeply grateful for the contribution of its former members Gutorm Thomas Høgåsen and Ronny Myhre.   

The pipeline is built on work for the [MoBaGenetics v.1.0](https://github.com/folkehelseinstituttet/mobagen/wiki/MoBaGenetics1.0) release, led by Øyvind Helgeland, Julius Juodakis, Jonas Bacelis, and Ronny Myhre.

It was extended thanks to help from the [PsychGen team](https://www.fhi.no/en/me/the-psychgen-centre-for-genetic-epidemiology-and-mental-health), especially Elizabeth Claire Corfield.

The pre-phasing and imputation modules have been built following the implementation of the [NORMENT team](https://www.med.uio.no/norment/english/) and the [PsychGen team](https://www.fhi.no/en/me/the-psychgen-centre-for-genetic-epidemiology-and-mental-health), especially Oleksandr Frei and Tahir Tekin Filiz.

We are grateful for the contribution of the community who helped, shared code, and provided support and feedback. We are especially grateful to the [PsychGen team](https://www.fhi.no/en/me/the-psychgen-centre-for-genetic-epidemiology-and-mental-health), the [NORMENT team](https://www.med.uio.no/norment/english/), the [Centre for Fertility and Health](https://www.fhi.no/en/ch/Centre-for-fertility-and-health/), the [Department of Obstetrics and Gynecology at the University of Gothenburg](https://www.gu.se/en/about/find-organisation/department-of-obstetrics-and-gynecology), and the [Mohn Center for Diabetes Precision Medicine](https://www.uib.no/en/diabetes).

The data are processed in the secure digital laboratories of the [HUNT Cloud at the Norwegian University of Science and Technology](https://www.ntnu.edu/mh/huntcloud). We are grateful for outstanding support from the HUNT Cloud community.

### License
The code and resources generated for the QC pipeline is freely available under the permissive [GNU General Public License v3.0](LICENSE). Other licenses might apply for external components of the pipeline. No copyright infringement intended.

### Code of Conduct
As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Convenant](https://www.contributor-covenant.org/) [Code of Conduct for Open Source Projects](CODE_OF_CONDUCT.md).
