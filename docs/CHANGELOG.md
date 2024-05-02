# Changelog


Versions are of the form MAJOR.MINOR.PATCH and this changelog tries to conform
to [Common Changelog](https://common-changelog.org)


## unreleased

### Changed

- Handle groups in hit extraction view. ([#84](https://github.com/metagenlab/zDB/pull/84)) (Niklaus Johner)
- Allow using groups to define phenotype in GWAS view. ([#82](https://github.com/metagenlab/zDB/pull/82)) (Niklaus Johner)
- Display form validation errors next to the corresponding fields. ([#83](https://github.com/metagenlab/zDB/pull/83)) (Niklaus Johner)
- Filter VF hits by SeqID and coverage and keep one hit per locus. ([#77](https://github.com/metagenlab/zDB/pull/77)) (Niklaus Johner)
- Improve layout for various views, making better use of available space. ([#70](https://github.com/metagenlab/zDB/pull/70)) (Niklaus Johner)
- Configure repository to publish the docs to [readthedocs](https://zdb.readthedocs.io)
- Port nextflow pipelines to DSL2 language. ([#26](https://github.com/metagenlab/zDB/pull/26)) (Niklaus Johner)
- Add CHANGELOG file. ([#25](https://github.com/metagenlab/zDB/pull/25)) (Niklaus Johner)
- Imported metagenlab_lib code into zdb (https://github.com/metagenlab/zDB/pull/33) (Bastian Marquis)
- Removed eutils module. ([#47](https://github.com/metagenlab/zDB/pull/47)) (Niklaus Johner)
- Show xaxis labels (as links) in heatmap view. ([#59](https://github.com/metagenlab/zDB/pull/59)) (Niklaus Johner)
- Remove unused html templates. ([#61](https://github.com/metagenlab/zDB/pull/61)) (Niklaus Johner)


### Added

- Allow defining groups of genomes in input file. ([#82](https://github.com/metagenlab/zDB/pull/82)) (Niklaus Johner)
- Add view to produce custom plots (phylogenetic trees and table). ([#78](https://github.com/metagenlab/zDB/pull/78)) (Niklaus Johner)
- Add AMRs and VFs to search index. ([#73](https://github.com/metagenlab/zDB/pull/73)) (Niklaus Johner)
- Add VFs to locus and orthogroup views. ([#71](https://github.com/metagenlab/zDB/pull/71)) (Niklaus Johner)
- Add new statistics to home view. ([#70](https://github.com/metagenlab/zDB/pull/70)) (Niklaus Johner)
- Add GWAS views. The associations are calculated using Scoary2. ([#62](https://github.com/metagenlab/zDB/pull/62)) (Niklaus Johner)
- Add VF details view. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add VF accumulation/rarefaction plot view. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add VF heatmap view. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add VF tabular comparison views. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add VF Venn diagram view. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add VF hit extraction view. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add VF comparison index view and link it in the menu and home page. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add VF entry list view. ([#58](https://github.com/metagenlab/zDB/pull/58)) (Niklaus Johner)
- Add virulence factor annotations. ([#56](https://github.com/metagenlab/zDB/pull/56)) (Niklaus Johner)
- Store reference DB and software versions in DB. ([#54](https://github.com/metagenlab/zDB/pull/54)) (Niklaus Johner)
- Add script to compute rendering times of various views. ([#51](https://github.com/metagenlab/zDB/pull/51)) (Niklaus Johner)
- Add possibility to dump HTML responses from webapp tests. ([#50](https://github.com/metagenlab/zDB/pull/50)) (Niklaus Johner)
- Add AMR accumulation/rarefaction plot view. ([#46](https://github.com/metagenlab/zDB/pull/46)) (Niklaus Johner)
- Add AMR heatmap view. ([#45](https://github.com/metagenlab/zDB/pull/45)) (Niklaus Johner)
- Add AMR details view. ([#44](https://github.com/metagenlab/zDB/pull/44)) (Niklaus Johner)
- Add AMR Venn diagram view. ([#41](https://github.com/metagenlab/zDB/pull/41)) (Niklaus Johner)
- Add AMR hit extraction view. ([#38](https://github.com/metagenlab/zDB/pull/38)) (Niklaus Johner)
- Add AMR entry list view. ([#37](https://github.com/metagenlab/zDB/pull/37)) (Niklaus Johner)
- Add AMR comparison index view and link it in the menu and home page. ([#35](https://github.com/metagenlab/zDB/pull/35)) (Niklaus Johner)
- Add AMR tabular comparison views. ([#24](https://github.com/metagenlab/zDB/pull/24), [#31](https://github.com/metagenlab/zDB/pull/31)) (Niklaus Johner)
- Add integration tests for pipelines and views. ([#29](https://github.com/metagenlab/zDB/pull/29)) (Niklaus Johner)
- Add unitesting for pipelines. ([#26](https://github.com/metagenlab/zDB/pull/26)) (Niklaus Johner)
- Add antimicrobial resistance annotations. ([#23](https://github.com/metagenlab/zDB/pull/23)) (Niklaus Johner)
- Add support for circular contigs in the genomic viewer. (https://github.com/metagenlab/zDB/pull/33) (Bastian Marquis)

### Fixed

- Fix displaying active tab in search results view. ([#73](https://github.com/metagenlab/zDB/pull/73)) (Niklaus Johner)
- Fix missing page title for kegg_mapp_ko view. ([#66](https://github.com/metagenlab/zDB/pull/66)) (Niklaus Johner)
- Handle empty result set in heatmap, venn and pan_genome views. ([#66](https://github.com/metagenlab/zDB/pull/66)) (Niklaus Johner)
- Fix tabular comparison views navigation tabs. ([#30](https://github.com/metagenlab/zDB/pull/30), [#31](https://github.com/metagenlab/zDB/pull/31)) (Niklaus Johner)
- Fix the representation of genes overflowing from a genomic region, that would appear shorter than they actually are. ([#52](https://github.com/metagenlab/zDB/pull/52)) (Bastian Marquis)
- Fix missing page title on error in blast view. ([#64](https://github.com/metagenlab/zDB/pull/64)) (Niklaus Johner)


## 1.2.1 - 2023-10-16

### Changed
- Bump checkm-genome version to 1.2.2. (Trestan Pillonel)


## 1.2.0 - 2023-10-13

### Changed

- Use mamba to build conda envs ([30bbb7530d75c3adf986aaf4c5052ad83ea301b8](https://github.com/metagenlab/zDB/commit/30bbb7530d75c3adf986aaf4c5052ad83ea301b8)) (Trestan Pillonel)

### Fixed

- Fix issue for genome file download ([#16](https://github.com/metagenlab/zDB/pull/16)) (Trestan Pillonel)
- Fix file extensions in annotation pipeline ([#17](https://github.com/metagenlab/zDB/pull/17)) (Thomas Kozusnik)


## 1.1.1 - 2023-05-16

### Fixed

- Fixed several issues with conda environments. (Bastian Marquis and Trestan Pillonel)


## 1.1.0 - 2023-03-30

### Added

- Add ```--resume``` option to ```zdb setup``` ([272b4b434279d725b5860db26458c430934fddd5](https://github.com/metagenlab/zDB/commit/272b4b434279d725b5860db26458c430934fddd5)) (Bastian Marquis)


## 1.0.8 - 2023-03-15

### Added
- Add support for docker and conda. ([#13](https://github.com/metagenlab/zDB/pull/13)) (Bastian Marquis)

### Fixed
- Several bugfixes


## 1.0.5 - 2022-12-08

### Added

- Add ```--resume``` option to ```zdb run``` ([#11](https://github.com/metagenlab/zDB/pull/11)) (Bastian Marquis)
- Add name field to the input csv file.
- Accept genbank files with other extensions in the pipeline.


## 1.0.4 - 2022-09-21

### Fixed

- Fix conda build ([#10](https://github.com/metagenlab/zDB/pull/10)) (Bastian Marquis)


## 1.0.3 - 2022-06-07

### Fixed

- Correctly set number of missing orthogroups. ([b9df07c2668f4d70cea2a345a7b3fdf2f6d55bfb](https://github.com/metagenlab/zDB/commit/b9df07c2668f4d70cea2a345a7b3fdf2f6d55bfb)) (Bastian Marquis)
- Fix bug with singularity. ([4d37463ed78789f26352c56d1c2af9c50e2b09b9](https://github.com/metagenlab/zDB/commit/4d37463ed78789f26352c56d1c2af9c50e2b09b9)) (Bastian Marquis)


## 1.0.2 - 2022-05-10

### Added

- Add support for conda. ([#9](https://github.com/metagenlab/zDB/pull/9)) (Bastian Marquis)


## 1.0.0 - 2022-04-01

_Initial release._

