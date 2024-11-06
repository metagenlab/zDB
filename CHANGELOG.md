# Changelog


Versions are of the form MAJOR.MINOR.PATCH and this changelog tries to conform
to [Common Changelog](https://common-changelog.org)

## 1.3.4 - 2024-11-06


- Fix annotation pipeline for new AMRFinder version 4. ([#141](https://github.com/metagenlab/zDB/pull/141)) (Niklaus Johner)
- Fix DB setup pipeline when using containers and an unmounted reference directory. ([#142](https://github.com/metagenlab/zDB/pull/142)) (Niklaus Johner)

### Changed

- Display Venn diagram as Edwards for 6 genomes, as numbers were missing in the classic representation. ([#135](https://github.com/metagenlab/zDB/pull/135)) (Niklaus Johner)
- Add contig accession and range to genomic regions in plot_region. ([#134](https://github.com/metagenlab/zDB/pull/134)) (Niklaus Johner)
- Improve error handling for missing input genomes and reference databases. ([#130](https://github.com/metagenlab/zDB/pull/130)) (Niklaus Johner)
- Optimize order of regions in plot_region. ([#137](https://github.com/metagenlab/zDB/pull/137)) (Niklaus Johner)
- Improve region edge representations in plot_region. ([#138](https://github.com/metagenlab/zDB/pull/138)) (Niklaus Johner)
- Add contig table to genome overview (extract_contig view). ([#139](https://github.com/metagenlab/zDB/pull/139)) (Niklaus Johner)
- Clean-up css files. ([#140](https://github.com/metagenlab/zDB/pull/140)) (Niklaus Johner)

### Added

- Add contig card to locusx view. ([#136](https://github.com/metagenlab/zDB/pull/136)) (Niklaus Johner)
- Add option to highlight certain leaves in custom plots. ([#131](https://github.com/metagenlab/zDB/pull/131)) (Niklaus Johner)


## 1.3.3 - 2024-08-08

### Fixed

- Fix Kegg based comparative analysis form. ([#123](https://github.com/metagenlab/zDB/pull/123)) (Niklaus Johner)
- Add missing assets when running webapp with nginx and gunicorn. ([#122](https://github.com/metagenlab/zDB/pull/122)) (Niklaus Johner)
- Fix checkm conda environment. ([#121](https://github.com/metagenlab/zDB/pull/121)) (Niklaus Johner)


## 1.3.2 - 2024-07-30

### Fixed

- Fix annotation pipeline for nextflow version > 24.03. ([#120](https://github.com/metagenlab/zDB/pull/120)) (Niklaus Johner)

## 1.3.1 - 2024-07-19

### Fixed

- Fix binding of folders when running with singularity. ([#116](https://github.com/metagenlab/zDB/pull/116)) (Niklaus Johner)
- Fix checkM conda environment. ([#114](https://github.com/metagenlab/zDB/pull/114)
- Fix AMRFinderPlus database update for Linux. ([#113](https://github.com/metagenlab/zDB/pull/113)
- Fix AMRFinderPlus database update for Linux. ([#113](https://github.com/metagenlab/zDB/pull/113)
- Update containers zDB release v 1.3. ([#112](https://github.com/metagenlab/zDB/pull/112)


## 1.3.0 - 2024-07-11 (broken release as containers are outdated)

### Changed

- Update the documentation. ([#110](https://github.com/metagenlab/zDB/pull/110)) (Niklaus Johner)
- Execute tblastn search against the fna (contigs) database. ([#106](https://github.com/metagenlab/zDB/pull/106)
- Handle groups when selecting genomes wherever pertinent. ([#84](https://github.com/metagenlab/zDB/pull/84) and [#85](https://github.com/metagenlab/zDB/pull/85)) (Niklaus Johner)
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

- Add button to blast protein sequence in the locus overview. ([#104](https://github.com/metagenlab/zDB/pull/104)) (Niklaus Johner)
- Integrate link to paper blast in the locus overview. ([#103](https://github.com/metagenlab/zDB/pull/103)) (Niklaus Johner)
- Add button to draw custom phylogeny from selection of entries in various tables. ([#92](https://github.com/metagenlab/zDB/pull/92)) (Niklaus Johner)
- Add buttons to download protein and DNA sequences in homologs tables. ([#89](https://github.com/metagenlab/zDB/pull/89)) (Niklaus Johner)
- Control data exports from tables by making rows selectable. ([#88](https://github.com/metagenlab/zDB/pull/88)) (Niklaus Johner)
- Add views to add, delete, and display groups. ([#86](https://github.com/metagenlab/zDB/pull/86)) (Niklaus Johner)
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

- Fix filtering of COGs in annotation pipeline. ([#108](https://github.com/metagenlab/zDB/pull/108)
- Fix blast view for hits with identifier containing a version number. ([#106](https://github.com/metagenlab/zDB/pull/106)
- Fix PFAM identifiers missing leading 0. ([#105](https://github.com/metagenlab/zDB/pull/105)) (Niklaus Johner)
- Fix product representation in GWAS view for orthogroups. ([#102](https://github.com/metagenlab/zDB/pull/102)) (Niklaus Johner)
- Fix Scoary-2 integration. ([#101](https://github.com/metagenlab/zDB/pull/101)) (Niklaus Johner)
- Fix search_bar results view when all loci miss any gene annotation. ([#99](https://github.com/metagenlab/zDB/pull/99)) (Niklaus Johner)
- Fix locusx view error for loci with homologs without any gene annotation. ([#99](https://github.com/metagenlab/zDB/pull/99)) (Niklaus Johner)
- Fix extract_contigs view error for contigs without any gene annotation. ([#99](https://github.com/metagenlab/zDB/pull/99)) (Niklaus Johner)
- Correctly handle group_0 in EntryIdParser. ([#93](https://github.com/metagenlab/zDB/pull/93)) (Niklaus Johner)
- Fix webapp conda environment for macOS. ([#87](https://github.com/metagenlab/zDB/pull/87)) (Niklaus Johner)
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

