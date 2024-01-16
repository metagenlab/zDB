# Changelog


Versions are of the form MAJOR.MINOR.PATCH and this changelog tries to conform
to [Common Changelog](https://common-changelog.org)


## unreleased

### Changed

- Configure repository to publish the docs to [readthedocs](https://zdb.readthedocs.io)
- Port nextflow pipelines to DSL2 language. ([#26](https://github.com/metagenlab/zDB/pull/26)) (Niklaus Johner)
- Add CHANGELOG file. ([#25](https://github.com/metagenlab/zDB/pull/25)) (Niklaus Johner)
- Imported metagenlab_lib code into zdb (https://github.com/metagenlab/zDB/pull/33) (Bastian Marquis)

### Added

- Add AMR hit extraction view. ([#38](https://github.com/metagenlab/zDB/pull/38)) (Niklaus Johner)
- Add AMR entry list view. ([#37](https://github.com/metagenlab/zDB/pull/37)) (Niklaus Johner)
- Add AMR comparison index view and link it in the menu and home page. ([#35](https://github.com/metagenlab/zDB/pull/35)) (Niklaus Johner)
- Add AMR tabular comparison views. ([#24](https://github.com/metagenlab/zDB/pull/24), [#31](https://github.com/metagenlab/zDB/pull/31)) (Niklaus Johner)
- Add integration tests for pipelines and views. ([#29](https://github.com/metagenlab/zDB/pull/29)) (Niklaus Johner)
- Add unitesting for pipelines. ([#26](https://github.com/metagenlab/zDB/pull/26)) (Niklaus Johner)
- Add antimicrobial resistance annotations. ([#23](https://github.com/metagenlab/zDB/pull/23)) (Niklaus Johner)
- Add support for circular contigs in the genomic viewer. (https://github.com/metagenlab/zDB/pull/33) (Bastian Marquis)

### Fixed

- Fix tabular comparison views navigation tabs. ([#30](https://github.com/metagenlab/zDB/pull/30), [#31](https://github.com/metagenlab/zDB/pull/31)) (Niklaus Johner)


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

