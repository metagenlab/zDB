# Changelog


Versions are of the form MAJOR.MINOR.PATCH and this changelog tries to conform
to [Common Changelog](https://common-changelog.org)


## 1.2.1 - 2023-10-16

### Changed
- Bump checkm-genome version to 1.2.2. (Trestan Pillonel)


## 1.2.0 - 2023-10-13

### Changed

- Use mamba to build conda envs ([30bbb7530d75c3adf986aaf4c5052ad83ea301b8](https://github.com/metagenlab/zDB/commit/30bbb7530d75c3adf986aaf4c5052ad83ea301b8)) (Trestan Pillonel)

### Fixed

- Fix issue for genome file download ([#16](https://github.com/metagenlab/zDB/pull/16)) (Trestan Pillonel)
- Fix file extensions in annotation pipeline ([#17](https://github.com/metagenlab/zDB/pull/17)) (Thomas Kozy)


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
