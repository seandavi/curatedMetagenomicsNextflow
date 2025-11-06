# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2025-10-30

### Added

- Complete nf-core refactoring of pipeline structure
- Modularized all processes following nf-core DSL2 conventions
- Added `workflows/` directory with main workflow
- Added `modules/local/` directory with individual process modules
- Added `conf/` directory with configuration files:
  - `base.config` for base process configuration
  - `modules.config` for module-specific settings
  - `test.config` for test profile
- Added `.nf-core.yml` configuration file
- Added GitHub Actions CI/CD workflows for linting and testing
- Added module metadata files (`meta.yml`) for documentation
- Updated README with comprehensive nf-core-style documentation
- Standardized parameter naming (`--input`, `--outdir`)
- Added `check_max()` function for resource management
- Added support for multiple container engines (Docker, Singularity, Podman)

### Changed

- Refactored monolithic `main.nf` into modular structure
- Updated `nextflow.config` to follow nf-core conventions
- Improved parameter handling with backwards compatibility
- Enhanced error handling and validation
- Updated manifest information
- Improved output organization with `pipeline_info` directory

### Improved

- Better separation of concerns with modular architecture
- Easier maintenance and updates
- More portable across different compute environments
- Better documentation and help messages
- Improved resource allocation with labels
- Enhanced container management

## [1.2.0] and earlier

Previous versions before nf-core refactoring. See git history for details.
