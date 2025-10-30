# nf-core Refactoring Summary

This document summarizes the complete nf-core refactoring of the curatedMetagenomicsNextflow pipeline.

## Overview

The pipeline has been successfully refactored from a monolithic structure to a modular nf-core-compatible architecture following DSL2 best practices.

## Changes Made

### 1. Directory Structure

**Before:**
```
curatedMetagenomicsNextflow/
├── main.nf (800+ lines)
├── nextflow.config
├── nextflow_schema.json
└── docker/
```

**After:**
```
curatedMetagenomicsNextflow/
├── main.nf (73 lines)
├── nextflow.config (enhanced)
├── nextflow_schema.json
├── .nf-core.yml (NEW)
├── workflows/
│   └── curatedmetagenomicsnextflow.nf (NEW)
├── modules/
│   └── local/
│       ├── fasterq_dump/
│       ├── local_fastqc/
│       ├── kneaddata/
│       ├── install_metaphlan_db/
│       ├── metaphlan_unknown_viruses_lists/
│       ├── metaphlan_unknown_list/
│       ├── metaphlan_markers/
│       ├── sample_to_markers/
│       ├── chocophlan_db/
│       ├── uniref_db/
│       ├── kneaddata_human_database/
│       ├── kneaddata_mouse_database/
│       └── humann/
├── subworkflows/
│   └── local/ (ready for future use)
├── conf/
│   ├── base.config (NEW)
│   ├── modules.config (NEW)
│   └── test.config (NEW)
├── assets/ (NEW)
├── bin/ (NEW)
├── .github/
│   └── workflows/
│       ├── ci.yml (NEW)
│       └── linting.yml (NEW)
├── CHANGELOG.md (NEW)
├── CITATIONS.md (NEW)
├── CODE_OF_CONDUCT.md (NEW)
└── README.md (completely rewritten)
```

### 2. Module Extraction

All processes have been extracted into individual modules:

1. **FASTERQ_DUMP** - Download FASTQ from SRA
2. **LOCAL_FASTQC** - Process local FASTQ files
3. **KNEADDATA** - Quality control and decontamination
4. **INSTALL_METAPHLAN_DB** - MetaPhlAn database installation
5. **METAPHLAN_UNKNOWN_VIRUSES_LISTS** - Taxonomic profiling with virus detection
6. **METAPHLAN_UNKNOWN_LIST** - Generate taxonomic profiles
7. **METAPHLAN_MARKERS** - Extract marker information
8. **SAMPLE_TO_MARKERS** - Generate StrainPhlAn markers
9. **CHOCOPHLAN_DB** - Download ChocoPhlAn database
10. **UNIREF_DB** - Download UniRef database
11. **KNEADDATA_HUMAN_DATABASE** - Download human reference
12. **KNEADDATA_MOUSE_DATABASE** - Download mouse reference
13. **HUMANN** - Functional profiling

Each module includes:
- `main.nf` - Process definition
- `meta.yml` - Module metadata and documentation (for key modules)

### 3. Configuration Improvements

#### conf/base.config
- Defines base process settings
- Implements resource labels (process_single, process_low, process_medium, process_high)
- Includes check_max() function for resource management
- Standard error handling strategies

#### conf/modules.config
- Module-specific configurations
- PublishDir settings per module
- Extension arguments support

#### conf/test.config
- Minimal test dataset configuration
- Resource limits for CI/CD
- Skip computationally expensive steps

### 4. Parameter Standardization

**New nf-core standard parameters:**
- `--input` (replaces `--metadata_tsv`)
- `--outdir` (replaces `--publish_dir`)
- `--publish_dir_mode` (replaces `--publish_mode`)
- `--max_cpus`, `--max_memory`, `--max_time`
- `--help`, `--version`

**Backwards compatibility maintained:**
- `--metadata_tsv` still works (maps to `--input`)
- `--publish_dir` still works (maps to `--outdir`)

### 5. Workflow Structure

**New workflow file:** `workflows/curatedmetagenomicsnextflow.nf`
- Clean separation of concerns
- Reusable helper functions
- Improved readability
- Better maintainability

**Updated main.nf:**
- Minimal entry point
- Parameter validation
- Help message
- Imports main workflow

### 6. Documentation

#### README.md
- Comprehensive usage guide
- Installation instructions
- Parameter documentation
- Profile descriptions
- Output structure
- Quick start examples

#### CHANGELOG.md
- Version history
- Detailed change log
- Migration guide

#### CITATIONS.md
- All tool citations
- Proper attribution
- DOI links

#### CODE_OF_CONDUCT.md
- Community standards
- Contributor guidelines

#### Module meta.yml files
- Input/output specifications
- Tool descriptions
- Keywords and authors

### 7. CI/CD

#### .github/workflows/ci.yml
- Automated testing
- Multiple Nextflow versions
- Stub run testing
- Profile testing

#### .github/workflows/linting.yml
- Code quality checks
- nf-core lint
- Pre-commit hooks
- Prettier formatting
- EditorConfig validation

### 8. Container Management

Enhanced container support:
- Docker profile
- Singularity profile
- Podman support
- Standardized container declarations

## Benefits

1. **Modularity**: Each process is independent and reusable
2. **Maintainability**: Easier to update individual components
3. **Portability**: Better support for different execution environments
4. **Testability**: Individual modules can be tested separately
5. **Documentation**: Better inline and external documentation
6. **Standards**: Follows nf-core best practices
7. **Compatibility**: Works with nf-core tooling
8. **Scalability**: Easier to add new processes or features

## Testing

The refactored pipeline can be tested with:

```bash
# Stub run (quick validation)
nextflow run . -profile test,docker -stub-run

# Full test run
nextflow run . -profile test,docker --outdir results
```

## Migration Guide

For users of the previous version:

### Old Command:
```bash
nextflow run main.nf --metadata_tsv samples.tsv --publish_dir results
```

### New Command (recommended):
```bash
nextflow run seandavi/curatedmetagenomicsnextflow --input samples.tsv --outdir results -profile docker
```

### Backwards Compatible:
```bash
nextflow run seandavi/curatedmetagenomicsnextflow --metadata_tsv samples.tsv --publish_dir results
```

## Files Modified

- `main.nf` - Complete rewrite (800+ → 73 lines)
- `nextflow.config` - Enhanced with nf-core conventions
- `README.md` - Complete rewrite with comprehensive documentation

## Files Added

- `.nf-core.yml`
- `workflows/curatedmetagenomicsnextflow.nf`
- `modules/local/*/main.nf` (13 modules)
- `modules/local/*/meta.yml` (5 metadata files)
- `conf/base.config`
- `conf/modules.config`
- `conf/test.config`
- `.github/workflows/ci.yml`
- `.github/workflows/linting.yml`
- `CHANGELOG.md`
- `CITATIONS.md`
- `CODE_OF_CONDUCT.md`

## Validation

The refactored pipeline:
- ✅ Maintains all original functionality
- ✅ Follows nf-core directory structure
- ✅ Uses nf-core parameter conventions
- ✅ Implements proper resource management
- ✅ Includes comprehensive documentation
- ✅ Has CI/CD pipelines
- ✅ Supports multiple container engines
- ✅ Is backwards compatible

## Next Steps

Future enhancements could include:
1. Creating subworkflows for related processes (e.g., METAPHLAN_PROFILE)
2. Adding more module metadata files
3. Implementing input validation schema
4. Adding more comprehensive tests
5. Submission to nf-core (if desired)

## Conclusion

This refactoring successfully transforms the pipeline into a modern, maintainable, and standards-compliant bioinformatics workflow that follows nf-core best practices while maintaining full backwards compatibility with the original implementation.
