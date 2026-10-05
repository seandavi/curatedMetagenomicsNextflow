# HUMAnN 4.0.0a1 image for the `humann4.0.0a1` bundle (conf/humann_bundles.config,
# ADR-0016, issue #87).
#
# HUMAnN 4.0.0a1 is on the biobakery conda channel only (not bioconda), and
# biobakery's Docker Hub images stop at 3.9, so this image is built here.
#
# The biobakery `humann=4.0.0a1` package depends on an unpinned `metaphlan`,
# which resolves to MetaPhlAn >= 4.2; those releases reject `--bowtie2out`,
# which HUMAnN 4a passes (humann/search/prescreen.py). MetaPhlAn 4.1.2 is no
# good either: its `--version` prints a second line ("No complete MetaPhlAn
# Bowtie2 database found" / "Installed databases: ..."), and HUMAnN reads the
# version from the last line (humann/config.py metaphlan_version "line": -1),
# so HUMAnN exits with "Can not call software version for metaphlan" even when
# given --taxonomic-profile. MetaPhlAn is therefore pinned to 4.1.1, the same
# release as the bundle's MetaPhlAn pass (metaphlan:4.1.1 biocontainer), whose
# profile reaches HUMAnN via --taxonomic-profile; this copy is never used to
# profile.
#
# The humann package also ships its own `bin/diamond` (2.0.15, the DIAMOND
# HUMAnN 4a is released with). It overwrites the binary of the bioconda
# `diamond` dependency, which is linked first because humann depends on it; no
# bioconda diamond 2.0.15 build is installable beside Python 3.12, so the
# dependency is pinned to 2.2.8 and the build asserts that the binary on PATH
# is HUMAnN's 2.0.15.
#
# Every package is pinned so a rebuild resolves the same environment.
#
# Build:  docker build --platform linux/amd64 -t seandavi/curatedmetagenomics:humann4.0.0a1 docker/humann4a
# Check:  docker run --rm <image> humann --version     # humann v4.0.0.alpha.1
#         docker run --rm <image> humann_test

FROM mambaorg/micromamba:2.3.3

RUN micromamba install -y -n base \
        -c conda-forge -c bioconda -c biobakery \
        python=3.12 \
        humann=4.0.0a1 \
        metaphlan=4.1.1 \
        diamond=2.2.8 \
        bowtie2=2.5.5 \
        glpk=5.0 \
        procps-ng=4.0.6 \
    && micromamba clean --all --yes

# Put the env on PATH for non-login shells: Nextflow runs `bash .command.run`
# with the image's entrypoint (which activates the env) bypassed, and so does
# a plain Dockerfile RUN.
ENV PATH=/opt/conda/bin:$PATH

# Fail the build if a pin stops holding (e.g. a different link order leaves
# bioconda's DIAMOND on PATH, or MetaPhlAn's `--version` gains the extra line
# HUMAnN cannot parse).
RUN humann --version 2>&1 | grep -qx 'humann v4.0.0.alpha.1' \
    && [ "$(metaphlan --version 2>&1)" = 'MetaPhlAn version 4.1.1 (11 Mar 2024)' ] \
    && diamond version 2>&1 | grep -qx 'diamond version 2.0.15' \
    && bowtie2 --version 2>&1 | grep -q 'version 2.5.5'
