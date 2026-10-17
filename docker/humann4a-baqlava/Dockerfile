# HUMAnN 4.0.0a1 + BAQLaVa 1.2.0 for the optional baqlava step (#117, #119).
#
# Built on the humann4.0.0a1 bundle image (docker/humann4a, ADR-0017) so the
# HUMAnN that BAQLaVa drives for bacterial depletion and its searches is the
# same one the bundle runs. BAQLaVa needs HUMAnN >= 4 (it checks only the major
# version), AnADAMA2 >= 0.10.0 and Python >= 3.10 (the base has 3.11).
#
# BAQLaVa 1.2.0 is not on PyPI (it stops at 0.5) and has no `v1.2.0` tag: there
# is a `v1.2.0` branch and a `v1.2.0-devel` release tag (629c744). The pin is the
# branch tip, 9a458b1, one commit past the release, which fixes task scripts
# invoking a bare `python`. Its setup.py takes the package version from
# configs/baqlava.cfg ('1.2.0-devel'), which current setuptools rejects as not
# PEP 440, so only setup.py's copy is rewritten to 1.2.0.dev0; `baqlava
# --version` still reads the untouched config and reports 1.2.0-devel.
#
# AnADAMA2 comes from bioconda, which ships the compiled leveldb binding (the
# PyPI sdist needs a C++ toolchain). BAQLaVa's workflow copies local files
# with `scp`; real scp exits 255 ("unknown user") when the container uid has no
# passwd entry, so `scp` is a symlink to `cp`. The base's conda pins are
# restated so the solver cannot move them (unpinned, it upgrades diamond to
# 2.1). Databases are not baked in; the pipeline downloads them with
# `baqlava_database` into store_dir.
#
# Published to ghcr.io/seandavi/humann-baqlava:4.0.0a1-1.2.0 by
# .github/workflows/humann4a-baqlava-image.yml.
FROM ghcr.io/seandavi/humann:4.0.0a1

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        anadama2=0.10.0 \
        python=3.11 \
        diamond=2.0.15 \
        bowtie2=2.5.4 \
        glpk=5.0 \
    && micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG BAQLAVA_COMMIT=9a458b1c3f88153be56f167cb85ebe828f15f6ba
RUN cd /tmp \
    && python -c "import urllib.request as u; u.urlretrieve('https://github.com/biobakery/baqlava/archive/${BAQLAVA_COMMIT}.tar.gz', 'baqlava.tar.gz')" \
    && tar xzf baqlava.tar.gz \
    && sed -i 's/^VERSION = .*/VERSION = "1.2.0.dev0"/' baqlava-${BAQLAVA_COMMIT}/setup.py \
    && pip install --no-cache-dir --no-deps ./baqlava-${BAQLAVA_COMMIT} \
    && rm -rf /tmp/baqlava* \
    && humann --version \
    && baqlava --version \
    && ln -s "$(command -v cp)" /opt/conda/bin/scp
