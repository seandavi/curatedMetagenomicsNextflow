# HUMAnN 4.0.0a1 for the humann4.0.0a1 bundle (conf/humann_bundles.config, #87).
#
# Not on bioconda or biobakery's Docker Hub, and biobakery's conda package
# resolves MetaPhlAn 4.2.x, which HUMAnN 4.0.0a1 cannot drive. MetaPhlAn 4.1.1
# is installed alongside because `humann` refuses to start without a metaphlan
# executable, even when given --taxonomic-profile; the bundle's MetaPhlAn pass
# itself runs from the metaphlan:4.1.1 biocontainer.
#
# Python 3.11, not 3.12: MinPath breaks under 3.12
# (https://forum.biobakery.org/t/troubleshooting-humann-in-2026-ce/9074).
# diamond/bowtie2/glpk come from conda; HUMAnN's own dependency download is
# bypassed. Published to ghcr.io/seandavi/humann:4.0.0a1 by
# .github/workflows/humann4a-image.yml.
FROM mambaorg/micromamba:2.3.3

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        pip \
        diamond=2.0.15 \
        bowtie2=2.5.4 \
        glpk=5.0 \
        procps-ng \
    && micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir --no-binary humann humann==4.0.0a1 metaphlan==4.1.1 \
    && humann --version \
    && metaphlan --version

ENV PATH=/opt/conda/bin:$PATH
