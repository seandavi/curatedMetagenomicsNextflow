# 0017. Publish self-built tool images to GHCR from GitHub Actions

- **Status:** Accepted
- **Date:** 2026-10-05
- **Deciders:** Sean Davis

## Context

ADR-0001 runs new tools in pinned upstream biocontainers, and the only image we
build ourselves is the base image (`seandavi/curatedmetagenomics`), built by a
manually triggered Google Cloud Build and pushed to Docker Hub with personal
credentials.

The `humann4.0.0a1` bundle (#87) needs HUMAnN 4.0.0a1. It is not on bioconda;
biobakery's Docker Hub images stop at 3.9; and biobakery's own conda package
resolves MetaPhlAn 4.2.x, which HUMAnN 4.0.0a1 cannot drive (it passes
`--bowtie2out`, removed in 4.2). No upstream image is usable, so we must build
one. Agents doing this work have no registry credentials, and the base image's
manual build doesn't leave a reproducible record of what was built from what.

## Decision

Self-built tool images live in `docker/<name>/` and are built and published by a
GitHub Actions workflow per image, to `ghcr.io/seandavi/<name>:<tool version>`
(plus a `<version>-<commit sha>` tag), using the workflow's `GITHUB_TOKEN`.
Pull requests build and smoke-test the image without pushing; merges to `main`
publish. Packages are made public once by a maintainer (#97) so clusters pull
without credentials. The first image is `docker/humann4a` →
`ghcr.io/seandavi/humann:4.0.0a1`.

Upstream biocontainers stay the default (ADR-0001); this applies only when no
usable upstream image exists. The base image is unchanged.

## Alternatives considered

- **Docker Hub via the existing manual Cloud Build** — needs personal
  credentials and a manual trigger for every change; not reproducible from the
  repo alone.
- **Use biobakery's conda package at run time (Nextflow `conda`)** — resolves
  the incompatible MetaPhlAn 4.2.x unless pinned (proposed upstream in
  biobakery/conda-biobakery#3), and solves an environment on every cluster.
- **Bake HUMAnN 4a into the base image** — re-couples HUMAnN to the base image,
  which ADR-0016 just separated.

## Consequences

- Images are rebuilt from `main` on every change to their directory, and the
  build log records the inputs. Bundles pin the version tag; the PR notes the
  digest.
- New GHCR packages start private; the one-time visibility change is manual
  (#97). Until then the image can be used on a cluster only as a locally built
  SIF.
- A second image provenance model (Actions → GHCR) now exists beside the base
  image's (Cloud Build → Docker Hub). Moving the base image to the same workflow
  is possible later but not decided here.

## References

- #87 (bundle and image), #97 (public visibility), #82 (plan).
- `.github/workflows/humann4a-image.yml`, `docker/humann4a/Dockerfile`.
- ADR-0001 (container strategy), ADR-0016 (HUMAnN bundles).
