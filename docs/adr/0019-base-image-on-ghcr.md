# 0019. Build the base image in GitHub Actions, publish to GHCR, drop HUMAnN from it

- **Status:** Accepted
- **Date:** 2026-10-08
- **Deciders:** Sean Davis

## Context

2.3.0 is the first full corpus run, and it moves the main taxonomy pass to
MetaPhlAn 4.2.6 with the vJan26 index (#105). The main pass runs in the base
image (ADR-0018: the default profile's container equals `process.container`),
so the base image has to be rebuilt.

The base image (`seandavi/curatedmetagenomics:metaphlan4.2.2`) was built by a
manually triggered Google Cloud Build and pushed to Docker Hub with personal
credentials. Nothing in the repo records what was built from what, and agents
have no registry credentials. ADR-0017 set up the Actions → GHCR route for
self-built tool images and explicitly left the base image out of it.

The base image also carried HUMAnN 4.0.0a1, which no process has used since
HUMAnN moved to bundles with their own containers (ADR-0016, #89).

## Decision

- The base image is built from `docker/Dockerfile` by
  `.github/workflows/base-image.yml`, as ADR-0017 does for tool images. Pull
  requests build and smoke-test it, and merges to `main` publish
  `ghcr.io/seandavi/curatedmetagenomics:metaphlan<version>` (plus a
  `-<commit sha>` tag). The tag is named after the MetaPhlAn version, as before.
- HUMAnN is removed from the base image.
- The default profile `mpa4.2.6_vJan26` and `process.container` both name
  `ghcr.io/seandavi/curatedmetagenomics:metaphlan4.2.6`, and CI checks that they
  stay equal.
- The Cloud Build files (`docker/cloudbuild.yaml`, `docker/docker_build.sh`) are
  deleted.

The rest of the toolchain stays at the 4.2.2 image's versions, on Debian 12
(`python:3.9-bookworm`), so the new image differs from the old one only in
MetaPhlAn and HUMAnN. MetaPhlAn 4.2.6 is not on PyPI, so it is installed from
the hash-pinned GitHub release tarball, the same source bioconda builds from.

## Alternatives considered

- **Keep the manual Cloud Build / Docker Hub route.** It needs personal
  credentials and a manual trigger, and leaves no record of the inputs.
- **Run the main pass in the MetaPhlAn 4.2.6 biocontainer.** ADR-0018 deferred
  this. KneadData and the other base-image tools would still need a base image,
  and it would split the main pass across two containers for no gain here.
- **Keep HUMAnN in the base image.** It's dead weight and a second, unused
  HUMAnN version.

## Consequences

- The base image is reproducible from the repo, and every change to the
  Dockerfile is built and smoke-tested on its PR.
- The new GHCR package starts private. A maintainer has to make it public once
  (as #97 did for `humann`) before clusters can pull it without credentials.
- The 2.2.x image stays on Docker Hub, and `mpa4.2.2_vJan25` still names it.
- MetaPhlAn 4.2.6 reports itself as 4.2.5, so the manifest records 4.2.5 for
  the default profile.

## References

- #105 (2.3.0), #89 (HUMAnN out of the base image), #97 (GHCR visibility).
- ADR-0001 (container strategy), ADR-0016 (HUMAnN bundles), ADR-0017 (GHCR for
  self-built images), ADR-0018 (MetaPhlAn profiles).
