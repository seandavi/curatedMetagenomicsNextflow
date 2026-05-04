set shell := ["bash", "-euo", "pipefail", "-c"]

[group('test')]
[default]
test:
    @just --list --unsorted

[group('test')]
test-config-local:
    nextflow config -profile local

[group('test')]
test-config-all:
    for profile in local anvil alpine google unitn; do \
      echo "==> $profile"; \
      nextflow config -profile "$profile" >/dev/null; \
    done

[group('test')]
test-stub:
    env NXF_DISABLE_CHECK_LATEST=true nextflow run . \
      -profile local \
      -stub-run \
      -c conf/test/disable-telemetry.config \
      --sample_id TEST_SAMPLE \
      --run_ids SRR000001

[group('test')]
test-nf-list:
    nf-test list

[group('test')]
test-nf:
    nf-test test

[group('test')]
test-nf-verbose:
    nf-test test --verbose

[group('test')]
test-clean:
    rm -rf .nf-test test-results null
