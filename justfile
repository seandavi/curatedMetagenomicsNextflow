set shell := ["bash", "-euo", "pipefail", "-c"]

nf_test := "./tools/bin/nf-test"

[group('test')]
[default]
test:
    @just --list --unsorted

[group('test')]
test-bootstrap:
    mkdir -p tools/bin
    if [[ ! -x {{nf_test}} ]]; then \
      echo "Installing nf-test into tools/bin"; \
      curl -fsSL https://get.nf-test.com | bash; \
      mv nf-test {{nf_test}}; \
      chmod +x {{nf_test}}; \
    else \
      echo "nf-test already installed at {{nf_test}}"; \
    fi

[group('test')]
test-nf-install: test-bootstrap

[group('test')]
test-nf-version: test-bootstrap
    {{nf_test}} version

[group('test')]
test-config-local:
    nextflow config -profile local

[group('test')]
test-config-all:
    for profile in local anvil alpine google unitn; do \
      echo "==> $profile"; \
      nextflow config -profile "$profile" >/dev/null; \
    done

# The config parses under the current (v2) parser, but the process scripts
# still target the v1 script grammar, so pin the parser for anything that
# compiles main.nf/modules until that migration happens.
[group('test')]
test-stub:
    env NXF_DISABLE_CHECK_LATEST=true NXF_SYNTAX_PARSER=v1 nextflow run . \
      -profile local \
      -stub-run \
      -c conf/test/disable-telemetry.config \
      --sample_id TEST_SAMPLE \
      --run_ids SRR000001

[group('test')]
test-nf-list: test-bootstrap
    {{nf_test}} list

[group('test')]
test-nf: test-bootstrap
    NXF_SYNTAX_PARSER=v1 {{nf_test}} test

[group('test')]
test-nf-verbose: test-bootstrap
    NXF_SYNTAX_PARSER=v1 {{nf_test}} test --verbose

[group('test')]
test-clean:
    rm -rf .nf-test test-results null
