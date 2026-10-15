#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python tests/make_fixture.py
python tests/unit.py
snakemake --configfile tests/config.yaml --cores 2 --dry-run
snakemake --configfile tests/config.yaml --cores 2
python tests/check_results.py
snakemake --configfile tests/config.yaml --cores 2 --dry-run
snakemake --configfile tests/config.empty.yaml --cores 2
python tests/check_empty.py
