#! /usr/bin/env bash
set -f
set -e
set -o pipefail

set -v

./prophyle/prophyle.py analyze ./_index_test _test_analyze _test_bam1.bam

