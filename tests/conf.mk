SHELL=/bin/bash -e -u -o pipefail

.SECONDARY:

PROP_DIR=../../prophyle

ASM=$(PROP_DIR)/prophyle_assembler/prophyle_assembler
IND=$(PROP_DIR)/prophyle_index/prophyle_index
BWA=$(PROP_DIR)/prophyle_index/bwa/bwa

TEST_NEWICK=$(PROP_DIR)/prophyle_validate_tree.py

F2K=$(PROP_DIR)/_fa_to_kmers.py
AK=$(PROP_DIR)/_all_kmers.py
NORM=$(PROP_DIR)/_fa_norm.py
1STEP=$(PROP_DIR)/1step_match.py

FQ=../simulation_bacteria.1000.fq
FA=index.fa

d=$(shell pwd)
$(info )
$(info Directory: $(d))
$(info )
