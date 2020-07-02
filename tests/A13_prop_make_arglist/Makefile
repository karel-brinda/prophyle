.PHONY: all clean _index _index_noprop

include ../conf.mk

K=5

all: _index _index_noprop

_index_noprop: lib.complete
	$(PROP) index -R -T -k $(K) -g lib tree.nw $@
	find $@/propagation -name '*.fa' | xargs -I {} bash -c 'diff -c <($(F2K) -i "{}" -k $(K) -m c -f fa) expected_noprop/`basename "{}"`'

_index: lib.complete
	$(PROP) index -T -k $(K) -g lib tree.nw $@
	find $@/propagation -name '*.fa' | xargs -I {} bash -c 'diff -c <($(F2K) -i "{}" -k $(K) -m c -f fa) expected/`basename "{}"`'

lib.complete:
	cd lib && for ((i=6; i<=10000; i++)); do ln -snf 5.fa "$$i.fa" ; done
	touch $@

clean:
	rm -rf _index* *.complete
	for ((i=6; i<=10000; i++)); do rm -f "lib/$$i.fa" ; done
