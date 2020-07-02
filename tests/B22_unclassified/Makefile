.PHONY: all clean

include ../conf.mk

K=3
tree=tree.nw
index=index

all: output_wo_header.sam output.kraken
	for x in $^; do test -s $$x; done;

output_wo_header.sam: index.complete
	$(PROP) classify $(index) -f sam to_be_unclassified.fa | samtools view > $@

output.kraken: index.complete
	$(PROP) classify $(index) -f kraken to_be_unclassified.fa > $@

index.complete:
	$(PROP) index -k $(K) -A $(tree) $(index)
	touch $@

clean:
	rm -rf $(index)
	rm -f *.complete *.sam *.kraken

