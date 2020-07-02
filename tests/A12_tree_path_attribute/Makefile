.PHONY: all clean diff index_fastapath index_path

include ../conf.mk

K=5

all: diff

diff: index_fastapath index_path
	diff -c index_fastapath/index.fa index_path/index.fa

index_fastapath:
	$(PROP) index -k $(K) tree_fastapath.nw $@

index_path:
	$(PROP) index -k $(K) tree_fastapath.nw $@

clean:
	rm -rf index_*
