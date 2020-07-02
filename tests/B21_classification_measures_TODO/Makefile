.PHONY: all clean

include ../conf.mk

K=3
tree=tree.nw
index=index

all: index.complete

index.complete:
	$(PROP) index -k $(K) -A $(tree) $(index)
	touch $@

clean:
	rm -rf $(index) $(res) $(diff_res) *.complete
