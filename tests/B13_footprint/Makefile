.PHONY: all clean footprint

include ../conf.mk

K=3
tree=tree.nw
index=index

all: footprint

footprint: index.complete
	$(PROP) footprint $(index)

index.complete:
	$(PROP) index -k $(K) -A $(tree) $(index)
	touch $@

clean:
	rm -rf $(index) $(res) $(diff_res) *.complete *.sam *.krak
