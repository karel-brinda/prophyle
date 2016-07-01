#! /usr/bin/env python3

import sys, os, argparse, re
from collections import deque
from ete3 import PhyloTree, NCBITaxa

def read_buf(file, buf_size=536870912):
	buf = file.read(buf_size)+file.readline()
	while buf:
		yield buf
		buf = file.read(buf_size)+file.readline()

def index_of(taxid, taxa_list):
	for i, seq in enumerate(taxa_list):
		if taxid == seq[3]:
			return i
	return -1

parser = argparse.ArgumentParser(
	description='Assigns taxids to the sequences in the input fasta file '
				'and builds a taxonomic tree in the newick format')

parser.add_argument('library_dir', help = 'directory containing sequences\' fasta files')
parser.add_argument('taxid_map_f', help = 'map of gis to taxid from NCBI database')
parser.add_argument('-o', default = 'taxonomic_tree.nw',
					help = 'output file (default: taxonomic_tree.nw)')
parser.add_argument('-e', default = 'error_tree.log',
					help = 'error log file (default: error_tree.log)')
parser.add_argument('--assign_only', action = 'store_true',
					help = 'only assign the taxids to the sequences and write them in taxids.dmp')
parser.add_argument('--build_tree', default = None,
					help = 'build the tree from the given file (previous output of --assign_only)')

args = parser.parse_args()
ncbi = NCBITaxa()
error = open(args.e, 'w')
gis = []
taxids = [-1]

if args.build_tree is not None:
	with open(args.build_tree, 'r') as seq_taxid:
	    for line in seq_taxid:
	        values = line.split("\t")
	        gis.append([int(values[0]),str(values[1]),str(values[2]),int(values[3])])
else:
	skipped = 0
	for f in os.listdir(args.library_dir):
		if f.endswith(".fna"):
			with open(args.library_dir+'/'+f, 'r') as fasta:
				seq = fasta.readline()
				(seqname,_,_) = seq.partition(" ")
				values = seqname.split("|")
				gi = values[2]
				try:
					gis.append([int(gi),str(seqname)[1:],f])
				except:
					if skipped == 0:
						error.write("NOT ACQUIRED:\n\n")
					error.write(f + str(seqname) + "\n")
					skipped += 1
					pass

	gis = sorted(gis, key = lambda x:x[0])
	gis_no = len(gis)
	print("Acquired " + str(gis_no) + " gis (" + str(skipped) + " skipped)")

	skipped = 0
	with open(args.taxid_map_f, 'r') as taxid_map:
		i = 0
		found = False
		for buf in read_buf(taxid_map):
			for line in buf.splitlines():
				(gi, _, taxid) = line.partition("\t")
				gi = int(gi)
				while i < len(gis) and gis[i][0] < gi:
					if not found:
						if skipped == 0:
							error.write("\n\nNOT ASSIGNED:\n\n")
						error.write(str(gis[i][2]) + str(gis[i][1])+"\n")
						del(gis[i])
						skipped += 1
					else:
						i += 1
						found = False
				if i < len(gis) and gis[i][0] == gi:
					found = True
					gis[i].append(taxid.strip())

	print("Assigned " + str(gis_no-skipped) + " gis (" + str(skipped) + " skipped)")

	if args.assign_only:
		if args.o == 'taxonomic_tree.nw':
			args.o = 'taxids_assignments.dmp'
		with open(args.o, 'w') as output:
			for i in gis:
				output.write("\t".join(map(str,i))+"\n")
		print("Assignments written to " + args.o)
		print("Launch again with --build_tree " + args.o +
				" to build a taxonomic tree from them")
		sys.exit(0)

gis = sorted(gis, key = lambda x:x[3])
prec = 0
for s in gis:
	tid = s[3]
	if tid != taxids[prec]:
		taxids.append(tid)
		prec += 1
del taxids[0]

topo = ncbi.get_topology(taxids)
new_id = 1
count = 0
lib_dir = args.library_dir[args.library_dir.find("library"):]+'/'
for node in topo.traverse("preorder"):
	node.name = new_id
	new_id += 1
	i = index_of(node.taxid, gis)
	if i != -1:
		gi = ""
		seqname = ""
		fastapath = ""
		infasta_seqnum = ""
		while i < len(gis) and gis[i][3] == node.taxid:
			if gi == "":
				gi = str(gis[i][0])
				seqname = gis[i][1]
				fastapath = lib_dir+gis[i][2]
				infasta_seqnum = '1'
			else:
				gi += "@"+str(gis[i][0])
				seqname += "@"+gis[i][1]
				fastapath += "@"+lib_dir+gis[i][2]
				infasta_seqnum += "@"+'1'
			i += 1
			count += 1
		node.add_features(gi = gi, seqname = seqname,
			fastapath = fastapath, infasta_seqnum = infasta_seqnum)

print("Built taxonomic tree for " + str(count) + " sequences")

topo.write(features=["lineage", "named_lineage", "seqname", "dist",
					"name", "support", "taxid", "rank", "fastapath",
					"sci_name", "common_name", "infasta_seqnum", "gi"],
			outfile=args.o)

error.close()
