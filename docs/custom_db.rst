.. _custom_db:

Building a custom database
==========================

.. contents::
	:depth: 3
	:local:
	:backlinks: none


Building a custom database with user-provided data
--------------------------------------------------





Building a custom database with RefSeq genomes and the NCBI taxonomy
--------------------------------------------------------------------


Download sequences from NCBI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Create a subdirectory in ProPhyle's main directory, and change the working directory to the newly created one
* Download an ``assembly_summary.txt`` file from NCBI's FTP server. For RefSeq's bacterial genomes, follow this `link <ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt>`_. More information about these files can be found `here <ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt>`_
* Select genomes of interest; the following command will select all complete genomes, but it is sufficient to add more conditions to ``awk`` to select a subset. Fields 1, 6 and 20 contain respectively the accession number, taxonomic identifier and ftp directory of each genome. ``acc2taxid.tsv`` will therefore be a tab separated file containing the taxid corresponding to each sequence's accession number. Please refer to `NCBI's genomes download FAQ <https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq>`_ for further information). Since strain taxids are deprecated, we are working on a solution to create artificial nodes for each genome. For the moment, you can associate the sequences to the species node, by selecting field 7 instead of 6 in the assembly summary, which corresponds to the species taxid.

.. code-block:: bash

	awk -F "\t" -v OFS="\t" '$12=="Complete Genome" && $11=="latest"\
		{split($1, a, "."); print a[1], $7, $20}' assembly_summary.txt >ftpselection.tsv
	cut -f 3 ftpselection.tsv | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}\
		{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' >ftpfilepaths.tsv
	cut -f 1,2 ftpselection.tsv >acc2taxid.tsv

* Download selected sequences using ``parallel --gnu -j 32 -a ftpfilepaths.tsv wget``. Number of jobs (``-j`` option) should be changed according to the number of available cores and bandwidth.
* Extract them using ``parallel --gnu -j 32 -a ftpfilepaths.tsv gunzip {/}``
* Create a `fasta index <http://www.htslib.org/doc/faidx.html>`_ for each file using ``find . -name '*.fna' | parallel -j 32 samtools faidx {}``


Build a tree
~~~~~~~~~~~~

Build a taxonomic tree for the downloaded sequences using

.. code-block:: bash

	prophyle_ncbi_tree.py <library_name> <library_main_dir>\
		<output_file> acc2taxid.tsv -l <log_file>

``library_name`` is the name of the subdirectory where the sequences are downloaded (e.g. bacteria, archaea, viral for the standard DBs, or any name of your choice), which will be prepended to the filenames in the tree. ``library_main_dir`` is the directory where all the sequences are downloaded, which usually is the parent of library_name. This should be used at indexing time (the -g parameter of prophyle index), unless the tree is placed in the very same directory.
Taxonomic identifiers are assigned to the sequences first, and then the tree is
built using `ETE Toolkit <http://etetoolkit.org/>`_ and saved with newick format
1. Necessary node attributes are:

* ``name``: unique node name (for RefSeq DB: the taxid of the node)
* ``path``: paths of the sequences' fasta files, separated by @ (relative paths from ProPhyle's home directory)

Other optional attributes are ``taxid``, ``sci_name``, ``named_lineage``, ``lineage``, ``rank`` (more info
`here <http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html#automatic-tree-annotation-using-ncbi-taxonomy>`_
). Using ETE library or modifying the ``prophyle_ncbi_tree.py`` script it is
easy to adapt any tree to match the requirements above.


Build the index
~~~~~~~~~~~~~~~

Run the standard command to build ProPhyle index

.. code-block:: bash

	prophyle index -k <kmer_length> [-g <library_main_dir>] <tree_1.nw> [<tree_2.nw> ...] <index_dir>
