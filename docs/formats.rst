.. _formats:

File formats
============

.. contents::
	:depth: 3
	:local:
	:backlinks: none






Input data formats
------------------


Node names should not contain the ``@`` character, and they should not be equal to ``A``, ``0``, and ``root`` (reserved names).


.. index:: Newick

Newick trees
^^^^^^^^^^^^

Introduction
""""""""""""

The Newick format can be used for index construction
in combination with the ``-A`` parameter.
Names of files with sequences will be inferred from the names of leaves
as ``[node_name].fa``.
If names of internal nodes are not specified in the original tree, they will be assigned automatically
as the lexigraphically minimal name of children's names with incremented ID.
Branch lenghts are ignored.

Specification
"""""""""""""

See specifications of Newick on the
`Phylip website <http://evolution.genetics.washington.edu/phylip/newicktree.html>`_
or on
`Wikipedia <https://en.wikipedia.org/wiki/Newick_format>`_.

Examples
""""""""

A Newick tree with named leaves::

	((n1,n2,n3),(n5,n6));


A Newick tree with named nodes::

	((n1,n2)o1,(n3,n4,n5)o2)p1;

A Newick tree with automatically assigned names of internal node names::

	((n1,n2,n3)n1-up1,(n4,n5)n4-up1)n1-up2;


.. index:: NHX

NHX trees
^^^^^^^^^

Introduction
""""""""""""

`New Hampshire X Format <https://sites.google.com/site/cmzmasek/home/software/forester/nhx>`_
is parsed using the `ETE3 library <http://etetoolkit.org/>`_  (see specification of `Format 1 <http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees>`_).

Specification
"""""""""""""

	.. list-table:: NHX attributes
	   :widths: 7 7 20
	   :header-rows: 1

	   * - Attribute
	     - Type
	     - Description
	   * - (name)
	     - string
	     - Node name (typically the TaxID of the node). The names should be unique and must not contain ``@``.
	   * - path
	     - string
	     - Files with genomic sequences, separated by ``@`` (relative paths from the directory of the tree). Only for leaves.
	   * - fastapath
	     - string
	     - Deprecated (use `path` instead).
	   * - rank
	     - string/int
	     - Taxonomic rank.
	   * - dist
	     - float
	     - To be ignored (an internal parameter of ETE3).
	   * - support
	     - float
	     - To be ignored (an internal parameter of ETE3).
	   * - kmers_full
	     - integer
	     - Number of k-mers associated with this node. Added automatically during index construction.
	   * - kmers_reduced
	     - integer
	     - Number of k-mers represented by this node. Added automatically during index construction.

Example
"""""""
Previous tree after autocompleting to NHX::

	(((n1:1[&&NHX:dist=1.0:path=n1.fa:support=1.0],n2:1[&&NHX:dist=1.0:path=n2.fa:support=1.0])o1:1[&&NHX:dist=1.0:support=1.0],(n3:1[&&NHX:dist=1.0:path=n3.fa:support=1.0],n4:1[&&NHX:dist=1.0:path=n4.fa:support=1.0],n5:1[&&NHX:dist=1.0:path=n5.fa:support=1.0])o2:1[&&NHX:dist=1.0:support=1.0])p1:0[&&NHX:dist=0.0:support=1.0])merge_root:1[&&NHX:dist=1.0:support=1.0];


Sequences
^^^^^^^^^

Input sequences can be provided in the FASTA or FASTQ formats. Any non-``ACGT`` characters are treated as
unknown nucleotides and k-mers containing them thus discarded.
Sequence names are ignored.














Assignments
-----------


.. index:: SAM, BAM

Read assignments in SAM/BAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Introduction
""""""""""""

ProPhyle uses `SAM/BAM <http://samtools.github.io/hts-specs/>`_ as
the main format for reporting the final assignments, i.e.,
the output of classification.

Specification
"""""""""""""

	.. list-table:: ProPhyle SAM headers
	   :widths: 1 3
	   :header-rows: 1

	   * - Tag
	     - Description
	   * - HD
	     - Version of SAM.
	   * - PG
	     - Version of ProPhyle.
	   * - SQ
	     - Description of a leaf. *SN:* Name of the node. *LN:* a fake length. *UR:* Name of the original FASTA file. *SP:* Name of the species (if present in the tree).

|

	.. list-table:: ProPhyle SAM fields
	   :widths: 3 3 20
	   :header-rows: 1

	   * - Column
	     - Name
	     - Description
	   * - 1
	     - QNAME
	     - Query name.
	   * - 2
	     - FLAG
	     - ``0`` if assigned, ``4`` otherwise.
	   * - 3
	     - RNAME
	     - Node name.
	   * - 4
	     - POS
	     - ``1`` if assigned, unused (``0``) otherwise.
	   * - 5
	     - MAPQ
	     - ``60`` if assigned, unused (``0``) otherwise.
	   * - 6
	     - CIGAR
	     - Coverage bit-mask encoded as a CIGAR string if assigned, unused (``*``) otherwise. For instance, `7=3X3=` means `1111111000111`.
	   * - 7
	     - RNEXT
	     - Unused (``*``).
	   * - 8
	     - PNEXT
	     - Unused (``0``).
	   * - 9
	     - TLEN
	     - Unused (``0``).
	   * - 10
	     - SEQ
	     - Sequence of bases if ``-P``, unused (``*``) otherwise.
	   * - 11
	     - QUAL
	     - Base qualities if ``-P``, unused (``*``) otherwise.

|

	.. list-table:: ProPhyle SAM tags
	   :widths: 3 3 15 7
	   :header-rows: 1

	   * - Tag
	     - Type
	     - Description
	     - Range
	   * - ln
	     - integer
	     - Read length.
	     -
	   * - h1
	     - integer
	     - Number of shared k-mers.
	     - :math:`\{1, \ldots, |query|-k+1\}`
	   * - h2
	     - float
	     - Proportion of hits in the reference.
	     - :math:`(0,1]`
	   * - hf
	     - float
	     - Proportion of hits in the query.
	     - :math:`(0,1]`
	   * - c1
	     - integer
	     - Number of covered positions in the query.
	     - :math:`\{k, \ldots, |query|\}`
	   * - c2
	     - float
	     - Normalized number of covered positions in the query.
	     - :math:`(0,1]`
	   * - cf
	     - float
	     - Proportion of covered positions in the query.
	     - :math:`(0,1]`
	   * - is
	     - int
	     - Number of reported assignments (nodes) for the query.
	     - :math:`\{1, \ldots, |leaves|\}`
	   * - ii
	     - int
	     - ID of the curent assignment.
	     - :math:`\{1, \ldots, is\}`
	   * - hc
	     - string
	     - Hit bit-mask encoded as a CIGAR string. For instance, `7=1X3=` means `11111110111`.
	     -

|

Read assignments in a Kraken-like format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Introduction
""""""""""""

ProPhyle uses a format similar to the `Kraken output <https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format>`_ for reporting k-mer matches by `ProPhyle Index <https://github.com/prophyle/prophyle_index>`_. It can also use this format
for reporting the final assignments.


Specification
"""""""""""""

	.. list-table:: Kraken-like format
	   :widths: 3 25
	   :header-rows: 1

	   * - Column
	     - Description
	   * - 1
	     - C / U (classified / unclassified)
	   * - 2
	     - Query name
	   * - 3
	     - Final assignments – a comma separated list of node names
	   * - 4
	     - Query length
	   * - 5
	     - K-mer mappings: a space-delimited lists of mappings. A single mapping is of the form ``comma_delimited_list_of_nodes:length``. Pseudo-nodes ``A`` and ``0`` are used for k-mers with a non-``ACGT`` nucleotide and without any mapping, respectively.



Examples
""""""""

Assigned k-mers, no sequences::

	U	read3	0	8	left,right:1 A:3 0:1 right:1


Assigned k-mers, version with sequences and base qualities::

	U	read3	0	8	left,right:1 A:3 0:1 right:1	CTTNGTTT	IGIIIIHI











Abundances estimates (experimental)
-----------------------------------

.. index:: Kraken report

Abundances in the Kraken report format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Introduction
""""""""""""

Specification
"""""""""""""

`kraken-report <https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports>`_ format:


	.. list-table:: Kraken report format
	   :widths: 5 20
	   :header-rows: 1

	   * - Column
	     - Description
	   * - 1
	     - Percentage of reads covered by the clade rooted at this taxon
	   * - 2
	     - Number of reads covered by the clade rooted at this taxon
	   * - 3
	     - Number of reads assigned directly to this taxon
	   * - 4
	     - A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
	   * - 5
	     - NCBI taxonomy ID
	   * - 6
	     - Indented scientific name



.. index:: MetaPhlAn report

Abundances in the MetaPhlAn2 report format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Introduction
""""""""""""

`MetaPhlAn2 <http://huttenhower.sph.harvard.edu/metaphlan2>`_ is a computational tool for profiling the composition of microbial communities from metagenomic sequencing data.


Specification
"""""""""""""

`MetaPhlAn2 report format <https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2#rst-header-output-files>`_

	.. list-table:: Metaphlan 2 report format
	   :widths: 5 20
	   :header-rows: 1

	   * - Column
	     - Description
	   * - 1
	     - Clades, ranging from taxonomic kingdoms (Bacteria, Archaea, etc.) through species
	   * - 2
	     - The taxonomic level of each clade is prefixed to indicate its level: Kingdom: ``k__``, Phylum: ``p__``, Class: ``c__``, Order: ``o__``, Family: ``f__``, Genus: ``g__``, Species: ``s__``



Since sequence-based profiling is relative and does not provide absolute cellular abundance measures, clades are hierarchically summed. Each level will sum to 100%; that is, the sum of all kindom-level clades is 100%, the sum of all genus-level clades (including unclassified) is also 100%, and so forth. OTU equivalents can be extracted by using only the species-level ``s__`` clades from this file (again, making sure to include clades unclassified at this level).


.. index:: Centrifuge report

Abundances in the Centrifuge report format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Introduction
""""""""""""

`Centrifuge <https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-summary-output-the-default-filename-is-centrifuge_report.tsv>`_ format.

Specification
"""""""""""""

	.. list-table:: Centrifuge format
	   :widths: 5 20
	   :header-rows: 1


	   * - Column
	     - Description
	   * - 1
	     - name of a genome, or the name corresponding to a taxonomic ID (the second column) at a rank higher than the strain (e.g., Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis).
	   * - 2
	     - taxonomic ID (e.g., 36870).
	   * - 3
	     - taxonomic rank (e.g., leaf).
	   * - 4
	     - number of k-mers propagated up to the node (e.g., 703004).
	   * - 5
	     - number of reads classified to this node including multi-classified reads (divided by the number of assignments, e.g., 5981.37)
	   * - 6
	     - number of reads uniquely classified to this genomic sequence (e.g., 5964)
	   * - 7
	     - unused


Example
"""""""

::

	#name                                                           taxID   taxRank    kmerCount   numReads   numUniqueReads   abundance
	Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis 36870   leaf       703004      5981.37    5964             0









Internal ProPhyle formats
-------------------------

.. index:: ProPhyle index


ProPhyle Index
^^^^^^^^^^^^^^

Introduction
""""""""""""

ProPhyle index directory contains a BWA index,
a k-LCP array and several auxiliary files.


Specification
"""""""""""""

	.. list-table:: ProPhyle index
	   :widths: 5 20
	   :header-rows: 1

	   * - File name
	     - Description
	   * - ``index.fa``
	     - Assembled contigs, name of sequences are of the following format: ``[node_name]@c[contig_id]``
	   * - ``index.fa.amb``
	     - List of ambiguous nucleotides, no values
	   * - ``index.fa.ann``
	     - List of contigs and their starting positions in the master string
	   * - ``index.fa.[k].klcp``
	     - k-LCP array
	   * - ``index.fa.bwt``
	     - Burrows-Wheeler Transform of the master string (merged sequences + reverse completement) + OCC table (BWA format)
	   * - ``index.fa.kmers.tsv``
	     - k-mer statistics, format: ``[node_name].[full|reduced].fa	[#kmers]``, where ``full`` refers to all associated k-mers and ``reduced`` to represented k-mers
	   * - ``index.fa.pac``
	     - Packed sequences (BWA format)
	   * - ``index.fa.sa``
	     - Sampled suffix array (BWA format)
	   * - ``index.json``
	     - Index parameters: k-mer size (``k``), ProPhyle version (``prophyle-version``, ``prophyle-revision``, ``prophyle-commit``)
	   * - ``log.txt``
	     - Log
	   * - ``tree.nw``
	     - Phylogenetic tree adjusted for classification
	   * - ``tree.preliminary.nw``
	     - Phylogenetic tree before adjusting



.. index:: ProPhyle compressed index

Compressed ProPhyle index for transmission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Introduction
""""""""""""

ProPhyle can create a ``.tar.gz`` archive with the a subset of the index files so
that the original index can be derived.

Specification
"""""""""""""

The archive contains the following subset of the original index files:

	.. list-table:: Compressed ProPhyle index
	   :widths: 5 20
	   :header-rows: 1

	   * - File name
	     - Description
	   * - ``index.fa.amb``
	     - Identical
	   * - ``index.fa.ann``
	     - Identical
	   * - ``index.fa.bwt``
	     - Burrows-Wheeler Transform *without* the OCC table (BWA format, before ``bwa bwtupdate``)
	   * - ``index.fa.kmers.tsv``
	     - Identical
	   * - ``index.json``
	     - Identical
	   * - ``tree.nw``
	     - Identical
	   * - ``tree.preliminary.nw``
	     - Identical
