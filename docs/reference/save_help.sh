#! /usr/bin/env bash

#set -e -f -o pipefail

cd $(dirname $0)

for x in main $(../../prophyle/prophyle --help | grep -E "^   " | perl -pe 's/^\s+//g' | cut -f1 -d ' ') ; do
	fn="$x.txt"
	echo "$fn"

	if [[ "$x" == "main" ]] ; then
		x=""
	fi

	echo "$ prophyle $x -h" > $fn
	echo >> $fn
	../../prophyle/prophyle.py $x -h >> $fn
done


# prophyle_ncbi_tree.py

fn="prophyle_ncbi_tree.py.txt"
com="../../prophyle/prophyle_ncbi_tree.py"

echo "$fn"
echo "$ `basename $com` -h" > $fn
echo >> $fn
$com -h 2>&1 >>$fn

# prophyle_index

fn="prophyle_index.txt"
com="../../prophyle/prophyle_index/prophyle_index"

echo "$fn"
echo "$ `basename $com` -h" > $fn
echo >> $fn
$com 2>> $fn

for sub_com in build query ; do
	sub_fn="${fn}_${sub_com}.txt"
	echo $sub_fn
	echo "$ `basename $com` $sub_com -h" > $sub_fn
	echo >> $sub_fn
	$com $sub_com 2>> $sub_fn
done

# prophyle_assembler

fn="prophyle_assembler.txt"
com="../../prophyle/prophyle_assembler/prophyle_assembler"

echo "$fn"
echo "$ `basename $com` -h" > $fn
echo >> $fn
$com -h 2>> $fn

# prophyle_assignment.py

fn="prophyle_assignment.py.txt"
com="../../prophyle/prophyle_assignment.py"

echo "$fn"
echo "$ `basename $com` -h" > $fn
echo >> $fn
$com -h >> $fn

# prophyle_assignment

fn="prophyle_assignment.txt"
com="../../prophyle/prophyle_assignment/prophyle_assignment"

echo "$fn"
echo "$ `basename $com`" > $fn
echo >> $fn
$com 2>> $fn

# prophyle_analyze

fn="prophyle_analyze.py.txt"
com="../../prophyle/prophyle_analyze.py"

echo "$fn"
echo "$ `basename $com` -h" > $fn
echo >> $fn
$com -h 2>&1 >>$fn
