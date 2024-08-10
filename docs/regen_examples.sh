#!/bin/bash

packname=qm;

cat $packname.texi | gawk -f ./regen_examples.awk > tmp.regen.1;

icount_list=$(cat tmp.regen.1 | gawk '($1~/rgen-icount/){printf("%s ",$2);}');
comm_list=$(cat tmp.regen.1 | \
		gawk '($1~/rgen-icount/){for(i=4;i<NF+1;i++){printf("%s ",$i);}printf("\n");}')

echo $icount_list
echo $comm_list

#maxima --batch-string="load(qm); brap(bra(a,g));"

