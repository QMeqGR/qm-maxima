#!/bin/bash

buildindex=/home/packages/SOURCE/maxima-code/doc/info/build_index.pl

packname=qm

makeinfo $packname.texi;
makeinfo --pdf $packname.texi;
makeinfo --html $packname.texi;
makeinfo --plaintext $packname.texi > ../README.md;

# build the index
$buildindex $packname.info > $packname-index.lisp;

rm $packname.aux $packname.fn $packname.fns $packname.log $packname.toc

## Create test suite from examples in $packname.texi
echo "#### creating examples.txt and the rtest file."
echo "display2d:false$" > examples.txt
echo "load($packname)\$" >> examples.txt
cat $packname.texi | grep "(%i" | \
    gawk '(NF>1){for(i=2;i<NF+1;i++){printf("%s",$i)};printf("\n")}' \
	 >> examples.txt;

echo "#### Running maxima on tests..."
maxima -q -b examples.txt > rtest.tmp.out;

echo "#### Constructing rtest file... the rtest file must be edited"
cat rtest.tmp.out | gawk '($1~/%i/){for(i=2;i<NF+1;i++){printf("%s",$i)};printf(";\n");}($1~/%o/){for(i=2;i<NF+1;i++){printf("%s ",$i)};printf("$\n\n");}' > ../rtest_$packname.mac;


