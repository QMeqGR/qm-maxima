#!/bin/bash

MAXIMA_ROOT=/home/packages/SOURCE/maxima-code
buildindex=$MAXIMA_ROOT/doc/info/build_index.pl

packname=qm

makeinfo $packname.texi;
makeinfo --pdf $packname.texi;
makeinfo --split=chapter --no-node-files --html \
	 -c OUTPUT_ENCODING_NAME=UTF-8 -e 10000 $packname.texi;
makeinfo --plaintext $packname.texi > ../README.md;

# build the .info index
$buildindex $packname.info > $packname-index.lisp;

# build the html index
maxima --no-init --no-verify-html-index  \
       --preload=$MAXIMA_ROOT/doc/info/build-html-index.lisp \
       --batch-string='build_and_dump_html_index("./qm_html/*.html", output_file="package-index-html.lisp",truenamep=true);';

if [ -f "package-index-html.lisp" ]; then
    echo "### Creation of HTML docs successful."
    mv -f package-index-html.lisp $packname-index-html.lisp;
else
    echo "### Warning: no maxima-index-html.lisp was created for html docs."
fi

# clean up temporary files from pdf creation
rm $packname.aux $packname.fn $packname.fns $packname.log $packname.toc


####################################################
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


