#!/bin/bash

debug=0;

packname=$(cat ../CONFIG | gawk -F'=' '($1=="package-name"){print $2}')
MAXIMA_ROOT=$(cat ../CONFIG | gawk -F'=' '($1=="max_src"){print $2}')
buildindex=$MAXIMA_ROOT/doc/info/build_index.pl

echo "### Running makeinfo..."
makeinfo $packname.texi;
echo "### Running makeinfo --pdf ..."
makeinfo --pdf $packname.texi;
echo "### Running makeinfo --html ..".
makeinfo --split=chapter --no-node-files --html \
	 -c OUTPUT_ENCODING_NAME=UTF-8 -e 10000 $packname.texi;
echo "### Running makeinfo --plaintext ..."
makeinfo --plaintext $packname.texi > ../README.txt;

echo "### Building .info index file..."
# build the .info index
$buildindex $packname.info > $packname-index.lisp;

# build the html index
echo "### Building .info html index file..."
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
if [ $debug -eq 0 ]; then
    rm $packname.aux $packname.fn $packname.fns $packname.log $packname.toc \
       $packname.vr $packname.vrs $packname.texi.toc build-html-index.log 
fi

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


