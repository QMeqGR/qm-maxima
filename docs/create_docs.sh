#!/bin/bash

buildindex=/home/packages/SOURCE/maxima-code/doc/info/build_index.pl

packname=qm

makeinfo $packname.texi;
texi2pdf $packname.texi;
texi2html $packname.texi;
texi2any --plaintext $packname.texi > ../README.md;

# build the index
$buildindex $packname.info > $packname-index.lisp;

rm $packname.aux $packname.fn $packname.fns $packname.log $packname.toc
