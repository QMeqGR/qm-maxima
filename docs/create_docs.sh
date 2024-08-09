#!/bin/bash

packname=qm

makeinfo $packname.texi;
texi2pdf $packname.texi;
texi2html $packname.texi;
texi2any --plaintext $packname.texi > ../README.md;

rm $packname.aux $packname.fn $packname.fns $packname.log $packname.toc
