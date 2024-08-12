#!/bin/bash

maxima_srcdir=/home/packages/SOURCE/maxima-code

maxima --no-init --no-verify-html-index  \
	--preload=$maxima_srcdir/doc/info/build-html-index.lisp \
 	--batch-string='build_and_dump_html_index("./*.html");'

# maxima --no-init --batch-string="quit();"
