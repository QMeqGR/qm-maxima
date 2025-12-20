#!/bin/bash

# defaults
debug=0;
help=0;
MAXIMA=maxima;
MAXIMA_SRC=""
packname=""

###############################################
###############################################
declare SWITCH
while getopts "Dhm:M:p:" SWITCH; do
    case $SWITCH in
	D) debug=1 ;;
	h) help=1 ;;
	m) MAXIMA=$OPTARG;;
	M) MAXIMA_SRC=$OPTARG;;
	p) packname=$OPTARG;;
    esac
done

if [ $help -eq 1 ]; then
    echo "Usage:"
    echo "create_docs.sh [-hD]  [-p package-name] [-m maxima] [-M maxima-src]"
    echo 
    echo "    -D --- debug (leaves temp files, default is OFF)"
    echo "    -h --- help (show this help)"
    echo "    -p --- set package name (default off, automatically determined)"
    echo "    -m --- maxima executable (default maxima)"
    echo "    -M --- maxima source root directory (automatically determined)"
    echo "           (something like /home/src/maxima-code if you cloned from"
    echo "            SourceForge, and if you build out-of-tree)"
    echo
    exit
fi

# get package name
n=$(ls *.texi | wc -l | awk '{print $1}');
if [ $n -gt 1 ]; then
    echo "Found more than one .texi file. There should be only one PKG.texi file."
    if [ -e regen.texi ]; then
	echo "Remove the regen.texi file before running create_docs.sh."
    fi
    exit
elif [ "$packname" == "" ]; then
    texifile=$(ls *.texi);
    packname=${texifile%.texi}
    echo "Found packname= "$packname;
fi

if [ ! -f "$packname.texi" ]; then
    echo "$packname.texi is not a file"
    exit 2
fi

if [ "$MAXIMA_SRC" == "" ]; then
    MAXIMA_SRC=$($MAXIMA -d | awk -F'=' '($1~/maxima-srcdir/){print $2}');
fi

if [ $debug -gt 0 ]; then
    echo "MAXIMA="$MAXIMA
    echo "MAXIMA_SRC="$MAXIMA_SRC
fi

buildindex=$MAXIMA_SRC/doc/info/build_index.pl
buildindex_lsp=$MAXIMA_SRC/doc/info/build-html-index.lisp

if [ ! -e $buildindex ]; then
    echo "Can't find build_index.pl using"
    echo "MAXIMA_SRC="$MAXIMA_SRC
    echo "Try using -M switch to set src directory manually."
    exit 2
fi


#################################################################
#################################################################
makeinfo $packname.texi;
makeinfo --pdf $packname.texi;
makeinfo --split=chapter --no-node-files --html \
	 -c OUTPUT_ENCODING_NAME=UTF-8 -e 10000 $packname.texi;
makeinfo --plaintext $packname.texi > ../README.txt;
# build the .info index
$buildindex $packname.info > $packname-index.lisp;

# build the html index
DIR="${packname}"'_html/*.html'
OUTPUT="${packname}"-index-html.lisp
BATCH=`echo build_and_dump_html_index\(\"$DIR\", output_file=\"$OUTPUT\", truenamep=true\)\;`

$MAXIMA --no-init --no-verify-html-index  \
       --preload=$buildindex_lsp \
       --batch-string="$BATCH"

# clean up temporary files from pdf creation
if [ $debug -eq 0 ]; then
    rm $packname.aux $packname.fn $packname.fns $packname.log $packname.toc \
       $packname.vr $packname.vrs build-html-index.log 
fi

echo "Generating test suite from examples..."
####################################################
## Create test suite from examples in $packname.texi
echo "(display2d:false,kill(all),load($packname),0);" > examples.txt
cat $packname.texi | grep "(%i" | \
    awk '(NF>1){for(i=2;i<NF+1;i++){printf("%s ",$i)};printf("\n"q);}' \
	 | grep -v "load(\"$packname.mac\")" >> examples.txt;

$MAXIMA -q -b examples.txt > rtest.tmp.out;

cat rtest.tmp.out | grep -v "examples.txt" | grep -v "Loading" \
    | awk '{for(i=2;i<NF+1;i++){printf("%s ",$i);printf(";\n");}}'\
	  > rtest_example_$packname.mac;

if [ $debug -eq 0 ]; then
    rm -f examples.txt rtest.tmp.out rtest.tmp2.out
fi
