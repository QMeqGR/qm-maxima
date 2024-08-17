########################################################
#        Maxima 3rd-party package instructions         #
########################################################

These instructions assume you have downloaded the package template.

This template holds skeleton documentation files for your package as
well as skeleton .mac and .lisp files for your source code.

The top-level directory will be named after your package, say "PKG".
Inside this directory you will find the source code template files
PKG.mac and PKG.lisp, as well as a file called CONFIG that will set
some parameters. In the top-level template directory you will find a
doc/ directory that will contain your PKG.info doc file and two
scripts that will help you build your documentation.

For your package documentation to be visible once you 'load()' the
package you must create the documentation in the doc/ folder using a
script file provided there. Detailed instructions follow.

For the documentation to work correctly you need to have the Maxima
source installed on your computer because some of the utilities that
generate the documentation index are only contained within the Maxima
source.

Steps:

1. Name your package. For this howto let's assume your package name is
"PKG". Edit the file named CONFIG in the top level directory to
contain your package name, and also point to the location of your
Maxima source code install. You will set the variables 'package-name'
and 'max_src' in CONFIG. You must do this for the scripts to
work. Your CONFIG file should look like (no spaces are allowed in your
package name, use an underscore or dash if you must):

$ cat CONFIG
package-name=PKG
max_src=/home/packages/SOURCE/maxima-code

2. Edit the LICENSE file to whatever license you like. The template
contains a GPL-v3 license.

3. You should name (or rename) your top level directory
"PKG-maxima". Now you can begin working on your package's source
code. Rename the source files PKG.mac and PKG.lisp appropriately. Put
all of your source code in these files. You may want additional .mac
or .lisp source files. If you do you should load all of them with
'load()' commands from the main PKG.mac file. Note: When you issue
'load(PKG)' it will only load the PKG.mac file, so if you have
PKG.lisp you need to load it from PKG.mac with a line like
'load("PKG.lisp")'.

4. In your PKG-maxima/doc/ directory rename PKG.texi appropriately and
edit the file. Heed the warnings in this file about keeping the
structure as it is, otherwise 'makeinfo' is likely to fail.

Examples in your PKG.texi must have the following format that must be
followed exactly for the scripts to work.

@example
@group
(%i1) command here;
(%i2) another command here;
@end group
@end example

If your example uses the output of one command for the next input
command, then both commands MUST be in the same '@group'.

Important note: You do *not* have to cut and paste output from a
Maxima session into your PKG.texi file. It is fine for your .texi file
to contain only the input commands for your examples. The output will
be generated automatically using the script file
"regen_examples.sh".

Once you are finished editing your PKG.texi file, you can generate all
of the examples inline by running:

$ ./regen_examples.sh

The script will process the different 'groups' in PKG.texi and place
the entire output of PKG.texi with the examples in
"regen.texi". WARNING: Before you overwrite your original PKG.texi
file YOU SHOULD DO:

$ diff PKG.texi regen.texi

and make sure that the only lines getting replaced are from the
examples.  If your PKG.texi file is missing "@end group" then
'regen_examples' will not process your .texi file correctly. BE
CAREFUL before you overwrite your PKG.texi file with the regen.texi
file!!

5. Now you are ready to generate the doc files. In the PKG-maxima/doc/
directory issue:

$ ./create_docs.sh

This should create PKG.info, PKG.pdf, a PKG_html/ directory containing
html-related doc files, a text file ../README.txt in the top level
directory, as well as some other index files.

From the main PKG.texi file it will also generate a file
../rtest_PKG.mac containing Maxima tests using the examples you have
in your PKG.texi file.

=============================== Done.

If everything worked correctly you should be able to do

(%i1)load(PKG);

and see the docs with ? and ??.

Troubleshooting:

1. Do you have your maxima path set up so that it can find your
package when you issue 'load(PKG);'?

2. Did everything build correctly? You can manually check the
*index.lisp index files, and the README.txt file.
