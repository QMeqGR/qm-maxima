#!/bin/bash

# script version 2.0

# Note: When maxima runs a batch file, the output contains the command line
# but without the last 'termination' character. So ; and $ are lost. Thus we
# have to do it here by saving those characters and re-inserting them. This
# is done below in the function restore_lastchar().

packname=""
debug=0


# Section A
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


function restore_lastchar () {
    count=0;
    macfile=$1;
    tmpfile=$2;
    n_macfile=$(wc -l $macfile | awk '{print $1}');
    n_tmpfile=$(wc -l $tmpfile | awk '{print $1}');
    cp $tmpfile .tmp.recurse;
    if [ $debug -ge 1 ]; then
	echo "mac file lines: "$n_macfile;
	echo "tmp file lines: "$n_tmpfile;
    fi
    # get last characters in the mac file, should be ; or $
    cat $macfile | awk '(NR>3){print $0}' |
	sed 's/[[:space:]]*$//;s/.*\(.\)/\1/' > .tmp.lastchar
    
    for line in $(seq 1 $n_tmpfile); do
	is_input=$(cat $tmpfile | awk -v LN=$line \
				      '(NR==LN){if($1~/\(%i[1-9]/){print 1}else{print 0}}');
	if [ $debug -ge 1 ]; then 
	    echo "is_input: "$is_input
	fi
	if [ $is_input -ge 1 ]; then
	    count=$((count+1)); # increment count and get this line from the lastchar file
	    lc=$(cat .tmp.lastchar | head -$count | tail -1);
	    if [ $debug -ge 1 ]; then
		echo "------------> last char lc is \"$lc\"        input line count: "$count "   line number: "$line
		echo "              input line: "$(cat $tmpfile | head -$line | tail -1)
	    fi
	    cat .tmp.recurse | awk -v LC=$lc -v LN=$line \
		   '(NR==LN && ($1~/\(%i[1-9]/)){printf("%s%s\n",$0,LC)}(NR!=LN){print $0}' > .tmp.junk
	    mv -f .tmp.junk .tmp.recurse
	fi
    done
    cat .tmp.recurse > .tmp.grp.$i.tmp4
}

###########################################
# awk script for generating markers
function make_awkfile (){
cat > tmp.awk <<EOF
BEGIN{
    debug=0;
    new_example=0;
    excnt=0;
    gcount=0;
    icount=0;
}

(\$1~/@example/){new_example=1; excnt++; print \$0;}
(\$1~/@group/){new_group=1; gcount++; print \$0;} 
(\$1~/@end/ && \$2~/example/){new_example=0; print \$0;}
(\$1~/@end/ && \$2~/group/){new_group=0; print \$0;}
(\$1~/\(%i[1-9]/ && new_example){
    icount++;
    if (debug) {
        printf("%s",substr(\$0,length(\$1)+1));
        printf("\n");
    }
    maxcom=substr(\$0,length(\$1)+1);
    printf("%%grpcnt= %d rgen-icount= %d -->>%s\n",gcount,icount,maxcom);
 }
(\$1 !~ /\(%i[1-9]/ &&
 \$1 !~/\(%o[1-9]/ &&
 (\$1 !~/@end/ && \$2 !~/example/ ) &&
 new_example==0 && new_group==0){print \$0}
(\$1 ~ /@end/ && \$2 !~ /example/ && \$2 !~ /group/ ){print \$0}

EOF
}

make_awkfile;

# The markers look like:
#                      %grpcnt= 2 rgen-icount= 3 -->> rvec(1,2,3);

# clean out examples and insert markers
cat $packname.texi | awk -f ./tmp.awk > tmp.regen.1;

# generate markers list
icount_list=$(cat tmp.regen.1 | awk '($1~/grpcnt/){printf("%s ",$4);}');
grpcnt_list=$(cat tmp.regen.1 | awk '($1~/grpcnt/){printf("%s ",$2);}');

if [ $debug -gt 0 ]; then
    echo "icount_list= "$icount_list
    echo
    echo "grpcnt_list= "$grpcnt_list
    echo
fi


# Section B
# All commands in the same group must be run as one batch file.

###################################################
# Loop over the groups and make the batch files
ngrps=$(echo $grpcnt_list | awk '{print $NF}')
if [ $debug -gt 0 ]; then echo "ngrps= "$ngrps; fi
for i in $(seq 1 $ngrps); do
    echo "display2d_unicode:false$" > .tmp.grp.$i.mac
    echo "load($packname)$" >> .tmp.grp.$i.mac
    echo "linenum:0$" >> .tmp.grp.$i.mac
    echo "################ extracting examples from group "$i "####################";
    cat tmp.regen.1 | awk -v G="$i" --source \
    '($1~/grpcnt/ && $2==G){for(i=6;i<NF+1;i++){printf("%s ",$i);}printf("\n");}' >> .tmp.grp.$i.mac;
    # now we must capture the last character of the input commands to append in the final output
    # now batch run the file
    maxima -q --batch .tmp.grp.$i.mac > .tmp.grp.$i.tmp1
    # throw out any tellsimp warnings
    cat .tmp.grp.$i.tmp1 | awk '($1 !~ "tellsimp" && $2 !~ "warning"){print $0}' > .tmp.grp.$i.tmp2
    # throw away the header and the last line
    nlines=$(wc -l .tmp.grp.$i.tmp2 | awk '{print $1}')
    # find which line contains 'linenum:0'
    linenumline=$(grep -n linenum .tmp.grp.$i.tmp2 | awk -F: '{print $1}')
    cat .tmp.grp.$i.tmp2 | awk -v N="$nlines" -v NLL="$linenumline"\
			       --source '(NR>NLL && NR<N){print $0}' > .tmp.grp.$i.tmp3
    # Now post process the output to include the lastcharacter from above...
    restore_lastchar .tmp.grp.$i.mac .tmp.grp.$i.tmp3;
    cat .tmp.grp.$i.tmp4 > .tmp.grp.$i.out
    #cat .tmp.grp.$i.tmp3 | awk '{if($1~/\(%i[1-9]/){printf("%s;\n",$0,$LC)}else{print $0}}' > .tmp.grp.$i.out
done

####################################################
# Insert the output in the appropriate places
# Loop over the groups
for i in $(seq 1 $ngrps); do
    if [ $i -eq 1 ]; then
	# spit out the file up to the first group
	cat tmp.regen.1 | awk --source \
	   'BEGIN{stop=0;}($1~/grpcnt/ && $2==1){exit;}(stop==0){print $0}' > tmp.regen.2;
	cat .tmp.grp.1.out >> tmp.regen.2;
    else
	# append for markers between i-1 and i
	im=$((i-1));
	cat tmp.regen.1 | awk -v I="$i" -v IM="$im" --source 'BEGIN{stop=1;}
					($1~/grpcnt/ && $2==IM){stop=0;}
					($1~/grpcnt/ && $2==I){stop=1;}
					(stop==0){print $0}' >> tmp.regen.2;
	cat .tmp.grp.$i.out >> tmp.regen.2;
    fi
    # spit out the rest of the file after the last group
    if [ $i -eq $ngrps ]; then
	cat tmp.regen.1 | awk -v LAST="$ngrps" --source \
			       'BEGIN{stop=1;}
			       ($1~/grpcnt/ && $2==LAST){stop=0;}
			       (stop==0){print $0}' >> tmp.regen.2;
    fi
done

####################################################
# Strip the markers from the file for final output
cat tmp.regen.2 | awk '{if($1!~/grpcnt/){print $0}}' > regen.texi;

echo "##### Output is in: regen.texi"
echo "#####"
echo "##### Warning: Check this file carefully with 'diff' before"
echo "##### you replace the original .texi file!!!"
echo "#####"
echo "##### try: diff $packname.texi regen.texi | less"
echo "#####"

############
# Clean up
if [ $debug -eq 0 ]; then
    rm -f .tmp.* tmp.regen.1 tmp.regen.2 tmp.awk
fi
