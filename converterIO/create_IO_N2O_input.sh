#!/bin/bash

nargs=3  # Arguments of the script: 
         #  time step 
         #  number parts 
         #  number of syncio files wanted

### Function die called when there is a problem 
function die() {
        echo -e "${1}"
        exit 1
}


#####################################################################
### End of functions - Beginning of the script  
#####################################################################


### Check that user entered proper number of args
if [ "${#}" -lt "$nargs" ] || [ "$1" == "-h" ]; then
        die "\n Check usage: \n  $0\n\n \$1: <time step>\n \$2: <number of parts>\n \$3: <number of SyncIO files>\n\n"
fi


### Read arguments of the script
N_steps=$1
N_parts=$2
N_files=$3

echo "Time step: $N_steps"
echo "Number of parts: $N_parts" 
echo "Number of SyncIO files: $N_files"
echo ""

dir=$N_parts-procs_case


### Do a couple of sanity check on the input parameters of the script
if [ ! -d $dir ]; then
	die "$N_parts-procs_case does not exist\n Aborting"		
fi

if [ ! -e $dir/restart-dat.$N_steps.1 ]; then
        die "Time step $N_steps does not exist in $N_parts-procs_case\n Aborting"
fi

nsynciofiles=`ls $dir/restart-dat.$N_steps.* | wc -l`
if [ "$nsynciofiles" -ne "$N_files" ]; then
	die "The number of SyncIO files requested does not match with what is found in the procs_case directory: $nsynciofiles - $N_files"
fi

resmodulo=$(($N_parts % $N_files))
if [ "$resmodulo" -ne "0" ]; then
        die "The number of SyncIO files requested $N_files is not a multiple of the number of parts $N_parts\n Aborting"
fi


### First, count the number of fields (blocks) in the restart files and list them
list_restart_fields=list_restart_fields.dat
grep -a ' : < ' $dir/restart-dat.$N_steps.1 | grep -v '< 0 >' | awk -F @ '{print $1}' | uniq -c > $list_restart_fields
restart_fields=`cat $list_restart_fields | wc -l`
echo "There are $restart_fields different block fields in the restart files:"
cat $list_restart_fields
echo ""

### Start to write now IO.N2O.input
file=IO.N2O.input
if [ -e $file ]; then
	rm $file
fi

# By default, the geombc files already exist under the posix format so no need to deconvert them.
# Focuss only on the restart files, which should by default only contain double blocks.
# In what follows, any single header with no block (characterized by '< 0 >') are also ignored.
N_geombc_fields_double=0
N_geombc_fields_integer=0
N_restart_fields_double=$restart_fields
N_restart_fields_integer=0

echo "N-geombc-fields-double: $N_geombc_fields_double;" >> $file
echo "N-geombc-fields-integer: $N_geombc_fields_integer;" >> $file
echo "N-restart-fields-double: $N_restart_fields_double;" >> $file
echo "N-restart-fields-integer: $N_restart_fields_integer;" >> $file
echo "N-steps: $N_steps;" >> $file
echo "N-parts: $N_parts;" >> $file
echo "N-files: $N_files;" >> $file

# Now add all the fields found in the restart files.
# Asumption always verified so far: all the fields in the restart files are double and their header includes 3 integer.
# If this is not the case any more, then more complex verifications should be added like in the other script.

echo "By default, they will all be specified in $file in order to be deconverted but this can be manually modified"
echo "This script also assumes that all blocks listed are double"
while read line
do
        field=`echo $line | awk '{$1="";print $0}' |  sed -e 's/  */ /g' | sed -e 's/^[ \t]*//g' | sed -e 's/[ \t]*$//g'`
	echo "restart, $field, double, block, 3;" >> $file 
done <  $list_restart_fields

### Some cleaning
rm $list_restart_fields

echo "$file generated for the converter"

