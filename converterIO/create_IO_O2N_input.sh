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


### Function check_field_geombc
function check_field_geombc() {
# All the critical fields that phasta needs in geombc are listed below. Update this list if needed for your application. Do not forget to update N_geombc_fields_double etc below
# The format is the following:  geombc, field name, double or integer, block or header, number of integer in the header after < >. 

        #field_fun="$@" # get all args
        #field_fun=$(echo $argfun | awk '{print $1}') #Name of the field
	#filefun=$1
	file_double_field_geombc_fun=$1
	file_integer_field_geombc_fun=$2
	field_fun=$3
	list_interior_tpblocks_fun=$4
	list_boundary_tpblocks_fun=$5

	### Double fields first (compulsary for the converter for now)

        teststring="co-ordinates"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, $field_fun, double, block, 2;" >> $file_double_field_geombc_fun
        fi

        teststring="boundary condition array"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, $field_fun, double, block, 1;" >> $file_double_field_geombc_fun
        fi

	# Trying to keep the order identical. 
	# After 'boundary condition array' (already treated above) comes usually 'nbc values', which is the last double field
        teststring="boundary condition array"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
	  while read line
	  do 
	  	field=`echo $line | awk '{$1="";$2="";$3="";print $0}'`

		echo "geombc, nbc values $field, double,   block, 8;" | sed -e 's/  */ /g'  >> $file_double_field_geombc_fun # sed remove extra space
	  done <  $list_boundary_tpblocks_fun
	fi

	### Integer fields next

        teststring="number of interior tpblocks"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, $field_fun, integer, header, 1;" >> $file_integer_field_geombc_fun
        fi

	# This will include all the interior topologies
        #teststring="connectivity interior"
        #if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
        #  echo "geombc, "$field_fun", integer, block, 7;" >> $file_integer_field_geombc_fun
        #fi

        teststring="number of interior elements"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, header, 1;" >> $file_integer_field_geombc_fun
        fi

	# This will include 
	# "number of nodes in the mesh", 
	# "number of nodes" 
	# "number of nodes with Dirichlet BCs" 
        teststring="number of nodes"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, header, 1;" >> $file_integer_field_geombc_fun
        fi
        
        teststring="maximum number of element nodes"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, header, 1;" >> $file_integer_field_geombc_fun
        fi

        teststring="number of modes"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, header, 1;" >> $file_integer_field_geombc_fun
        fi

        teststring="bc mapping array"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, block, 1;" >> $file_integer_field_geombc_fun
        fi
        
        teststring="bc codes array"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, block, 1;" >> $file_integer_field_geombc_fun
        fi

        teststring="number of boundary tpblocks"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, header, 1;" >> $file_integer_field_geombc_fun
        fi

        #teststring="connectivity boundary"
        #if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
        #  echo "geombc, "$field_fun", integer, block, 8;" >> $file_integer_field_geombc_fun
        #fi

        #teststring="nbc codes"
        #if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
        #  echo "geombc, "$field_fun", integer, block, 8;" >> $file_integer_field_geombc_fun
        #fi

        teststring="number of boundary elements"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, header, 1;" >> $file_integer_field_geombc_fun
        fi

        teststring="size of ilwork array"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, header, 1;" >> $file_integer_field_geombc_fun
        fi

        teststring="ilwork"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, block, 1;" >> $file_integer_field_geombc_fun
        fi

        teststring="periodic masters array"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "geombc, "$field_fun", integer, block, 1;" >> $file_integer_field_geombc_fun
        fi

        # NSpre and others produces a field named "number of shapefunctions soved on processor"
        # Phasta-SyncIO requires now a field named "number of shape functions@partID" read from readnblk.f. 
	# Therefore, rename this field accordingly here
	# Beware: in some version of NSpre, the historic spelling mistake "soved" has been replaced by "solved"
        teststring1="number of shapefunctions soved on processor"
        teststring2="number of shapefunctions solved on processor"
        if [ "${field_fun:0:${#teststring1}}" == "$teststring1" ] || [ "${field_fun:0:${#teststring2}}" == "$teststring2" ]; then
          echo "geombc, number of shape functions, integer, header, 1;" >> $file_integer_field_geombc_fun
	  echo "INFO: 'number of shapefunctions soved on processor' is renamed 'number of shape functions' for readnblk.f"
        fi

	# Trying to keep the order identical. 
	# For phParAdapt, after 'number of nodes in the mesh' (already treated above) comes usually 'connectivity interior'.
        # But there is no such field generated with NSpre so we use instead "number of shapefunctions soved on processor" which is common to both.
        # Moreover, in readnblk.f, the connectivity is read just after "number of shapefunctions soved on processor" so this ordering is coherent.
	# NOTE that we also add a new field called 'total number of interior tpblocks'. This field will be saved in the new syncIO geombc files.
        #teststring="number of nodes in the mesh"
        teststring1="number of shapefunctions soved on processor"
        teststring2="number of shapefunctions solved on processor"
        if [ "${field_fun:0:${#teststring1}}" == "$teststring1" ] || [ "${field_fun:0:${#teststring2}}" == "$teststring2" ]; then
	  
	  echo "geombc, total number of interior tpblocks, integer, header, 1;" >> $file_integer_field_geombc_fun
	  echo "INFO: a new field called 'total number of interior tpblocks' will be added in the new geombc files"

	  while read line
	  do 
	  	field=`echo $line | awk '{$1="";$2="";$3="";print $0}'`
		echo "geombc, connectivity interior $field, integer, block, 7;" | sed -e 's/  */ /g'  >> $file_integer_field_geombc_fun # sed remove extra space
	  done <  $list_interior_tpblocks_fun
        fi

	# Trying to keep the order identical. 
	# After 'number of nodes in the mesh' (already treated above) comes usually 'connectivity boundary' and 'nbc code'
	# NOTE that we also add a new field called 'total number of boundary tpblocks'. This field will be saved in the new syncIO geombc files.
        teststring="number of boundary tpblocks"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string

	  echo "geombc, total number of boundary tpblocks, integer, header, 1;" >> $file_integer_field_geombc_fun
	  echo "INFO: a new field called 'total number of boundary tpblocks' will be added in the new geombc files"

	  while read line
	  do 
	  	field=`echo $line | awk '{$1="";$2="";$3="";print $0}'`

	  	echo "geombc, connectivity boundary $field, integer, block, 8;" | sed -e 's/  */ /g'  >> $file_integer_field_geombc_fun # sed remove extra space
	  	echo "geombc, nbc codes $field, integer, block, 8;" | sed -e 's/  */ /g'  >> $file_integer_field_geombc_fun # sed remove extra space
	  done <  $list_boundary_tpblocks_fun
	fi
}


### Function check_field_restart
function check_field_restart() {
# All the critical fields that phasta needs in restart are listed below. Update this list if needed for your application.
# The format is the following:  restart, field name, double or integer, block or header, number of integer in the header after < >. 

        #field_fun="$@" # get all args
        #field_fun=$(echo $argfun | awk '{print $1}') #Name of the field
	#filefun=$1
	file_double_field_restart_fun=$1
	file_integer_field_restart_fun=$2
	field_fun=$3

	### Double fields first (compulsary for the converter for now)

        teststring="solution"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "restart, $field_fun, double, block, 3;" >> $file_double_field_restart_fun
        fi

	### Integer fields next

        teststring="byteorder magic number"
        if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
          echo "restart, $field_fun, integer, block, 1;" >> $file_integer_field_restart_fun
        fi

        #teststring="number of modes" # Already in the geombc files and read by phasta from the geombc files
        #if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
        #  echo "restart, $field_fun, integer, header, 1;" >> $file_integer_field_restart_fun
        #fi

        #teststring="number of variables" # Was required for the converter. Info now read from the solution field
        #if [ "${field_fun:0:${#teststring}}" == "$teststring" ]; then #we compare the first letters of the string
        #  echo "restart, $field_fun, integer, header, 1;" >> $file_integer_field_restart_fun
        #fi

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

if [ ! -e $dir/restart.$N_steps.1 ]; then
        die "Time step $N_steps does not exist in $N_parts-procs_case\n Aborting"
fi

resmodulo=$(($N_parts % $N_files))
if [ "$resmodulo" -ne "0" ]; then
        die "The number of SyncIO files requested $N_files is not a multiple of the number of parts $N_parts\n Aborting"
fi


### First, count the interior and boundary topology blocks
list_interior_tpblocks=list_interior_tpblocks.dat
grep -aH ' : < ' $dir/geombc.dat.* | grep 'connectivity interior' | awk -F : '{print $1,$2}' | awk '{$1=""; print $0}' | sort | uniq -c  > $list_interior_tpblocks
interior_tpblocks=`cat $list_interior_tpblocks | wc -l`
echo "There are $interior_tpblocks different interior tp blocks in all the geombc files:"
cat $list_interior_tpblocks
echo ""

list_boundary_tpblocks=list_boundary_tpblocks.dat
grep -aH ' : < ' $dir/geombc.dat.* | grep 'connectivity boundary' | awk -F : '{print $1,$2}' | awk '{$1=""; print $0}' | sort | uniq -c  > $list_boundary_tpblocks
boundary_tpblocks=`cat $list_boundary_tpblocks | wc -l`
echo "There are $boundary_tpblocks different boundary tp blocks in all the geombc files:"
cat $list_boundary_tpblocks
echo ""


### Grep all the fields from geombc.dat.1 (posix file)
file_grep_geombc=grep_geombc_posix.dat
grep -a ' : < ' $dir/geombc.dat.1 > $file_grep_geombc


### Get only the fields name from geombc.dat.1 and remove any extra unwanted space
file_field_geombc=field_geombc_posix.dat
cat $file_grep_geombc | awk -F : '{print $1}' |  sed -e 's/  */ /g' | sed -e 's/^[ \t]*//g' | sed -e 's/[ \t]*$//g' > $file_field_geombc


### Now check which fields are critical for phasta and save the double fields in $file_double_field_geombc and the integer fields in $file_integer_field_geombc for future use
file_double_field_geombc=file_double_field_geombc.dat
if [ -e $file_double_field_geombc ]; then
        rm $file_double_field_geombc
fi

file_integer_field_geombc=file_integer_field_geombc.dat
if [ -e $file_integer_field_geombc ]; then
        rm $file_integer_field_geombc
fi

while read line
do
        field=`echo $line`
        check_field_geombc "$file_double_field_geombc" "$file_integer_field_geombc" "$field" "$list_interior_tpblocks" "$list_boundary_tpblocks"
done <  $file_field_geombc

### Double check that the topologies have been correctly found 
interior_tpblocks2=`grep 'connectivity interior' $file_integer_field_geombc | wc -l`
if [ "x$interior_tpblocks" != "x$interior_tpblocks2" ]; then
	echo ""
	echo "ERROR: Missing connectivity interior, probably due to a spelling issue ($interior_tpblocks2 != $interior_tpblocks)"
        echo "ERROR: Check that \"number of shapefunctions soved on processor\" or \"number of shapefunctions solved on processor\""
	echo "ERROR: is present in your geombc files"
	echo ""
fi

boundary_tpblocks2=`grep 'connectivity boundary' $file_integer_field_geombc | wc -l`
if [ "x$boundary_tpblocks" != "x$boundary_tpblocks2" ]; then
	echo ""
	echo "ERROR: Missing connectivity boundary ($boundary_tpblocks2 != $boundary_tpblocks)"
	echo ""
fi

### Grep all the fields from restart.##.1 (posix file)
file_grep_restart=grep_restart_posix.dat
grep -a ' : < ' $dir/restart.$N_steps.1 > $file_grep_restart


### Get only the fields name from restart.##.1 and remove any extra unwanted space
file_field_restart=field_restart_posix.dat
cat $file_grep_restart | awk -F : '{print $1}' |  sed -e 's/  */ /g' | sed -e 's/^[ \t]*//g' | sed -e 's/[ \t]*$//g' > $file_field_restart


### Now check which fields are critical for phasta and save the double fields in $file_double_field_geombc and the integer fields in $file_integer_field_geombc for future use
file_double_field_restart=file_double_field_restart.dat
if [ -e $file_double_field_restart ]; then
        rm $file_double_field_restart
fi

file_integer_field_restart=file_integer_field_restart.dat
if [ -e $file_integer_field_restart ]; then
        rm $file_integer_field_restart
fi

while read line
do
        field=`echo $line`
        check_field_restart "$file_double_field_restart" "$file_integer_field_restart" "$field" 
done <  $file_field_restart

echo ""


### Start to write now IO.O2N.input
file=IO.O2N.input
if [ -e $file ]; then
	rm $file
fi

N_geombc_fields_double=`cat $file_double_field_geombc | wc -l`
N_geombc_fields_integer=`cat $file_integer_field_geombc | wc -l`
N_restart_fields_double=`cat $file_double_field_restart | wc -l`
N_restart_fields_integer=`cat $file_integer_field_restart | wc -l`

echo "N-geombc-fields-double: $N_geombc_fields_double;" >> $file
echo "N-geombc-fields-integer: $N_geombc_fields_integer;" >> $file
echo "N-restart-fields-double: $N_restart_fields_double;" >> $file
echo "N-restart-fields-integer: $N_restart_fields_integer;" >> $file
echo "N-steps: $N_steps;" >> $file
echo "N-parts: $N_parts;" >> $file
echo "N-files: $N_files;" >> $file

# The converter expects first the geombc double in IO.O2N.input, then geombc integer, then restart double, then restart integer.
cat $file_double_field_geombc >> $file
cat $file_integer_field_geombc >> $file
cat $file_double_field_restart >> $file
cat $file_integer_field_restart >> $file

### Some cleaning
rm $list_boundary_tpblocks
rm $list_interior_tpblocks
rm $file_grep_geombc
rm $file_grep_restart
rm $file_field_geombc
rm $file_field_restart
rm $file_double_field_geombc
rm $file_integer_field_geombc
rm $file_double_field_restart
rm $file_integer_field_restart

echo "$file generated for the converter"

