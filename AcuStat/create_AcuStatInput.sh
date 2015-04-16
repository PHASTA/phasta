#!/bin/bash

nargs=6

function die() {
echo -e "${1}"
exit 1
}

# check that user entered proper number of args
if [[ "${#}" -ne "$nargs" ]]
then
die "\n Check usage: \n  $0\n <Number of cores (for procs_case)>\n <Starting time step for the average> \n Accumulate error field? <0/1>\n Accumulate phase average fields? <0/n_avg_fields> \n Copy vorticity field from the last time step? <0/1>\n Copy dwal field from the last time step? <0/1>"
fi

cores=$1
start_ts=$2
accu_error=$3
accu_phavg=$4
copy_vort=$5
copy_dwal=$6

outputfile='AcuStat_input.dat'


ls -lrt $cores-procs_case/restart-dat.*.1 | grep -v '^lr'> tmp
cat tmp  | awk -F dat '{print $2}' | sed 's/^\.//g' | sed 's/\.1//g' | sort -n > tmp2
rm tmp


if [ -f tmp3 ] ; then
	rm tmp3 
fi

iline=0
while  read ts; do
	iline=$((iline+1))		
	if [ "$iline" -gt "1" ] ; then
		diffts=$((ts-tsold))
		#echo $iline $ts $tsold $diffts
                if [ "$ts" -ge "$start_ts" ] ; then
			echo $ts $diffts >> tmp3
		fi 
        fi
        tsold=$ts
done < tmp2
rm tmp2

nlines=`wc -l tmp3 | awk '{print $1}'`

# Print first the info for the AcuStat executable
notice='# Number of restart files to read'
echo "$nlines $notice" > $outputfile
notice='# Accumulate error field? (0/1)'
echo "$accu_error $notice" >> $outputfile
notice='# Accumulate phase_average field? (0/n phase average fields)'
echo "$accu_phavg $notice" >> $outputfile
notice='# copy vorticity field? (0/1)'
echo "$copy_vort $notice" >> $outputfile
notice='# copy dwal? (0/1)'
echo "$copy_dwal $notice" >> $outputfile
cat tmp3 >> $outputfile 
rm tmp3

# Print now the usage of the script
echo ' ' >> $outputfile
echo '### What follows is a short manual and is not read by AcuStat ###' >> $outputfile
echo "This Script is primary used to accumulate the statistics from the ybar field and write the solution from the last time step." >> $outputfile
echo "This script can also handle additional fiels. In summary, it writes the following fields: " >> $outputfile
echo "        1) solution of the last time step;" >> $outputfile
echo "        2) accumulated ybar field;" >> $outputfile
echo "        3) If requested (optional), accumulated error field;" >> $outputfile
echo "        4) If requested (optional), accumulated phase average fields;" >> $outputfile
echo "        5) If requested (optional), vorticity field (from the last time step);" >> $outputfile
echo "        6) If requested (optional), dwal field (from the last time step)." >> $outputfile
echo "The next lines contain the time steps series for restart files along with the associated number of time steps inside the ybar, error and yphavg fields." >> $outputfile

mv $outputfile $cores-procs_case


