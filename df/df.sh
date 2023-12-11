#!/bin/bash

#write by hand the snapshot numbers avalaible in your directory 
declare -a snapshots=(1 34 98)
nl=$(cat ../lc/planes_list.txt | wc -l)
declare -a x
declare -a y
for i in $(seq 1 $nl)
do
    #plane number
    x[i]="$(cat ../lc/planes_list.txt | awk -v p="$i" '{if(NR==p) print $1}')"
    #snapshot number
    y[i]="$(cat ../lc/planes_list.txt | awk -v p="$i" '{if(NR==p) print $6}')"
done


for (( i = 0; i < $nl; i++ ))
do
    for (( j=0; j < ${#snapshots[@]}; j++ ))
    do
	if [[ "${snapshots[$j]}" -eq "${y[$i]}" ]];
	   then

	   echo "${x[$i-1]}" > ${x[$i-1]}.d
	   ./executable_df ${y[$i]} ${x[$i-1]}
	   echo "from bash: snap = ${y[$i]}"
	   rm ${x[$i-1]}.d
	fi
    done        
done
