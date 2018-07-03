#! /bin/bash

while getopts ":e:h" opt; do
	case $opt in
		e) END="$OPTARG" ;;
		h) echo "-e is the ending number of the sequence" ;;
		\?) echo "invalid option -$OPTARG" >&2 ;;
	esac
done

#for i in $(seq 1 $END); do echo $i; done

for i in $(seq 1 $END);
do
	dd if=/dev/zero of=/dev/null &
	#stress --cpu 1 --timeout 60 &

done

# use: killall dd
# to kill all instances of dd
