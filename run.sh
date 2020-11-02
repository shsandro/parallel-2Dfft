#!/bin/bash

helpFunction() {
   echo "Usage: $0 -n -p"
   echo -e "\t-n <number> number of times that the programs will run"
   echo -e "\t-p <n° of processes> sets the number of processes"
}

if [ $# -eq 0 ] 
    then
    helpFunction
    exit 1
fi

while getopts "n:p:h" opt
do
   case "$opt" in
      n ) n_repeticoes="$OPTARG" ;;
      p ) n_processes="$OPTARG" ;;
      h ) helpFunction; exit 1 ;; # Print helpFunction in case parameter is non-existent
      * ) exit 1
   esac
done

seq_executable="seq_2dfft.out"
par_executable="par_2dfft.out"
sizes=(2500 5000 10000 15000)

seq_time=()
par_time=()
par_time_comms=()

echo "Execução para ${n_processes} processos:"
echo
echo

for size in ${sizes[@]}; do
    echo "Size ${size}"
    ./${seq_executable} -n ${size} > /dev/null
    for i in $(eval echo "{1..$n_repeticoes}"); do
        seq_time+=($(./${seq_executable} -n ${size} | grep "Time elapsed" | grep -Eo "?([0-9]*[.])?[0-9]+"))
        par_run=$(sudo perf stat -o output/perf_${size}_${n_processes}_${i} -B -e cache-references,cache-misses,cycles,instructions,branches,branch-misses,faults mpirun --host ${HOSTNAME}:${n_processes} ./${par_executable} -n ${size})
        par_time+=($(echo -e ${par_run} | grep "Time elapsed" | grep -Eo "?([0-9]*[.])?[0-9]+" | head -n1))
        par_time_comms+=($(echo -e ${par_run} | grep "Time with comms" | grep -Eo "?([0-9]*[.])?[0-9]+" | tail -n1))
    done

    python3 speedup.py $size $n_processes ${n_repeticoes} ${seq_time[@]} ${par_time[@]} > output/speedup${n_processes}_size_${size}
    python3 speedup.py $size $n_processes ${n_repeticoes} ${seq_time[@]} ${par_time_comms[@]} > output/speedup_comm${n_processes}_size_${size}

    seq_time=()
    par_time=()
    par_time_comms=()
done 