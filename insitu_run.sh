#!/bin/bash

# Run batch processes for insitu_stats library

global_size=( 4 8 16 32 64 128 256 512)
local_size=( 4 8 16 32 64 128 256 512)
mode=( 64 32 16 8 4 1 )

#run_command="mpiexec"
run_command="runjob"
LOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"
stand_input="--verbose=INFO :"
run_file="./3D_Grid"
out_file="run_results.txt"

make
rm $out_file
touch $out_file

let "count= 0"
# Run file and write to stdout and to file
for i in "${local_size[@]}"
do
  l_param=$(echo "-l $i/x$i/x$i" | tr -d "/")
  
  for j in "${global_size[@]}"
  do 
    g_param=$(echo "-g $j/x$j/x$j" | tr -d "/")
    let "div= $j / $i"
    if [ $div -gt 0 ]; then
      let "blocks= $div * $div * $div"
      #num_proc="-n $num"
      num_proc="--np $blocks"
      #let "k=2"
      #per_node="-p $k"
      
      for k in "${mode[@]}"
      do
        let "p= $blocks % $k"
      	let "nodes=$blocks/$k"
        #echo "np:$num k:$k p:$p"
        if [ $p == 0 ] && [ $nodes -lt 1025 ]; then 
          per_node="-p $k"
          #break
        #fi
      #done
      let "count+=1"
      #echo "******* RUN:$count: **************" >> $out_file
      #echo "$j $i $k $nodes" 
      #echo "******* RUN:$count: **************" | sudo tee -a $out_file
      #echo "Nodes:$nodes" >> $out_file
      #echo "Processors:$num:" >> $out_file
      #echo "***** Processors:$num: *********" | sudo tee -a $out_file
      #echo "$run_command $num_proc $run_file $g_param $l_param | tee /A $out_file"

      #$run_command $num_proc $run_file $g_param $l_param >> $out_file
      #$run_command $num_proc $run_file $g_param $l_param | sudo tee -a $out_file
      #of the form: runjob --np # -p # --block $COBALT_PARTNAME --verbose=INFO : run1.exe arg1
      #from: http://www.alcf.anl.gov/user-guides/running-jobs#submitting-a-script-job
      echo "$run_command $LOCARGS $num_proc $per_node $stand_input $run_file $g_param $l_param $nodes >> $out_file"
        fi
      done
    fi
  done
done
