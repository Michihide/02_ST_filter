#!/bin/sh

# ---- Set case name ----
  case_name='20191205'
# -----------------------


# ---- read 00_input_file_list ----
time_out=()
while read line || [ -n "${line}" ]
do
    time1=`echo ${line} | cut -c 1-4` # if file name is n digit, you should change number 4 to n
    time_out=(${time_out[@]} $time1)    
done < 00_input_file_list/${case_name}.csv
# -----------------------


# ---- read directory of mocm_vtk ----
cd ..
cd mocm_vtk/${case_name}/ 
time_in=()
time2=`ls -v | sort -n`
for i in ${time2[@]}; do
    time_in=(${time_in[@]} $i)
done
# -----------------------


# # -----------------------
# i=0
# for stage in ${time_out[@]}; do
#     echo $i ${time_in[i]} ${time_out[i]}
#     i=$(( $i + 1 ))
# done
# # -----------------------


# ---- change directory name of mocm_vtk ----
i=0
for stage in ${time_out[@]}; do
    mv -i ${time_in[i]} ${time_out[i]}
    sleep 1
    echo $i ${time_in[i]} ${time_out[i]}
    i=$(( $i + 1 ))
done
# -----------------------
