#!/bin/sh

# ---- Set case name ----
case_name='20191205'
# -----------------------


# ---- read file list of mocm_vtk ----
file=()
time2=`ls -v mocm_vtk/${case_name}/`
i=1
for i in ${time2[@]}
do
    file=(${file[@]} $i)
done
# -----------------------


# ---- make directory of mocm_vtk ----
nend=8
for j in `seq $nend`; do # you should change numbers which you wanna use numbers of cpu
    mkdir mocm_vtk/${case_name}_$j
done
# -----------------------


# ---- move file list to another directory of mocm_vtk ----
i=0
k=0
n=1
for stage in ${file[@]}; do
    if [ $k -le $(( ${#file[@]} / $(($nend)))) ]; then
        sleep 0.5
        mv mocm_vtk/${case_name}/${file[i]} mocm_vtk/${case_name}_$n/${file[i]}
        if [ $k -eq $(( ${#file[@]} / $(($nend)))) ]; then
            n=$(( $n + 1 ))
        fi
        if [ $k -eq $(( ${#file[@]} / $(($nend)))) ]; then
            k=1
        else
            k=$(( $k + 1 ))
        fi
    else
        k=1
    fi
    i=$(( $i + 1 ))
done
# -----------------------


# ---- do 01_Trim_median filter ----
for j in `seq $((nend+1))`; do
    echo ${case_name}_$(($j-1)) >Target_file.txt &&
    sleep 5
    echo ${case_name}_$(($j-1))
    cd 01_Trim_Median/
    make && ./median.o -db & cd .. && continue
done
# -----------------------

wait

# ---- move file list to original directory of 00_input_vtk ----
mkdir 00_input_vtk/${case_name}
mkdir 00_input_vtk/${case_name}/vtk_schalar_median
mkdir 00_input_vtk/${case_name}/vtk_vector_median

i=0
k=0
n=1
for stage in ${file[@]}; do
    if [ $k -le $(( ${#file[@]} / $(($nend)))) ]; then
        mv 00_input_vtk/${case_name}_$n/vtk_schalar_median/${file[i]}_trm.vtk 00_input_vtk/${case_name}/vtk_schalar_median/${file[i]}_trm.vtk
        mv 00_input_vtk/${case_name}_$n/vtk_vector_median/${file[i]}_trm.vtk 00_input_vtk/${case_name}/vtk_vector_median/${file[i]}_trm.vtk
        if [ $k -eq $(( ${#file[@]} / $(($nend)))) ]; then
            n=$(( $n + 1 ))
        fi
        if [ $k -eq $(( ${#file[@]} / $(($nend)))) ]; then
            k=1
        else
            k=$(( $k + 1 ))
        fi
    else
        k=1
    fi
    i=$(( $i + 1 ))
done
# -----------------------

wait

# ---- remove another directory of mocm_vtk ----
for j in `seq $nend`; do
    rm -r 00_input_vtk/${case_name}_$j
done
# -----------------------


# ---- move file list to original directory of mocm_vtk ----
i=0
k=0
n=1
for stage in ${file[@]}; do
    if [ $k -le $(( ${#file[@]} / $(($nend)))) ]; then
        mv mocm_vtk/${case_name}_$n/${file[i]} mocm_vtk/${case_name}/${file[i]}
        if [ $k -eq $(( ${#file[@]} / $(($nend)))) ]; then
            n=$(( $n + 1 ))
        fi
        if [ $k -eq $(( ${#file[@]} / $(($nend)))) ]; then
            k=1
        else
            k=$(( $k + 1 ))
        fi
    else
        k=1
    fi
    i=$(( $i + 1 ))
done
# -----------------------


# ---- remove another directory of mocm_vtk ----
for j in `seq $nend`; do
    rm -r mocm_vtk/${case_name}_$j
done
# -----------------------
