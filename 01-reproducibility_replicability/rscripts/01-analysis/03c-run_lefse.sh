#!/bin/bash 

IFS='_|.'
read -a strarr <<< "$1"

in="${strarr[0]}_${strarr[1]}_${strarr[2]}_${strarr[3]}.txt"
out="${strarr[1]}_${strarr[2]}_${strarr[3]}.out"
res="${strarr[1]}_${strarr[2]}_${strarr[3]}.res"

lefse_format_input.py "${in}" "${out}" -c 1 -u 2
## For Tsementzi threshold for LDA score is 3 and for Chao is default
lefse_run.py "${out}" "${res}" -l 3