#!/bin/bash

for d in *; do
    echo "$d"
    ../bin/tpl-mcce2free.py "$d" > tmp.ftpl
    ../bin/vdw-complete.py tmp.ftpl > ../test/"$d"
done

