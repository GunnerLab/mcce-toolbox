#!/bin/bash

for d in param04/*; do
    echo $d
    fname=$(basename $d)
    base=${fname%.tpl}
    ../bin/tpl-mcce2free.py "$d" > a.ftpl
    ../bin/rxn-complete.py a.ftpl 4 > b.ftpl
    ../bin/vdw-complete.py b.ftpl > test/"$base.ftpl"
done

