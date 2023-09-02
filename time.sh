#! /bin/bash

ipara=$(awk '/iparallel/{print $1}' AFSSH.inp)

for((i=1;i<=$ipara;i++));do
  echo -n $i
  tail $i/output|grep 'Total time'
done

