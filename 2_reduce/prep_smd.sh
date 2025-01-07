#!/bin/bash

joblist=$1
topext=$2
trajext=$3
suffix=$4

while IFS="" read -r p || [ -n "$p" ]
do
  cd $p/
  rm reduce.in
  echo "parm *.${topext}
  trajin *.qmmm.${trajext}
  autoimage
  center
  strip :WAT,Na+,Cl- parmout top_${suffix}.top
  trajout traj_${suffix}.nc" >> reduce.in
  cpptraj < reduce.in
  cp *.txt smd_${suffix}.txt
  cd ../
done < ${joblist}
