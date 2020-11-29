#!/bin/bash

set -o xtrace

LIST_PDBS=(1yv3 1ywr 1mzc)

for ipdb in ${LIST_PDBS[@]}; do
  echo " "
  #make PDB=$ipdb eval
  make PDB=$ipdb OMP=YES eval
done
