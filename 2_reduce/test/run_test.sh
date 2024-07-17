#!/bin/bash

rm test_fluxes.png test_nres.png testflux.dat
python ../reduce.py -qmmm_list qmmm_list.txt -sufix lmpdt_dIno -nres 307 -cr 306 86 273 -cutoff 10 -ncpus 5 -out test
