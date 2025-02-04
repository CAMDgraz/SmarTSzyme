# SmarTSzyme
![alt text](./cover_temp2.png)
**SmarTSzyme** is a command-line interface (CLI) designed for selecting key residues to manipulate enzyme activity in the (re)design of enzymes. By simulating the reaction mechanism of the desired enzyme activity at the active site of the target enzyme using quantum mechanics/molecular mechanics (QM/MM) methods, **SmarTSzyme** identifies residues that (de)stabilize the enzyme-transition state (ETS) complex relative to the enzyme-substrate (ES) complex. By guiding the mutational landscape, **SmarTSzyme** facilitates enzyme activity manipulation with atomic-level attention to the catalytic mechanism, ultimately saving time and reducing costs in protein engineering campaigns.

**SmarTSzyme** consist of three modules:

1. SIMULA, for the preparation and simulation of the catalyzed reaction using QM/MM MD simulations (

## Installation
### SmarTSzyme
To run **SmarTSzyme**,
```bash
git clone https://github.com/CAMDGraz/SmarTSzyme.git
```
```bash
cd /path/to/SmarTSzyme
conda env create -f environment.yml
conda activate smartszyme_env
```

## Preparing the steered Molecular Dynamics files
For preparing the steered Molecular Dynamics files we provide a bash script (./2_reduce/prep_smd.sh) that can be use as follows:
```bash
./2_reduce/prep_smd.sh <joblist> <topext> <trajext> <suffix>
# <joblist> => file with the path to the qmmm folders (an example can be found in ./example/)
# <topext>  => extension of the topology file (by default the script is gonna load the file matching *.<topext>)
# <trajext> => extension of the trajectory file (by default the script is gonna load the file matching *.qmmm.<trajext>)
# <suffix>  => suffix for the output files (top_<suffix>.parm7, traj_<suffix>.nc and smd_<suffix>.txt)
```
The script will remove the waters, Na+ and Cl- and create the input files for reduce.py with the following structure:
```bash
top_<suffix>.parm7
traj_<suffix>.nc
smd_<suffix>.txt
```
Modify the script according to your systems deleting everything but the protein, substrate and cofactor if any. Please make sure that a single txt file (the one resulting from the sMD) is present in the paths. 

## Basic Usage
You can display the help of **SmarTSzyme** in the command-line as follows:
```bash
python 2_reduce/reduce.py -h

********************************************************************************
* SmarTSzyme-reduce:                                                           *
*      Selection of important residues for enzyme engineering                  *
********************************************************************************

usage: Reduce [-h] -qmmm_list QMMM_LIST -sufix SUFIX -nres
              NRESIDUES -cr CATALYTIC_RESIDUES
              [CATALYTIC_RESIDUES ...] -cutoff CUTOFF
              [-ncpus NCPUS] -out OUTPUT

Reduce: Identification of key residues in order to reduce the the
mutational landscape

options:
  -h, --help            show this help message and exit

Input options:
  -qmmm_list QMMM_LIST  List of QMMM jobs to analyze.
  -sufix SUFIX          Sufix for the top_, traj_ and smd_ files
  -nres NRESIDUES       Number of residues
  -cr CATALYTIC_RESIDUES [CATALYTIC_RESIDUES ...]
                        Catalytic residues
  -cutoff CUTOFF        Maximum distance between residues to be
                        consider in the pairwise interactions (in
                        A)
  -ncpus NCPUS          Number of CPUs to use [default: 1]

Output options:
  -out OUTPUT           prefix for the outputs
```

In the the example folder, you will find a topology, trajectory and steered MD output file to run a test. For reproducing the results type the following command in the console:

```bash
cd example
python ../2_reduce/reduce.py -qmmm_list job_list.txt -sufix mhet -nres 562 -cr 183 450 486 562 -cutoff 10 -ncpus 1 -out out_reduce
```
## License
**SmarTSzyme** is licensed under GNU General Public License v3.0.

## Citation
The corresponding publication is under preparation

## Contact
**Computed-Aided Molecular Design**

Division of Medicinal Chemistry\
Medical University of Graz\
Neue Stiftingstalstraße 6/III\
A-8010 Graz, Austria

Head of the Group: Ass.-Prof. Dr. Pedro A. Sánchez Murcia
 
In case of questions and/or suggestions you can contact us at: pedro.murcia@medunigraz.at and daniel.platero-rochart@medunigraz.at
