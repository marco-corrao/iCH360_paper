This folder contains the pipeline used to construct the reduced model variant *i*CH360, as well as counting and enumerating its elementary flux modes

#### `./iCH360red_construction`
Contains files and code required to prune the relevant reactions from *i*CH360 and export the reduced model variant.
#### `./efm_enumeration`
contains the relevant scripts to count, enumerate and filter EFMs. 

### A note on EFM enumeration/counting
For the purpose of this paper, we counted *all* EFMs of *i*CH360red, but only explicitely enumerated a subset of them, which we refer to as *filtered*, that exclude a number of physiologically unrealistic EFMs from the total, (see the paper or more details on this).

the **EFM counting** script (`./efm_enumeration/count_efms.m`) counts *all* EFMs, but doesn't explicitely store them in memory. At the time of writing, we were unable to perform use EFMtools with the count-only option from the python API, so that the counting script has to be run in Matlab. To run the script, 1) download the EFMTool Matlab toolbox (https://csb.ethz.ch/tools/software/efmtool.html); 2) Run the script `./efm_enumeration/create_matlab_structs_for_efm_counting.ipynb` to construct and export the condition-specific Matlab models to be with EFMtool; 3) move `count_efms.m` and the folder `matlab_models` in the root of the EFMtool package directory; 4) run `count_efms.m`: this will reproduce the count file `efm_counts/efm_counts_unfiltered.csv` which we included in this directory.

On the other hand, the **EFM enumeration** script (`./efm_enumeration/enumerate_efms.ipynb`) runs directly on python and requiresthe `efmtools` python package (https://pypi.org/project/efmtool/).  This will explicitely enumerate *filtered* EFMs and save them to file, which can be memory intensive for some growth conditions ( ~$10^6$ filtered modes for growth on glucose). We could, however, run the enumeration script on a standard laptop with 16 GB of RAM. The output of this script is stored in the `efm/` subfolder and is not included in this repository due to the large size of the files.
