# CC-CH-pip-ana
Analysis code for single charged pion production by muon neutrinos in the MINERvA CH detector and in the "medium energy" NuMI beam era.

The main branch now uses MasterAnaDev (MAD) tuples as the input, general-reconstruction event record. And it also makes nearly complete use of the MINERvA Analysis Toolkit (MAT) to perform the full cross section calculation with systematics.

## Reproducing Cross Section Results:
1. **Create MC histograms** -- Loop over MC tuple files to create all needed truth and reco-MC histograms (~1-3 hours).
2. **Create data histograms** -- Loop over data tuple files to create all needed data histograms (~30 minutes).
3. **Combine MC and data histograms into a cross section** --  background, unfolding, efficiency and normalization (~ seconds).
4. **Plot results** -- (~seconds).

Step 1 involves running the script `xsec/makeCrossSectionMCInputs.C`. This script could be run over all playlists on a single machine over the course of days. Instead, `ProcessCCPiMacro.py` parallelizes the task by submitting to the Fermilab grid `makeCrossSectionMCInputs` jobs over single input tuples. This way, all jobs start and finish in a few hours, and they only need to be merged after using `MAT/macros/madd.cxx` (which is like ROOTs `hadd` except it merges root files containing `MnvHND` histograms.

Steps 2 and 3 are combined into a single script: `crossSectionDataFromFile.C`. This script first performs step 2, which takes less than an hour on a single machine, and then does step 3 immediate after. Input is the MC input file resulting from step 1.

Step 4 takes the input file from step 3 and just takes a few seconds to run on a single machine. All cross section steps -- including event selection, sideband tune, background subtraction, migration matrix, unfolding, efficiency, and cross sections -- are plotted along with systematics in many variables.

## Reproducing Cross Section Diagnostic Results:
Other plots, not strictly apart of the cross section calculation pipeline, are made by running other scripts, mostly in the studies `studies/` folder. They include cut-by-cut event selection truth breakdowns (`runCutVars`), efficiency-purity counts table (`runEffPurTable`), sideband breakdowns (`runSidebands`), and background truth breakdowns (`runBackgrounds`). There is also `xsec/binningStudy`. The closure test is a three-script process (1. modified `xsec/makeCrossSectionMCInputs.C`, 2. `runXSecLooper.cpp` and 3. `xsec/GXSEClosure.C`) that will be described in a separate section. The warping study to test unfolding just involves running the cross section pipeline with modifications to `xsec/makeCrossSectionMCInputs.C`.

## Setup

* Install [MAT-MINERvA](https://github.com/MinervaExpt/MAT-MINERvA) and [UnfoldUtils](https://github.com/MinervaExpt/UnfoldUtils).
  * [GENIEXSecExtract](https://github.com/MinervaExpt/GENIEXSecExtract) is only needed if you will be performing closure tests.
  * [MAT](https://github.com/MinervaExpt/MAT) is automatically installed by `MAT-MINERvA`
```
> tree your_working_dir/ -L 1
├── CC-CH-pip-ana
├── GENIEXSecExtract
├── MAT
├── MAT-MINERvA
├── opt
└── UnfoldUtils
```
* Set your `$TOPDIR` environment variable `export TOPDIR=/path/to/your_working_dir/`
* Run other environmental setup: `source permissions.sh; source setSL7.sh`

## Notes
* lots more details to be added here.
* legacy instructions are [here](https://docs.google.com/document/d/1fC22eg0K82p5Jc_UwUKmOgB-SSq5zDVDBe0p1iXLvxQ/edit?usp=sharing). They still give accurate direction for running scripts, but they're setup section is out of date.
