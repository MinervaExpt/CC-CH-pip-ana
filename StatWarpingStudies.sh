

#!/bin/bash

#Usage: runTransWarp.sh runEventLoopMC.root warped.root

#declare -a VARIABLE=("mixtpi" "wexp" "pmu" "ptmu" "pzmu" "q2" "thetamu_deg")
declare -a VARIABLE=("mixthetapi_deg")
#eclare -a warps=("WARP1" "WARP2" "WARP3" "WARP4")
#declare -a warps=("NOMINAL")
#declare -a VARIABLE=("thetamu_deg")
#declare -a STAT_SCALE=("1" "2" "3" "4" "4.689984" "5" "6" "7" "8" "9" "10")
#declare -a STAT_SCALE=("4.689984" "5" "6" "7" "8" "9" "10")
#declare -a STAT_SCALE=("10.2" "10.4" "10.6" "10.8" "11" "11.5" "12" "12.5" "13" "13.5" "14" "15")
#declare -a STAT_SCALE=("10.2")
#declare -a STAT_SCALE=("1.1" "1.2" "1.3" "1.4" "1.5" "1.6" "1.7" "1.8" "1.9")
#declare -a STAT_SCALE=("2.1" "2.2" "2.3" "2.4" "2.5" "2.6" "2.7" "2.8" "2.9")
#declare -a STAT_SCALE=("4.2" "4.4" "4.6" "4.8" "5.2" "5.4" "5.6" "5.8")
declare -a STAT_SCALE=("4.689984")
declare -a warps=("NOMINAL")
#Ximaxaxis=(300)
DATE="20241231"
Plist="ALL"
#unfac="100000000"
unfac="6.2"
#TRUE_HIST=effnum_${VARIABLE}
#RECO_HIST=selection_mc_${VARIABLE}

#OUTFILE_NAME="/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/WarpingStudies/TpiWeightFixed/"
OUTFILE_NAME="/minerva/data/users/granados/WarpingStudies/Warping/1DWarping/mixedp4tpiestDec2023AllPlist/StatVariationsmixtpi/"
#OUTFILE_NAME=$(basename $2)
counter=0
for st_scale in "${STAT_SCALE[@]}"; do
  for TAG in "${warps[@]}"; do
    for v in "${VARIABLE[@]}"; do
      ls ${OUTFILE_NAME}Warping_${Plist}_mix_thetapi_CCPiWeightNoSys_stat_scale${st_scale}_${DATE}_${TAG}_${v}.root
      python ../MAT/macros/ProcessMCSampleSizeScan.py --n_uni 500 --iters 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,50,100,200 --input ${OUTFILE_NAME}Warping_${Plist}_mix_thetapi_CCPiWeightNoSys_stat_scale${st_scale}_${DATE}_${TAG}_${v}.root --uncfactor ${unfac} --f_option_used_transwarp ${st_scale}

      #echo "Variable ${v} Warp ${TAG} Xi Y Axis ${Ximaxaxis[${counter}]}"
  #    cd ../MAT/macros/
  #    python PrintWarpingStudy.py -i ${OUTFILE_NAME}Warping_${TAG}_${v}.root -o ${OUTFILE_NAME}${v}_${TAG} -L
  #    cd -
  #let counter=counter+1
    done # vars
  done #warps
done #Statistical scale variation

