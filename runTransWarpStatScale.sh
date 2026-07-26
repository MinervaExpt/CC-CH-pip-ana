
#!/bin/bash

#Usage: runTransWarp.sh runEventLoopMC.root warped.root

#declare -a VARIABLE=("mixtpi" "enu" "pmu" "ptmu" "pzmu" "q2" "wexp")
declare -a VARIABLE=("mixthetapi_deg")
#declare -a VARIABLE=("thetamu_deg")
#declare -a warps=("WARP2" "WARP3" "WARP4" "WARP1")
declare -a warps=("NOMINAL")
#declare -a VARIABLE=("mixtpi")
declare -a STAT_SCALE=("1" "2" "3" "4" "4.689984" "5" "6" "7" "8" "9" "10")
#declare -a STAT_SCALE=("10.2" "10.4" "10.6" "10.8" "11" "11.5" "12" "12.5" "13" "13.5" "14" "15")
#declare -a STAT_SCALE=("18.5" "18.5" "19" "19.5" "20")
#declare -a STAT_SCALE=("1.1" "1.2" "1.3" "1.4" "1.5" "1.6" "1.7" "1.8" "1.9")
#declare -a STAT_SCALE=("4.2" "4.4" "4.6" "4.8" "5.2" "5.4" "5.6" "5.8")
#declare -a warps=("WARP1")
Ximaxaxis=(1500 10000 1000 100 650 75 4000 100 100 50 550 275 800 10000 10000 200 100 2200 75 4000 100 50 50 850 150 150 1100 10000 1200 100 450 75 4500 100 75 50 200 50 550 1000 10000 1400 100 450 75 5000 100 250 50 500 100 5000)
#Ximaxaxis=(300)
DATE="20241231"
Plist="ALL"
#MIGRATION_FILE=$1
#TRUE_HIST=effnum_${VARIABLE}
#WARPED_FILE=$2
#RECO_HIST=selection_mc_${VARIABLE}

#OUTFILE_NAME="/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/WarpingStudies/TpiWeightFixed/"
OUTFILE_NAME="/minerva/data/users/granados/WarpingStudies/Warping/1DWarping/mixedp4tpiestDec2023AllPlist/StatVariationsmixtpi/"
#OUTFILE_NAME=$(basename $2)
counter=0
for st_scale in "${STAT_SCALE[@]}"; do
  for TAG in "${warps[@]}"; do
    for v in "${VARIABLE[@]}"; do
      MIGRATION_FILE="MCXSecInputs_${DATE}_${Plist}_mix_thetapi_CCPiWeighNoSys_p4_NOMINAL.root"
      WARPED_FILE="MCXSecInputs_${DATE}_${Plist}_mix_thetapi_CCPiWeighNoSys_p4_${TAG}.root"
  
      TransWarpExtraction --output_file ${OUTFILE_NAME}Warping_${Plist}_mix_thetapi_CCPiWeightNoSys_stat_scale${st_scale}_${DATE}_${TAG}_${v}.root --data effnum_${v} --data_file $WARPED_FILE --data_truth effnum_${v}_true --data_truth_file $WARPED_FILE --migration migration_${v} --migration_file $MIGRATION_FILE --reco effnum_${v} --reco_file $MIGRATION_FILE --truth effnum_${v}_true --truth_file $MIGRATION_FILE --num_uni 500 --step_chi2 0.5 --num_iter 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,50,100,200 -C 0.5 --stat_scale ${st_scale} # --max_chi2 ${Ximaxaxis[${counter}]}
      #echo "Variable ${v} Warp ${TAG} Xi Y Axis ${Ximaxaxis[${counter}]}"
  #    cd ../MAT/macros/
  #    python PrintWarpingStudy.py -i ${OUTFILE_NAME}Warping_${TAG}_${v}.root -o ${OUTFILE_NAME}${v}_${TAG} -L
  #    cd -
  #let counter=counter+1
    done # vars
  done #warps
done #Statistical scale variation

