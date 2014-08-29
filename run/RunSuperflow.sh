#!/bin/bash
# Submit from a clean environment!
# The scripts set up the environment!
R_ANA_NAME='SuperflowAna'
R_START_DIR='/gdata/atlas/suneetu/Documents/LFV_Higgs2014/S10_Superflow_production/'
R_OUTPUT_DIR='/gdata/atlas/suneetu/Documents/LFV_Higgs2014/output/R3_August_22_patch/'
R_LOG_DIR=${R_OUTPUT_DIR}"logs/"

mkdir -p ${R_LOG_DIR}

R_WORK_BASE='/scratch/suneetu/'

R_GEN='/gdata/atlas/suneetu/Documents/LFV_Higgs2014/generation/'
R_LIST_DIR=${R_GEN}"filelists/"
R_LIST_PREFIX='file_list_'
R_LIST_POST='.txt'

R_SAMPLES_LIST=( "ttbar" "Higgs" "WW" "ZV" "ZPlusJets" )

# SUBMIT DATA
# SUBMIT DATA
# SUBMIT DATA

FILE_NAME=${R_LIST_DIR}'file_list_Data.txt'
FILE_LINES=`cat $FILE_NAME`

for line in $FILE_LINES ; do
	export S_ANA_NAME=${R_ANA_NAME}
	export S_MODE='c' # central for data
	export S_STARTDIR=${R_START_DIR}
	
	sleep 0.13 # to regenerate entropy
	RANNUM=$RANDOM$RANDOM$RANDOM
	
	export S_WORKDIR=${R_WORK_BASE}${RANNUM}
	export S_IN_DIRECTORY=${line%?}'/'
	export S_SYSTEMATIC='NONE'
	export S_OUTPUT_DIR=${R_OUTPUT_DIR}
	
	RED_line=${line%?}
	lFileName=$(basename $RED_line)
	strip_one=${lFileName#user.*.}
	strip_two=${strip_one%.SusyNt*}
	strip_one=${strip_two#group*y.}
	
	echo $strip_one
	
	sbatch -J 'Superflow '${strip_one} -o ${R_LOG_DIR}${lFileName}_slurm-%j.log ${R_GEN}Superflow.sh
done

# SUBMIT NOMINAL + ALL_SYST
# SUBMIT NOMINAL + ALL_SYST
# SUBMIT NOMINAL + ALL_SYST

for file_ in ${R_SAMPLES_LIST[@]}; do
	FILE_NAME=${R_LIST_DIR}${R_LIST_PREFIX}${file_}${R_LIST_POST}
	FILE_LINES=`cat $FILE_NAME`
	
	for line in $FILE_LINES ; do
		export S_MODE='a' # all systematics
		export S_STARTDIR=${R_START_DIR}
		
		sleep 0.13 # to regenerate entropy
		RANNUM=$RANDOM$RANDOM$RANDOM
		
		export S_WORKDIR=${R_WORK_BASE}${RANNUM}
		export S_IN_DIRECTORY=${line%?}'/'
		export S_SYSTEMATIC='NONE'
		export S_OUTPUT_DIR=${R_OUTPUT_DIR}
		
		RED_line=${line%?}
		lFileName=$(basename $RED_line)
		strip_one=${lFileName#user.*.}
		strip_two=${strip_one%.SusyNt*}
		
		echo $strip_two
		
		sbatch -J 'Superflow '${strip_two} -o ${R_LOG_DIR}${lFileName}_slurm-%j.log ${R_GEN}Superflow.sh
	done
done