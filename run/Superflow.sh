#!/bin/bash
#SBATCH -p atlas_all
#SBATCH --distribution=cyclic
#SBATCH -N 1 -n 2
#SBATCH --mem-per-cpu=1920
#SBATCH --time=8:00:00

# cd to working directory

echo 'S_ANA_NAME: '${S_ANA_NAME}
echo 'S_MODE: '${S_MODE}
echo 'S_IN_DIRECTORY: '${S_IN_DIRECTORY}  
echo 'S_WORKDIR: '${S_WORKDIR}  
echo 'S_STARTDIR: '${S_STARTDIR}  
echo 'S_SYSTEMATIC: '${S_SYSTEMATIC}  
echo 'S_OUTPUT_DIR: '${S_OUTPUT_DIR}

echo ""
echo ""
echo ""

cd ${S_STARTDIR}
sleep 0.25
cd ${S_STARTDIR}

ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

asetup AtlasProduction 19.1.1.3,gcc48
source RootCore/scripts/setup.sh

echo ""

pwd

echo ""

date
date
date

echo ""

mkdir -p ${S_WORKDIR}
cd ${S_WORKDIR}

pwd

echo ""

${S_ANA_NAME} /i ${S_IN_DIRECTORY} /${S_MODE}

echo ""

sleep 2

cp -fv *.root ${S_OUTPUT_DIR}

echo ""

sleep 1

THIS_DIR=$(basename `pwd`)
cd ..
rm -rf ${THIS_DIR}
ls -l

date
date
date
