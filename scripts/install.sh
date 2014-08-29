#!/bin/bash

#
# This script is obsolete! See installMinimalSUSYTools.sh instead, 
# and the general instructions in the README of this package
#

# It is better to run kinit before this script

# Run this from your work dir, not within this package
asetup AtlasProduction,17.2.9.1,setup,here

# For now, we will need SUSYTools and Mt2 and DGTriggerReweighting
# Do we really still need all of these? Commenting them out for now
#svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/SUSYTools-00-02-08 SUSYTools
#svn co svn+ssh://svn.cern.ch/reps/atlasphys/Physics/SUSY/Analyses/WeakProduction/Mt2/tags/Mt2-00-00-01 Mt2
#svn co svn+ssh://svn.cern.ch/reps/atlasphys/Physics/SUSY/Analyses/WeakProduction/DGTriggerReweight/tags/DGTriggerReweight-00-00-25 DGTriggerReweight
#python SUSYTools/python/install.py

# For BTag calibration, need to check out CalibrationDataInterface
calibTag="svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/JetTagging/JetTagPerformanceCalibration/CalibrationDataInterface/tags/CalibrationDataInterface-00-03-04"
svn co $calibTag/cmt $calibTag/CalibrationDataInterface $calibTag/Root $calibTag/src CalibrationDataInterface

# Configure RootCore
cd RootCore
./configure
source scripts/setup.sh
cd ..

# Compile everything
$ROOTCOREDIR/scripts/build.sh
