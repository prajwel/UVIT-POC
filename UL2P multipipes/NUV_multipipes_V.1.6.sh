#!/bin/csh

which UVIT_DriverModule

# Random name generation 
set random_name = "`openssl rand -hex 4`"

# To set local environment variables and what not.
setenv PFILES "$HOME/multipipe_$random_name/paramfiles"
setenv GLOG_log_dir "$HOME/multipipe_$random_name/log"
mkdir -p $PFILES $GLOG_log_dir
cp -R $AS1/uvit/paramfiles/* $PFILES

#To copy VIS driver to paramfiles.
cp $HOME/switch_drivers/NUV_UVIT_DriverModule.par $PFILES/UVIT_DriverModule.par

# Run the driver module
UVIT_DriverModule n LEVL1AS1UVT*tar_V*

# To copy the driver module param file
cp $PFILES/UVIT_DriverModule.par .

# To remove the temporary log and paramfiles directory.
rm -rf $HOME/multipipe_$random_name

# To remove unwanted files
rm -rf AUX4*
rm -rf *_level*
rm -f *_EDITED_L2
rm -f *PCONLY*
rm -f Differences.txt
rm -f Driver_reference.txt
rm -f Driver_Total_NUV_SCIENCEDATAFILE.txt
rm -f xytheta.txt
rm -f ZeroCentroid.txt
rm -f core.*
find . -name 'DataIngest_*' -type d -exec rm -rf {} +
cd output
find . -name '*uvtFrameIntegration*' -type d -exec rm -rf {} +

# To remove the L2tars. (IMPORTANT!)
rm -rf /data/pipeline/prajwel/L2tar/*



