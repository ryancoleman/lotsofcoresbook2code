#Set 
MIC=mic0
PATH_IOMP=/opt/intel/composer_xe_2015.3.187/compiler/lib/mic


############################
#Copy library/data/script if they don't exists on Xeon Phi already

#Copy OpenMP library to Xeon Phi
ssh $MIC "test -e libiomp5.so"
if [ $? -eq 1 ]; then
  scp $PATH_IOMP/libiomp5.so $MIC:.
fi

#Copy data to Xeon Phi
ssh $MIC "test -e data"
if [ $? -eq 1 ]; then
  scp -r data $MIC:.
fi

#Copy a script to Xeon Phi
ssh $MIC "test -e script_to_be_run_on_xeon_phi.sh"
if [ $? -eq 1 ]; then
  scp script_to_be_run_on_xeon_phi.sh $MIC:.
fi

############################


#copy executables to Xeon Phi
scp gsfcgce_fortran_MIC gsfcgce_fortran_opt_MIC $MIC:.


ssh $MIC 'cd $HOME; ./script_to_be_run_on_xeon_phi.sh'
