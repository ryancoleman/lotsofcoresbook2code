#export MIC_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH:/home1/02857/prash/intel-img/flann-phi/build/lib/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home1/02857/prash/intel-img/flann-1.8.4-cpu/build/lib/
export MIC_ENV_PREFIX=PHI
export PHI_KMP_AFFINITY=balanced
export PHI_KMP_PLACE_THREADS=60c,2t
export PHI_OMP_NUM_THREADS=120
export PROJ_DIR=/home1/02857/prash/intel-img/visualsearch_orig
echo Done
./server 5304 $(PROJ_DIR)/alldb.txt 0 8 $(PROJ_DIR)/data $(PROJ_DIR)/input.txt 1 1 &
sleep 3
../clientdir/client localhost 5304 $(PROJ_DIR)/input.txt
