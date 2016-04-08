scp /opt/intel/composer_xe_2015.3.187/compiler/lib/mic/libiomp5.so mic0:/tmp
scp stream.pf mic0:
ssh mic0 "export LD_LIBRARY_PATH=/tmp; ./stream.pf"
