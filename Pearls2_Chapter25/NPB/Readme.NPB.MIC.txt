Register and download the NASA NAS parallel benchmarks
http://http://www.nas.nasa.gov/publications/npb.html

Build to run natively on Intel Xeon Phi

mkdir NPB_knc
cd NPB_knc
tar -xzf <PathTo NASA tarfile>/NPB3.3.tar.gz 
cd NPB3.3/NPB3.3-OMP
cp config/make.def.template config/make.def
#set F77 to ifort and CC plus UCC variables to icc. Specify -mmic in the options

#build
make ep CLASS=B VERSION=VEC

#copy measuringPowerLab to mic0
scp -r ../../../measuringPowerLab mic0:

#copy ep.B to measuringPowerLab
scp bin/ep.B mic0:measuringPowerLab/ep.B.knc

# login to mic0
ssh mic0
# cd to measuringPowerLab
cd measuringPowerLab
#set the path to find libiomp5.so
sh collectData.MIC.sh










