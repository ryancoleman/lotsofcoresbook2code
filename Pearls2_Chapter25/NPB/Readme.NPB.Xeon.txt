Register and download the NASA NAS parallel benchmarks
http://http://www.nas.nasa.gov/publications/npb.html

Build to run on Xeon

mkdir NPB_xeon
cd NPB_xeon
tar -xzf <PathTo NASA tarfile>/NPB3.3.tar.gz 
cd NPB3.3/NPB3.3-OMP
cp config/make.def.template config/make.def
#set F77 to ifort and CC plus UCC variables to icc. 

#build
make ep CLASS=B VERSION=VEC

#copy ep.B to measuringPowerLab
cp bin/ep.B ../../../measuringPowerLab/ep.B

cd ../../measuringPowerLab
sh collectData.sh










