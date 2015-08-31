# FILE LICENSE TAG: SAMPLE

if [ $# -gt 0 ] && [ -d $1 ]; then
    bin_dir=$1
else
    bin_dir=../../bin
fi

mpss_dir=`echo /opt/mpss/3.?*/sysroots/k1om-mpss-linux/usr/lib64/`
export SINK_LD_LIBRARY_PATH=$SINK_LD_LIBRARY_PATH:$MKLROOT/lib/mic/:$MKLROOT/../compiler/lib/mic
if [ -d $bin_dir/host ]; then
    pushd $bin_dir/host
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`:/usr/lib64
    if [ -d ../dev ]; then
        pushd ../dev
        export SINK_LD_LIBRARY_PATH=$SINK_LD_LIBRARY_PATH:`pwd`:$mpss_dir
        popd
    else
        export SINK_LD_LIBRARY_PATH=$SINK_LD_LIBRARY_PATH:$mpss_dir
    fi
    popd
else
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
    export SINK_LD_LIBRARY_PATH=$SINK_LD_LIBRARY_PATH:$mpss_dir
fi
# This can have a very large positive impact on card-side MKL performance
# Make this 4G for an 8G card
export MKL_MIC_MAX_MEMORY=8G
# This can have a 10-15% impact on buffer transfer speeds
export MIC_USE_2MB_BUFFERS=64K
