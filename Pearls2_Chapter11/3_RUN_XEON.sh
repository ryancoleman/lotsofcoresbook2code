PORT=5327
echo "$PWD/pics/s23/6.pgm" > input.txt

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/../flann-1.8.4-src/build/lib/

serverdir/server $PORT $PWD/alldb.txt 0 8 $PWD/data/ $PWD/input.txt 1 1 &

echo "Script is sleeping for a short while"
sleep 10
cd clientdir
./client localhost $PORT ../input.txt 
cd ..
#now kill the server process
killall server



