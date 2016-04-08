#This script assumes that the home directory is automounted
PORT=5425
echo "$PWD/pics/s23/6.pgm" > input.txt

serverdir_phi/server $PORT $PWD/serverdir/alldb.txt 0 8 $PWD/data/ $PWD/input.txt 1 1 &
echo "Script is sleeping for a short while"
sleep 10
cd clientdir
./client localhost $PORT ../input.txt 
cd ..
#now kill the server process
killall server



