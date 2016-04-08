killall -9 server
./server 5304 ~/intel-img/visualsearch_orig/alldb.txt 0 8 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt 2 2>&1 > server.log &
