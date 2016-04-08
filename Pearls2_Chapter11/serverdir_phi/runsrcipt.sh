killall -9 server
./server 5300 ~/intel-img/visualsearch_orig/alldb.txt 0 1 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1 1 &
./server 5301 ~/intel-img/visualsearch_orig/alldb.txt 0 2 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1  1 &
./server 5302 ~/intel-img/visualsearch_orig/alldb.txt 0 4 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1  1 &
./server 5303 ~/intel-img/visualsearch_orig/alldb.txt 0 8 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1  1 &
./server 5304 ~/intel-img/visualsearch_orig/alldb.txt 0 16 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1  1 &
./server 5305 ~/intel-img/visualsearch_orig/alldb.txt 0 32 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1  1 &
./server 5306 ~/intel-img/visualsearch_orig/alldb.txt 0 64 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1  1 &
./server 5307 ~/intel-img/visualsearch_orig/alldb.txt 0 120 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt $1  1 &
