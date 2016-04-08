#!/bin/bash
for tid in { 1..32 }
do 
./server 5600 ~/intel-img/visualsearch_orig/alldb.txt 0 1 ~/intel-img/visualsearch/data/ ~/intel-img/visualsearch_orig/input.txt 1 &
../clientdir/client localhost 5600 ~/intel-img/visualsearch_orig/input.txt
done
