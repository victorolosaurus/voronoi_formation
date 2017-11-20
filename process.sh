#!/bin/bash

#convert the csv into whitespace separated file (also just print tag_id and position, theres more data in the set)
#note: one might want to grep out tagid 12 here, looks like a bench player thats often detected as being on field
awk -F, 'BEGIN{last="";num=0;}(last!=$1){print ""; last=$1;num=num+1}(num%20==0){print $2,$3,$4}' 2013-11-03_tromso_stromsgodset_first.csv > tromso.txt

#split into separate files (called whatever-????
awk -v RS= '{print > ("whatever-" NR ".txt")}' tromso.txt

#process all whatever files
for i in whatever-*.txt; do j=${i%%.txt}; python voronoi_formation2.py $i $j; done
