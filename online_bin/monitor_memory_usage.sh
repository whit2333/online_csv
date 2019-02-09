#!/bin/bash

while true; do
#ps -C hcana -o pid=,%mem=,vsz= >> /tmp/mem.log
ps -C hcana -o pid=,%mem=,vsz= | sed '1{N;s/\n/ /;}'  | sed '1{N;s/\n/ /;}' >> /tmp/mem.log
gnuplot show_mem.plt
sleep 1
done 

#ps --pid $1 -o pid=,%mem=,vsz= >> /tmp/mem.log
