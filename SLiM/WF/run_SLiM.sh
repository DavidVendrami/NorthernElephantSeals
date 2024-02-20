#!/bin/bash

for i in {1..100}
do
slim INPUT_rep$i.txt > rep_$i/F_$i.txt
done
