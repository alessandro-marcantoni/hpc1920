#!/bin/bash

for f in *.in
do
	inputName=$f
	hullName=outputs/${f:0:-3}.hull
	
	echo -e $f > outputs/${f:0:-3}.txt
	./convex-hull < $f >> outputs/${f:0:-3}.txt
done
