#!/bin/bash

for ((COUNT = 1; COUNT <= 19; COUNT ++)); do 
	printf "=================================\n"
	printf "=== test$COUNT\n"
	printf "=================================\n"
	diff --ignore-all-space test$COUNT.out benchmark/
done
