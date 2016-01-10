#!/bin/bash

for f in *.rst
do 
	CMD="rst2html.py -v $f ${f%.rst}.html --stylesheet-path=css/voidspace.css --embed-stylesheet"
	#CMD="rst2html.py -v $f ${f%.rst}.html"
	echo $CMD
	$CMD
done


