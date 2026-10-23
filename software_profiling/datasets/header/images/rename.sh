#!/bin/bash

declare -i x=25
for i in {1..22}; do
	mv ./dataset_2/image_data_$i.h ./dataset_2/image_data_$x.h
	x+=1


done
