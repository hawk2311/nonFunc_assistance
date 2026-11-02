#!/bin/bash


convert -resize 512X512 $1 "512x512_${1}"

convert -resize 720X720 $1 "720x720_${1}"

#convert -resize 1280X720 $1 "1280x720_${1}"

#convert -resize 1080X1080 $1 "1080x1080_${1}"

#convert -resize 720X512 $1 "720x512_${1}"
