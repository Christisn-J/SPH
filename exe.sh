#!/bin/bash 

make all
./bin/exe
./tools/visualisation/2D.py -d ./out/
ffmpeg -framerate 30 -pattern_type glob -i "out/*.png" -c:v libx264 -pix_fmt yuv420p ./out/movie.mp4