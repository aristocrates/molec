#!/bin/bash

echo "cd $1"
cd $1
echo "mkdir workspace"
mkdir workspace
echo "cp *.png workspace"
cp *.png workspace
echo "cd workspace"
cd workspace
echo "mogrify -resize 480x360 *.png"
mogrify -resize 480x360 *.png
echo "convert -delay 1 -loop 0 *.png $1.gif"
convert -delay 1 -loop 0 *.png $1.gif
echo "cp $1.gif ../../animated_gif/$1.gif"
cp $1.gif ../../animated_gif/$1.gif
echo "cd ../.."
cd ../..
