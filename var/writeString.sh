#!/bin/bash

str=$1

while test -n "$str"
do
  xdotool key ${str:0:1}
  str=${str:1}
done
