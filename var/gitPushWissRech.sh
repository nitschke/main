#!/bin/bash

MYGITPATH=~/git/WissRech/main

cd $MYGITPATH
git add -vA
git commit -m $1
git push -v
cd -
