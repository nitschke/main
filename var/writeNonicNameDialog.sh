#!/bin/bash

if ! IN=$(zenity --entry --title "nonic ID -> nonic name" --text "Enter ID")
then
  exit;
fi

STR=`nonicName.py --id $IN`

writeString.sh $STR
