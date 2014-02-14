#!/bin/bash
#this will find all lines of code in the cur dir
#containing <cite> and delete any preceding text
grep '<cite>' *.m|sed 's/^.*<cite>\(.*\)<\/cite>/\1/'

