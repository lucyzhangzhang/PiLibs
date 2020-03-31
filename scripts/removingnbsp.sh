#!/bin/bash
#remove that fucking non-breaking white space
cat ../data/pnas_0505266102_05266Table5.csv  | head -n 20 | sed 's/\xc2\xa0//g'
