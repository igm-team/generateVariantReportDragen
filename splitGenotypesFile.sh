#!/bin/bash

awk 'NR==FNR{a[$2]=1;next}FNR==1{split($0,b,",");for(i=1;i<=length(b);i++){if(b[i] == "Sample Name")col=i};print;FS=","}{if(a[$col])print}' $1 $2 # 11 and 168
