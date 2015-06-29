#!/bin/bash
wget https://github.com/pratas/smash-phylogenomics/archive/master.zip
unzip master.zip
cd smash-phylogenomics-master/src1
cmake . ; make
cp smash-phylog ../../
cd ../../
rm -rf master.zip smash-phylogenomics-master
