#!/bin/bash
#scp pratas@sapiens.ieeta.pt:/home/pratas/birds5/FB* .
wget https://github.com/pratas/smash-phylogenomics/archive/master.zip
unzip master.zip
cd smash-phylogenomics-master/src1
make clean; make;
cp smash-phylog ../../
cd ../../
(time ./smash-phylog -v -l 14 -t 1.5 -n 4 FB1:FB2:FB3:FB4:FB5:FB6:FB7:FB8:FB9:FB10:FB11:FB12:FB13:FB14:FB15:FB16:FB17:FB18:FB19:FB20:FB21:FB22:FB23:FB24:FB25:FB26:FB27:FB28:FB29:FB30:FB31:FB32:FB33:FB34:FB35:FB36:FB37:FB38:FB39:FB40:FB41:FB42:FB43:FB44:FB45:FB46:FB47:FB48 ) &> REPORT1 &
