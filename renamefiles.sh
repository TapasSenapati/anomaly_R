#!/bin/bash

unzip ./as-733.zip -d ./as-733
unzip ./enron.zip -d ./enron
unzip ./p2p-Gnutella.zip -d ./p2p-Gnutella
unzip ./reality_mining_voices.zip -d reality-mining-voices

for FILE in ./as-733/* ./enron/* ./p2p-Gnutella/* ./reality-mining-voices/*
do
# sed -i 1d $FILE
 newname=`echo $FILE | cut -d'_' -f1`
 mv $FILE $newname
done

