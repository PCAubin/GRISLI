#!/usr/bin/env bash

# Download reference database
curl http://www.regulatorynetworks.org/results/networks/mm9/networks.v12032013.tgz -o networks.tgz
tar xvzf networks.tgz

# Combine all cell type-specific regulations together
for d in buffer.5000.mm9-120313/*/ ; do
	cat $d/genes.regulate.genes | cut -f 4,5 >> tempnet
done

sort -u tempnet > regnet

# Create list of transcription factors
cut -f 1 regnet > tflisttmp
cut -f 2 regnet >> tflisttmp
sort -u tflisttmp > tflist

# Clean up
rm tflisttmp
rm tempnet
rm -rf buffer.5000.mm9-120313
rm networks.tgz
