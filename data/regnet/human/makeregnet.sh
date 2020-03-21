#!/usr/bin/env bash

# Download reference database
curl http://www.regulatorynetworks.org/results/networks/hg19/networks.v09162013.tgz -o networks.tgz
tar xvzf networks.tgz

# Combine all cell type-specific regulations together
for d in net/lebowski/vol1/work/sjn/papers/regnetworks.resubmission.official/create.networks.per.celltype/results.regulators.combine.all/buffer.5000/*/ ; do
	cat $d/genes.regulate.genes >> tempnet
done

sort -u tempnet > regnet

# Create list of transcription factors
cut -f 1 regnet > tflisttmp
cut -f 2 regnet >> tflisttmp
sort -u tflisttmp > tflist

# Clean up
rm tflisttmp
rm tempnet
rm -rf net
rm networks.tgz
