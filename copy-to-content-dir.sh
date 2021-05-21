#!/usr/bin/env bash
#
# copy result figures to docs directory
#
cd docs
for date in `ls ../archives/`
do
  mkdir -p "archives/${date}/Figures/"
  cp ../archives/${date}/Figures/*.png archives/${date}/Figures/
  for region in `ls -d ../archives/${date}/Figures/*/ | xargs -n 1 basename`
  do
    mkdir -p "archives/${date}/Figures/${region}"
    cp ../archives/${date}/Figures/${region}/*.png archives/${date}/Figures/${region}/
  done
done
