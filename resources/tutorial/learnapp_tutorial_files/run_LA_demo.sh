# !/bin/bash

# reset files
mv input_LA/zipped/*.gz input_LA/
rm -rf learn apply 


# make run directories
mkdir -p learn/input apply/input

# copy annotations and inputs to local dir
cp -r ../../annotations learn 
cp -r ../../annotations apply 


cp inputs/UP000317332_2567861.fasta inputs/UP000319088_2600309.fasta inputs/UP000319639_2592816.fasta inputs/UP000319825_1880.fasta inputs/UP000480288_2708300.fasta inputs/UP000481030_1602942.fasta inputs/UP000317977_2528013.fasta inputs/UP000319374_2585119.fasta learn/input

cp inputs/UP000319776_92403.fasta inputs/UP000319897_1715348.fasta inputs/UP000480297_2681465.fasta inputs/UP000482960_1076125.fasta apply/input


# run snekmer cluster on examples using provided config.yaml and return to base dir
cd learn
snekmer learn --configfile=../../../LA_config.yaml
cd ..

# make new directories required for apply input
mkdir -p apply/counts apply/confidence

# move learn model and confidence files to above directories
cp learn/output/learn/kmer-counts-total.csv apply/counts
cp learn/output/eval_conf/global-confidence-scores.csv apply/confidence

# run snekmer model on examples using provided config.yaml and return to base dir
cd apply
snekmer apply --configfile=../../../LA_config.yaml
cd ..
