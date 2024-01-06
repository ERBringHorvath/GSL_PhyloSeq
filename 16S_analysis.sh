#!/bin/bash

GSS = Great Salt Lake Sediment

mkdir /srv/data/16S/projects/GSS
cd /srv/data/16S/projects/GSS

forward=$(find /srv/data/16S/raw_data/20190708_16S-V4_PE250/ -maxdepth 1 -type f -regex ".*BRB.*_R1_.*")
reverse=$(find /srv/data/16S/raw_data/20190708_16S-V4_PE250/ -maxdepth 1 -type f -regex ".*BRB.*_R2_.*")
for file in ${forward[@]}; do from=$file; IFS="_ " read -a array <<< $(basename $file); sample=${array[0]}; to="/srv/data/16S/projects/GSS/${sample}.20.forward.fastq.gz"; ln -sv $from $to; done
for file in ${reverse[@]}; do from=$file; IFS="_ " read -a array <<< $(basename $file); sample=${array[0]}; to="/srv/data/16S/projects/GSS/${sample}.20.reverse.fastq.gz"; ln -sv $from $to; done

forward=$(find /srv/data/16S/raw_data/20190708_16S-V4_PE250/ -maxdepth 1 -type f -regex ".*Marina.*_R1_.*")
reverse=$(find /srv/data/16S/raw_data/20190708_16S-V4_PE250/ -maxdepth 1 -type f -regex ".*Marina.*_R2_.*")
for file in ${forward[@]}; do from=$file; IFS="_ " read -a array <<< $(basename $file); sample=${array[0]}; to="/srv/data/16S/projects/GSS/${sample}.20.forward.fastq.gz"; ln -sv $from $to; done
for file in ${reverse[@]}; do from=$file; IFS="_ " read -a array <<< $(basename $file); sample=${array[0]}; to="/srv/data/16S/projects/GSS/${sample}.20.reverse.fastq.gz"; ln -sv $from $to; done

mkdir GSL-sediment
cd GSL-sediment
ln -svi /srv/data/16S/projects/GSS/* ./

srun readstats.py --csv -o readstats.csv ./*.20.* &
less readstats.csv

#! /bin/sh
#SBATCH -J GSS.decontam
reads_path="/home/bbrazelton/GSL-sediment/";
files=$(find -L ${reads_path} -type f -regex ".*\.20.forward\.fastq.*");
for file in ${files[@]}; do
   file=$(basename ${file});
   IFS=". " read -a array <<< ${file};
   sample=${array[0]};
   forward="${reads_path}/${sample}.20.forward.fastq.gz";
   reverse="${reads_path}/${sample}.20.reverse.fastq.gz";
   cutadapt --discard-trimmed --error-rate 0.10 -a ATTAGAWACCCVHGTAGTCCGGCTGACTGACT -A TTACCGCGGCMGCTGGCACACAATTACCATA -g ^TTAGAWACCCVHGTAGTCCGGCTGACTGACT -G ^TACCGCGGCMGCTGGCACACAATTACCATA -o ${sample}.20.forward.decontam.fastq.gz -p ${sample}.20.reverse.decontam.fastq.gz ${forward} ${reverse} > ${sample}.cutadapt.log
done

sbatch decontam.sh

grep "Pairs written (passing filters):" *.log > passed_filter.txt
less passed_filter.txt
