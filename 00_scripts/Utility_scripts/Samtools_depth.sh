#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/EUostrea
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N bam_depth
#PBS -e 98_log_files/Covstat/__BASE__depth.err
#PBS -o 98_log_files/Covstat/__BASE__depth.out
#PBS -l nodes=1:ppn=8:thinnode
#PBS -l walltime=100:00:00
#PBS -l mem=16gb
#PBS -m n
#PBS -r n


#module
module load tools
module load ngs
module load samtools/1.14

#move to present working dir
cd $PBS_O_WORKDIR

base=__BASE__

#variables
WORKDIR=/home/projects/dp_00007/people/hmon/EUostrea
BAMDIR=/home/projects/dp_00007/people/hmon/Bamfile_EUostrea

samtools depth -aa $BAMDIR/"$base".nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam | cut -f 3 | gzip > $WORKDIR/02_data/Depth/"$base"_depth.gz
