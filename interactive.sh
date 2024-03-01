srun -J interactive  --cpus-per-task 16 --ntasks 1 -t 7-0:00 --mem=64000 --partition guru --pty bash
module load samtools/1.9 sratoolkit/2.9.6-1 fastqc/0.11.9 trimmomatic/0.39 reditools/2.0 jacusa/2.0.2 bcftools/1.10 hisat2/2.1.0 bwa/0.7.17 star/2.7.9a

srun -J interactive  --cpus-per-task 1 --ntasks 1 --nodes=1 -t 3-0:00 --mem=64000 --partition mahaguru --nodelist chela-g01 --gres=gpu:tesla:1 --pty bash

