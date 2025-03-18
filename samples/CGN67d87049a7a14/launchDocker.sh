#!/bin/csh
#$ -q local.q
#$ -cwd
#$ -N CGN67d87049a7a14
#$ -e /data/Web/CGN67d87049a7a14/error.log
#$ -o /data/Web/CGN67d87049a7a14/output.log

hostname > hostname.out

docker run --rm -v workflow_data:/mnt -v workflow_scripts:/app/Scripts --cpus "8.00" --memory "8G" workflow_image sh /mnt/Web/CGN67d87049a7a14/launch.sh
