#!/bin/csh
#$ -q local.q
#$ -cwd
#$ -N CGN67d9bc14bc723
#$ -e /data/Web/CGN67d9bc14bc723/error.log
#$ -o /data/Web/CGN67d9bc14bc723/output.log

hostname > hostname.out

docker run --rm -v workflow_data:/mnt -v workflow_scripts:/app/Scripts --cpus "8.00" --memory "8G" workflow_image sh /mnt/Web/CGN67d9bc14bc723/launch.sh
