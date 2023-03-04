
#Specify job name
jobname='m'
#Specify directory for the logs
logdir='/directory/for/logs'
#Specify into how many parallel jobs each chr will be conducted with
noBlocks=500

# 
for chr in {1..22} ; do

#Count the number of lines in the bim file, each chr has their own bim file
noLines=$(wc -l < /location/of/bim/file/fileName_c${chr}.bim)
blockSize=$((noLines / noBlocks))

for startVal in $(seq 1 $blockSize $noLines)
do
for endVal in $(seq $blockSize $blockSize $noLines)
do

#conduct only these where endval-startval+1=blocksize
if [ $((endVal - startVal + 1)) == $blockSize ] 
then

output=marg_${chr}_${startVal}_$endVal

echo '#!/bin/bash' > $logdir/$output.sh
echo "
#SBATCH --job-name=${jobname}${chr}_${startVal}_${endVal}
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=12G
#SBATCH --time 1-00:00:00
#SBATCH --output=$logdir/marg_${chr}_${startVal}_${endVal}.log

module add R

cmd=\"Rscript /location/of/scripts/coxExecute.R $chr $startVal $endVal ;\"
echo
echo \$cmd
echo 
\$cmd" >> $logdir/$output.sh
sbatch $logdir/$output.sh



#call one more time if it's the last one
if [ $((endVal + blockSize )) -gt $noLines ]
then
startValLast=$((endVal + 1))

output=marg_${chr}_${startValLast}_$noLines

echo '#!/bin/bash' > $logdir/$output.sh
echo "
#SBATCH --job-name=${jobname}${chr}_${startValLast}_${noLines}
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=12G
#SBATCH --time 1-00:00:00
#SBATCH --output=$logdir/marg_${chr}_${startValLast}_${noLines}.log

module add R

cmd=\"Rscript /location/of/scripts/coxExecute.R $chr $startValLast $noLines ;\"
echo
echo \$cmd
echo 
\$cmd" >> $logdir/$output.sh
sbatch $logdir/$output.sh

fi

fi

done

done

done
