set -euox pipefail

prefix=$1
num_pc=$2
num_sd=$3
LASER_PATH=$4
Rscript_PATH=$5
workflow_PATH=$6
results_PATH=$7

mkdir -p Results
cp $prefix.geno Results/
cp $prefix.site Results/
cd Results
touch removedSamples.txt

for i in {1..30}
do

	# Run PCA for given geno and site file

	$LASER_PATH -g $prefix.geno -k ${num_pc} -pca 1
	# Identify Outliers based on samples having more than X sd in the top K PC's 
	$Rscript_PATH ${workflow_PATH}/removeOutliers.R ${num_pc} ${num_sd} ${results_PATH}
	# Making directory for moving results of ith run
	mkdir -p $i
	cp laser* $i/
	# Break if no more outliers present
	if [ ! -s "remove.txt" ] 
	then
        	break
	fi
	mv laser* $i/
	awk '{print $1}' remove.txt > removeS.txt
	awk 'NR==FNR{a[$0];next} !($2 in a)' removeS.txt $1.geno > $1_tmp.geno

	cp $1_tmp.geno $1.geno
	mv $1_tmp* $i/

	cat removedSamples.txt remove.txt > tmp.txt
	mv tmp.txt removedSamples.txt
	mv remove.txt $i/
	mv removeS.txt $i/

done

mv $1.geno $1.filtered.geno
mv $1.site $1.filtered.site

