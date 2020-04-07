mkdir -p wrapper
cd wrapper

tar -xvvzf subset-bam-1.0-x86_64-linux.tar.gz

BAM=/scratch/bioinfo/thomas/Knoblich.GRP/Abel.Vertesy/114593/114593_premRNA/outs/possorted_genome_bam.bam
CLUSTER=/scratch/bioinfo/thomas/Knoblich.GRP/Abel.Vertesy/114593/114593_premRNA/outs/analysis/clustering/kmeans_6_clusters/clusters.csv

mkdir analysis
awk -vFS="," '{print $1 > "analysis/"$2".txt"}' $CLUSTER
rm analysis/Cluster.txt

for BC in analysis/*.txt; do
    echo $BC
    subset-bam-1.0-x86_64-linux/subset-bam -b $BAM -c $BC -o ${BC/txt/bam}
done