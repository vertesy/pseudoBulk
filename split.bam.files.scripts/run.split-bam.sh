# https://github.com/10XGenomics/subset-bam
# alternative scripts: https://divingintogeneticsandgenomics.rbind.io/post/split-a-10xscatac-bam-file-by-cluster/

# --bam (-b): Input 10x Genomics BAM. This BAM must have the CB tag to define the barcodes of cell barcodes (or the tag defined by --bam-tag). Must also have an index (.bai) file. REQUIRED.
# --cell-barcodes (-c): A cell barcodes file as produced by Cell Ranger that defines which barcodes were called as cells. One barcode per line. In Cell Ranger runs, this can be found in the sub-folder outs/filtered_gene_bc_matrices_mex/${refGenome}/barcodes.tsv where ${refGenome} is the name of the reference genome used in your Cell Ranger run. This file can be used as column labels for the output matrix. REQUIRED.
# --out-bam (-o): A path to write the subsetted BAM file to. REQUIRED.
# --cores: Number of parallel cores to use. DEFAULT: 1.
# --log-level: One of info, error or debug. Increasing levels of logging. DEFAULT: error.
# --bam-tag: Change this to use an alternative tag to the default CB tag. This can be useful for subsetting BAMs from LongRanger.


"on server"

tmux 
cdd /groups/knoblich/users/abel/Data/

# lib="101147"
lib="101146"

bamdir="bam.files.cellranger/"$lib
bamfile=$bamdir/$(ls $bamdir | grep ".bam$")
outdir=$bamdir"/bam.per.cl/"
mkdir $outdir
BCdir="CBCs/"$lib"/"

for BC in $(ls $BCdir | grep ".csv")
do
   BCfile=$BCdir$BC
   echo $BCfile
   echo $bamfile
   outbam=$outdir$BC".bam"
   echo $outbam
   subset-bam --cores 8 --bam $bamfile --cell-barcodes $BCfile --out-bam $outbam --log-level debug
done

# samtools view /Volumes/abel/Data/bam.files.cellranger/101146/bam.per.cl/Cl.1.csv.bam | head
# samtools view /Volumes/abel/Data/bam.files.cellranger/101146/Oli.d110.101146.WT.bam | head -100

