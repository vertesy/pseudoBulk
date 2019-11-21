# mv 'Galaxy11-[CellRanger_count_BAM_file_for_data_5,_data_4,_and_data_7].bam' Oli.d110.101147.TSC2.bam 
# mv 'Galaxy15-[CellRanger_count_BAM_file_for_data_3,_data_2,_and_data_6].bam' Oli.d110.101146.WT.bam 

cd /Volumes/groups/knoblich/users/abel/Data/bam.files.cellranger/

subset-bam --cores 4 --bam 101146/Oli.d110.101146.WT.bam --cell-barcodes 101146/BCs/BCs.101146.cl1.tsv --out-bam 101146/101146.cl1.bam --log-level info

# --bam (-b): Input 10x Genomics BAM. This BAM must have the CB tag to define the barcodes of cell barcodes (or the tag defined by --bam-tag). Must also have an index (.bai) file. REQUIRED.
# --cell-barcodes (-c): A cell barcodes file as produced by Cell Ranger that defines which barcodes were called as cells. One barcode per line. In Cell Ranger runs, this can be found in the sub-folder outs/filtered_gene_bc_matrices_mex/${refGenome}/barcodes.tsv where ${refGenome} is the name of the reference genome used in your Cell Ranger run. This file can be used as column labels for the output matrix. REQUIRED.
# --out-bam (-o): A path to write the subsetted BAM file to. REQUIRED.
# --cores: Number of parallel cores to use. DEFAULT: 1.
# --log-level: One of info, error or debug. Increasing levels of logging. DEFAULT: error.
# --bam-tag: Change this to use an alternative tag to the default CB tag. This can be useful for subsetting BAMs from LongRanger.

# https://github.com/10XGenomics/subset-bam
# alternative scripts: https://divingintogeneticsandgenomics.rbind.io/post/split-a-10xscatac-bam-file-by-cluster/

tmux 

"on server"
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
   subset-bam --cores 4 --bam $bamfile --cell-barcodes $BCfile --out-bam $outbam --log-level debug
done

# samtools view /Volumes/abel/Data/bam.files.cellranger/101146/bam.per.cl/Cl.1.csv.bam | head
# samtools view /Volumes/abel/Data/bam.files.cellranger/101146/Oli.d110.101146.WT.bam | head -100

# samtools view -X 101146.cl1.bam | less -S
# samtools view -X 101146.cl1.bam | tail
cd bam.files.cellranger
subset-bam --cores 4 --bam 101146/Oli.d110.101146.WT.bam --cell-barcodes 101146/BCs/BCs.101146.cl1.tsv --out-bam 101146/101146.cl1.bam --log-level debug
subset-bam --cores 4 --bam 101146/Oli.d110.101146.WT.bam --cell-barcodes ../CBCs/101146/Cl.9.csv --out-bam 101146/101146.cl1.bam --log-level debug


subset-bam --cores 4 --bam $bamfile --cell-barcodes 101147/BCs/BCs.101147.cl1.tsv --out-bam 101147/101147.cl1.bam --log-level info
