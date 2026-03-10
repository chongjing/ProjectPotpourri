# Stage-specific transcriptome of Bursaphelenchus xylophilus

This project is to generate expression counts from Bursaphelenchus xylopilus transcriptome. Data from [Tanaka et al. (2019)](https://www.nature.com/articles/s41598-019-42570-7; https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=DRP002610&o=acc_s%3Aa) 

### 1. Reference

#### 1.1 Genome
```bash
cd /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus

# genome
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/bursaphelenchus_xylophilus/PRJEA64437/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.genomic.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/bursaphelenchus_xylophilus/PRJEA64437/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa.gz
# annotation
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/bursaphelenchus_xylophilus/PRJEA64437/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.annotations.gff3.gz

# get GO terms for each protein, using interproscan6
cd /home/cx264/project/00.ref/002.Bursaphelenchus_xylophilus_PRJEA64437.WBPS19
export JAVA_HOME=/home/cx264/program/jdk-25.0.2/
nohup nextflow run ebi-pf-team/interproscan6 \
  -profile singularity \
  --input bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa  \
  --datadir interproscan --goterms --pathways 1>interproscan.log 2>&1 &

#re-format GO annotation
awk -F"\t" '$9 ~ "Ontology_term=" {print $1"\t"$9}' bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa.gff3 | sed -e 's/Name=.*.Ontology_term=//g' -e 's/;type=.*//g' | sort -k1,1n | sed -i 's/\.1//g' | awk '!seen[$1,$2]++' > bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.GO.1.tsv
python3 00.re-format.GO.py bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.GO.1.tsv bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.GO.tsv

# index
/data/pathology/program/STAR/bin/Linux_x86_64/STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/index \
  --genomeFastaFiles /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.genomic.fa \
  --sjdbGTFfile /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.annotations.gff3 \
  --sjdbOverhang 99

/data/pathology/program/braker3/braker3.0.8/opt/ETP/tools/gffread bursaphelenchus_xylophilus.PRJEA64437.WBPS19.annotations.gff3 -T -F -o bursaphelenchus_xylophilus.PRJEA64437.WBPS19.annotations.gtf
```

#### 1.2 Annotation
```bash
# SignalP 6.0
nohup /data/pathology/program/Miniforge3/envs/python3.7.12/bin/signalp6 --fastafile bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa --output_dir ./signalp  --format all --organism eukarya --mode slow --write_procs 16 --model_dir /data/pathology/program/SignalP/signalp6_slow_sequential/signalp-6-package/models/ 1>01.signalp.log 2>&1 &

# DeepTMHMM
nohup /data/pathology/program/Miniforge3/envs/python3.7.12/bin/biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa > DeepTMHMM 1>02.DeepTMHMM.log 2>&1 &

# TMHMM2.0
/data/pathology/program/tmhmm-2.0c/bin/tmhmm bursaphelenchus_xylophilus.PRJEA64437.WBPS19.protein.fa > TMHMM2
```

### 2. Data

```bash
cd /data/pathology/cxia/projects/Sebastian/06.Bursaphelenchus_xylophilus/01.data/

for i in {029450..029457} {141203..141222}; do
        echo "Processing DRR${i} "
        /data/pathology/cxia/program/SRAToolkit/sratoolkit.3.0.10-ubuntu64/bin/prefetch DRR${i} && /data/pathology/cxia/program/SRAToolkit/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files --gzip DRR${i} && rm -rf DRR${i}.sra
        output_fwd_paired="DRR${i}_1P.fq.gz"
        output_fwd_unpaired="DRR${i}_1U.fq.gz"
        output_rev_paired="DRR${i}_2P.fq.gz"
        output_rev_unpaired="DRR${i}_2U.fq.gz"
        java -jar /data/pathology/program/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 -summary "DRR${i}.summary" "DRR${i}_1.fastq.gz" "DRR${i}_2.fastq.gz" \
        "$output_fwd_paired" "$output_fwd_unpaired" "$output_rev_paired" "$output_rev_unpaired" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
        rm -rf DRR${i}_1.fastq.gz DRR${i}_2.fastq.gz
        echo "Done DRR${i}"

done
```

### 3. Alignment

```bash
#!/bin/bash
cd /data/pathology/cxia/projects/Sebastian/06.Bursaphelenchus_xylophilus/02.mapping

# we are using STAR for alignment
output_dir="/data/pathology/cxia/projects/Sebastian/06.Bursaphelenchus_xylophilus/02.mapping"

# Run STAR iterate through each sample
for i in {029450..029457} {141203..141222}; do
    cd "${output_dir}"

    # Extract sample name
    sample_name="DRR${i}"

    # Set directory
    mkdir -p "${sample_name}"
    echo "Directory $sample_name created"
    # Change to the newly created directory
    cd "${sample_name}"
    echo "Changed directory to $sample_name"

    # Run STAR
    echo "Processing ${sample_name}"
    ulimit -n 20480
    /data/pathology/program/STAR/bin/Linux_x86_64/STAR \
      --runThreadN 32 \
      --genomeDir /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/index \
      --readFilesIn /data/pathology/cxia/projects/Sebastian/06.Bursaphelenchus_xylophilus/01.data/DRR${i}_1P.fq.gz /data/pathology/cxia/projects/Sebastian/06.Bursaphelenchus_xylophilus/01.data/DRR${i}_2P.fq.gz\
      --readFilesCommand zcat \
      --outFilterMultimapNmax 1 \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${sample_name}

    samtools index ${sample_name}Aligned.sortedByCoord.out.bam

    # Run mapinsights bamqc
    /data/pathology/program/mapinsights/mapinsights/mapinsights bamqc -r /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.genomic.fa -i ${sample_name}Aligned.sortedByCoord.out.bam -o ./
    # Run QualiMap
    /data/pathology/program/QualiMap/qualimap_v2.3/qualimap rnaseq -bam ${sample_name}Aligned.sortedByCoord.out.bam -gtf /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.annotations.gtf --java-mem-size=32G -outdir ./qualimap/ -outformat pdf

    # Get expression counts
    /data/pathology/program/Miniforge3/bin/htseq-count --type transcript --counts_output 3.sorted.bam.count.tsv --nprocesses 32 --max-reads-in-buffer 1000000 ${sample_name}Aligned.sortedByCoord.out.bam /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.annotations.gtf
    # Get expression FPKM TPM
    /data/pathology/program/stringtie-3.0.0.Linux_x86_64/stringtie -p 32 -G /data/pathology/cxia/projects/0.ref/11.Bursaphelenchus_xylophilus/bursaphelenchus_xylophilus.PRJEA64437.WBPS19.annotations.gtf -e -B -A 4.sorted.FPKM.tsv ${sample_name}Aligned.sortedByCoord.out.bam
    
   cd "${output_dir}"

   echo "${sample_name} finished"
done

echo "STAR alignment process completed."
```

### 5. Gene Expression Clustering
#### 5.1.KD0224400_vs_Gfp0224400
```bash
cd /data/pathology/cxia/projects/Sebastian/09.X204SC25076967-Z01-F001/X204SC25076967-Z01-F001_01/06.Clustering_Madalena.AllGene
cat ../06.Clustering_Madalena/KD0224400_vs_Gfp0224400.csv ../06.Clustering_Madalena/KD0795900_vs_Gfp0795900.csv | sort -t',' -k1,1 | uniq -f0 > 00.AllGenes.csv
conda activate R4.3.2
/data/pathology/program/Miniforge3/envs/R4.2.3/bin/R

# standalone node: ClusterGVis in /home/cx264/program/anaconda3/envs/R4.5.1/
```

##### 5.1.2 Clustering for expression matrix
```R
library(ClusterGVis)
library(Biobase)
library(Mfuzz)

FPKM_counts <- read.csv("00.AllGenes.csv", sep=",", row.names = 1)
head(FPKM_counts)
#define a suitable cluster numbers
pdf("01.getClusters.pdf")
getClusters(FPKM_counts)
cm <- clusterData(FPKM_counts, cluster.method = "mfuzz", cluster.num = 8)
ct <- clusterData(FPKM_counts, cluster.method = "TCseq", cluster.num = 8)
ck <- clusterData(FPKM_counts, cluster.method = "kmeans", cluster.num = 8)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ct, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ck, plot.type = "line")
dev.off()
svg("01.getClusters.svg")
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
dev.off()
jpeg("01.getClusters.jpeg", quality=100)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
dev.off()

#line plot
pdf("02.line_heatmap.pdf")
visCluster(object = cm,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
visCluster(object = ct,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
visCluster(object = ck,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
dev.off()
jpeg("02.line_heatmap_cm.jpeg", quality = 100)
visCluster(object = cm,
       plot.type = "both",
       ms.col = c("green","orange","red"),
       column_names_rot = 45)
dev.off()

# save cluster information and membership (if exit)
write.csv(cm$wide.res, "03.cm.8clusters.csv",row.names = TRUE, quote = F)
write.csv(ct$wide.res, "03.ct.8clusters.csv",row.names = TRUE, quote = F)
write.csv(ck$wide.res, "03.ck.8clusters.csv",row.names = TRUE, quote = F)
```

#### 5.2.KD0795900_vs_Gfp0795900
```bash
cd /data/pathology/cxia/projects/Sebastian/09.X204SC25076967-Z01-F001/X204SC25076967-Z01-F001_01/05.Clustering/02.KD0795900_vs_Gfp0795900
awk -F"," 'function abs(x) { return (x < 0) ? -x : x } abs($3)>="0.05" && $7 < "0.05" {print $1}' ../../04.DE_analysis/02.KD0795900_vs_Gfp0795900/03.P-cutoff-0.05.csv | sed 's/"//g' > 001.DE.gene.list
while read line; do awk -v line=$line -F"\t" '$1 == line {print $0}' ../07.FPKM.genes4clustering.tsv >> 002.DEgenes.FPKM.tsv; done < 001.DE.gene.list
conda activate R4.3.2
/data/pathology/program/Miniforge3/envs/R4.2.3/bin/R
```

##### 5.2.2 Clustering for expression matrix
```R
library(ClusterGVis)
library(Biobase)
library(Mfuzz)

FPKM_counts <- read.csv("002.DEgenes.FPKM.tsv", sep="\t", row.names = 1)
FPKM_counts <- FPKM_counts[,c(4,5,6,10,11,12)]
#define a suitable cluster numbers
pdf("003.getClusters.pdf")
getClusters(FPKM_counts)
cm <- clusterData(FPKM_counts, cluster.method = "mfuzz", cluster.num = 4)
ct <- clusterData(FPKM_counts, cluster.method = "TCseq", cluster.num = 4)
ck <- clusterData(FPKM_counts, cluster.method = "kmeans", cluster.num = 4)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ct, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ck, plot.type = "line")
dev.off()


#line plot
pdf("004.line_heatmap.pdf")
visCluster(object = cm,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
dev.off()

# save cluster information and membership (if exit)
write.csv(cm$wide.res, "005.cm.4clusters.csv",row.names = TRUE, quote = F)
write.csv(ct$wide.res, "005.ct.4clusters.csv",row.names = TRUE, quote = F)
write.csv(ck$wide.res, "005.ck.4clusters.csv",row.names = TRUE, quote = F)
```
#### 5.3.AllGenes
```bash
cd /data/pathology/cxia/projects/Sebastian/09.X204SC25076967-Z01-F001/X204SC25076967-Z01-F001_01/05.Clustering/03.AllDEgenes
cat ../01.KD0224400_vs_Gfp0224400/001.DE.gene.list ../02.KD0795900_vs_Gfp0795900/001.DE.gene.list | sort | uniq > 001.DE.AllGenes.list
while read line; do awk -v line=$line -F"\t" '$1 == line {print $0}' ../07.FPKM.genes4clustering.tsv >> 002.DEgenes.FPKM.tsv; done < 001.DE.AllGenes.list
conda activate R4.3.2
/data/pathology/program/Miniforge3/envs/R4.2.3/bin/R
```

##### 5.2.2 Clustering for expression matrix
```R
library(ClusterGVis)
library(Biobase)
library(Mfuzz)

FPKM_counts <- read.csv("00.AllGenes.csv", sep=",", row.names = 1)
# FPKM_counts <- FPKM_counts[,c(4,5,6,10,11,12)]
#define a suitable cluster numbers
pdf("01.getClusters.pdf")
getClusters(FPKM_counts)
cm <- clusterData(FPKM_counts, cluster.method = "mfuzz", cluster.num = 8)
ct <- clusterData(FPKM_counts, cluster.method = "TCseq", cluster.num = 8)
ck <- clusterData(FPKM_counts, cluster.method = "kmeans", cluster.num = 8)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ct, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ck, plot.type = "line")
dev.off()
svg("01.getClusters.svg")
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
dev.off()
png("01.getClusters.png", quality=100)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
dev.off()
jpeg("01.getClusters.jpeg", quality=100)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
dev.off()


#line plot
pdf("02.line_heatmap.pdf")
visCluster(object = cm,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
visCluster(object = ct,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
visCluster(object = ck,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
dev.off()

# save cluster information and membership (if exit)
write.csv(cm$wide.res, "03.cm.8clusters.csv",row.names = TRUE, quote = F)
write.csv(ct$wide.res, "03.ct.8clusters.csv",row.names = TRUE, quote = F)
write.csv(ck$wide.res, "03.ck.8clusters.csv",row.names = TRUE, quote = F)
jpeg("02.line_heatmap_cm.jpeg", quality = 100)
visCluster(object = cm,
       plot.type = "both",
       ms.col = c("green","orange","red"),
       column_names_rot = 45)
dev.off()

```

![Clustering_heatmap](https://github.com/chongjing/ProjectPotpourri/blob/main/01.Bursaphelenchus_RNAseq/Clustering/014.AllGenes.line_heatmap.jpeg).

### 6. GO enrichment
```bash
cd /home/cx264/project/09.BX_Madalena/X204SC25076967-Z01-F001_01/07.GO_enrichment
awk -F"," '$21 > 0.6 {print $2"\t"$20"\t"$21}' ../06.Clustering_Madalena.AllGene/03.cm.8clusters.csv > 000.AllClusters.tsv

#get gene list for each cluster
for i in {1..8}; do awk -v cluster="$i" '$2 == cluster {print $1}' 000.AllClusters.tsv > Cluster.${i}.list; done

# GO emrichment performed using TBtools

# for plot:
sed -e 's/Molecular\ function/MF/g' -e 's/Cellular\ component/CC/g' -e 's/Biological\ process/BP/g' -e 's/corrected\ p-value(BH\ method)/p_adj/g' Repressed\ BXY_0795900.enrich.temp.GO.Enrichment.final.xls | awk -F"\t" 'NR == 1 || ($12 + 0) < 0.05 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$12}' > 01.Repressed_BXY_0795900.GO_enrich.tsv
```

#### 6.1 GO enrichment visualization
```R
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/home/cx264/project/09.BX_Madalena/X204SC25076967-Z01-F001_01/07.GO_enrichment/01.two_TFs")
df <- read.delim("03.Repressed_BXY_0795900.GO_enrich.tsv", header = TRUE, stringsAsFactors = FALSE)
head(df)
df$p_adj[df$p_adj == 0] <- 1e-5
df$neg_log_pval <- -log10(df$p_adj)
df <- df %>%
  group_by(Class) %>%
  arrange(desc(EnrichmentScore)) %>%
  mutate(GO_Name = factor(GO_Name, levels = rev(GO_Name))) %>%  # Reverse for top-to-bottom ordering
  ungroup()
main_plot <- ggplot(df, aes(
  x = EnrichmentScore,
  y = GO_Name,
  size = HitsGenesCountsInSelectedSet,
  color = neg_log_pval
)) +
  geom_point(alpha = 0.85, stroke = 0.3) +
  scale_color_gradient(
    low = "#228B22",   # Green for low significance
    high = "#D73027",  # Red for high significance
    name = "-log10(p_adj)",
    limits = c(0, max(df$neg_log_pval, na.rm = TRUE))
  ) +
  scale_size(
    name = "Gene count",
    range = c(3, 10),
    breaks = pretty(range(df$HitsGenesCountsInSelectedSet), n = 4)
  ) +
  labs(x = "Enrichment score", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),                     # Remove ALL internal grid lines
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),  # Keep outer frame
    axis.text.y = element_text(size = 9, face = "italic"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  )
#main_plot
ggsave("03.Activated_BXY_0795900.GO_enrich.padj_0.05.pdf",
       main_plot,
       width = 10,
       height = 20, limitsize = FALSE,
       dpi = 1200,
       bg = "white")

jpeg("03.Activated_BXY_0795900.GO_enrich.padj_0.05.jpeg", height = 10000, width = 10000, res=1200)
main_plot
dev.off()
```
GO Enrichment plot of Activated BXY_0795900:
<img src="https://github.com/chongjing/ProjectPotpourri/blob/main/01.Bursaphelenchus_RNAseq/GO_enrichment/2_TFs/03.Activated_BXY_0795900.GO_enrich.padj_0.05.jpeg" alt="Image 1" width="600"/>
