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

### 3. Aligment

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