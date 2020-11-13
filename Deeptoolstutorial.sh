# DeepTools: tools for exploring deep sequencing data: RNAseq, ChipSeq and MNaseseq
# Install:
pip3 install deeptools 
deeptools --version
# Using pip3 install python packages
# Take the samples file on Github
git clone https://github.com/HGNAIG2020TanTaoUniversity/deeptools.git
# Sample files are take from: GSE113139 for Chipseq bw, Galaxy course for Bam and bed file.
# Download for sample files
https://drive.google.com/drive/folders/1mdgZSzb3XPwJh60depY8VcQfoVb6HsGR?usp=sharing
########################################################################################################################################

# 1.1 MultiBamSummary:
# Sorted and index bam file
for i in *bam
do 
samtools index $i 
done
# Input: indexed bam files- Output:npz (binary) or text (bedgraph), bins (10kb), bed(depend)
# Example
multiBamSummary bins --bamfiles wt_H3K27me3_rep1.bam  wt_H3K4me3_rep1.bam wt_H3K27me3_rep2.bam  wt_H3K4me3_rep2.bam -o results.npz --outRawCounts bin.bedgraph
    # Number of bins found: 272566
multiBamSummary BED-file --BED chr11.bed wt_H3K27me3_rep1.bam  wt_H3K4me3_rep1.bam wt_H3K27me3_rep2.bam  wt_H3K4me3_rep2.bam -o results.npz --outRawCounts scores_per_bin.tab
    # Number of bins found: 2812
# Look at out put with the valid value on each bedgraph
head bin.bedgraph|column -t
    # #'chr'  'start'  'end'  'h3k4me1.bw'         'h3k4me3.bw'          'h3k27ac.bw'
    # chr1    0        10000  0.0915995924346301   0.010728255733366934  0.3699630377970951
    # chr1    10000    20000  7.4212148502272015   1.1757552155525453    7.027449007619939
    # chr1    20000    30000  3.6208818589813854   22.297666320407007    8.761509724017936
    # chr1    30000    40000  1.1712912672393176   2.0783475125220514    2.5744355315302463
    # chr1    40000    50000  0.1226002886013109   0.1416485154674899    0.12407720412133447
    # chr1    50000    60000  0.4560585424423218   0.5194957308369299    0.5843986997684962
    # chr1    60000    70000  0.22172189468072387  0.12271779233870968   0.10268626314082617
    # chr1    70000    80000  0.6058275445743482   0.34281951884608114   0.3365096338997424
    # chr1    80000    90000  0.4944642887349031   0.1764895275485131    0.22226487460069252
head bed.bedgraph|column -t
    # #'chr'  'start'  'end'   'h3k4me1.bw'           'h3k4me3.bw'         'h3k27ac.bw'
    # chr11   116986   121920  0.34431737376708255    0.20028218089398653  0.4964965109887558
    # chr11   183079   184573  0.0007198272819021141  0.0                  0.0
    # chr11   186760   190258  4.349577947602712      22.325290570332704   16.119647024017254
    # chr11   186760   190258  4.349577947602712      22.325290570332704   16.119647024017254
    # chr11   192923   197422  1.866245985488604      0.27205421225604387  1.5285499857303486
    # chr11   192923   197422  1.866245985488604      0.27205421225604387  1.5285499857303486
    # chr11   197510   205175  1.8101928261605906     0.3200449753025351   1.4946615195967201
    # chr11   197510   205175  1.8101928261605906     0.3200449753025351   1.4946615195967201
    # chr11   199335   199406  1.257184517215675      1.1557799577713013   0.6886899622393327
# 1.2. MultiBigwigSummary:
cd /home/giang/Desktop/Deeptools/1.2.multiBigwigSummary
# Example: bin/bed similar to 1.1
multiBigwigSummary bins -b h3k4me1.bw h3k4me3.bw h3k27ac.bw -o results.npz --outRawCounts bin.bedgraph
    # Number of bins found: 313762
multiBigwigSummary BED-file --BED chr11.bed -b h3k4me1.bw h3k4me3.bw h3k27ac.bw -o results.npz --outRawCounts bed.bedgraph
    #
head bin.bedgraph
        #'chr'	'start'	'end'	'h3k4me1.bw'	'h3k4me3.bw'	'h3k27ac.bw'
        # chr1	0	10000	0.0915995924346301	0.010728255733366934	0.3699630377970951
        # chr1	10000	20000	7.4212148502272015	1.1757552155525453	7.027449007619939
        # chr1	20000	30000	3.6208818589813854	22.297666320407007	8.761509724017936
        # chr1	30000	40000	1.1712912672393176	2.0783475125220514	2.5744355315302463
        # chr1	40000	50000	0.1226002886013109	0.1416485154674899	0.12407720412133447
        # chr1	50000	60000	0.4560585424423218	0.5194957308369299	0.5843986997684962
        # chr1	60000	70000	0.22172189468072387	0.12271779233870968	0.10268626314082617
        # chr1	70000	80000	0.6058275445743482	0.34281951884608114	0.3365096338997424
        # chr1	80000	90000	0.4944642887349031	0.1764895275485131	0.22226487460069252
# 1.3 correctGCBias: When the coverage is higher than expected value
# Sort bam files
for i in *bam
do 
samtools sorted -b $i
done
# Example
computeGCBias --numberOfProcessors 3 -b wt_H3K4me3_rep1.bam --effectiveGenomeSize 2913022398 -g mm9.2bit -l 200 --GCbiasFrequenciesFile freq.txt
correctGCBias --numberOfProcessors 3 -b wt_H3K4me3_rep1.bam --effectiveGenomeSize 2913022398 -g mm9.2bit --GCbiasFrequenciesFile freq.txt -o gc_corrected.bam
# Change chromosome notation: 1-> chr1
samtools view -H sorted.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - sorted.bam > test_chr.bam
computeGCBias --numberOfProcessors 3 -b test_chr.bam --effectiveGenomeSize 2913022398 --genome hg38.2bit  -freq freq_test.txt --region X --biasPlot test_gc.png
correctGCBias --numberOfProcessors 4 -b test_chr.bam --effectiveGenomeSize 2913022398 --genome hg38.2bit -freq freq_test.txt -o gc_corrected.bam
        # Sam for chr1 0 184150643 
        # Sam for chr15 0 104043685 
        # Sam for chr9 0 124595110 
        # Sam for chr16 0 98207768 
# --effectiveGenomeSize, define by the mappable genome, depend on references genome for instance:GRCh37	2864785220 GRCh38	2913022398
# --genome: to 2 bit format
# --GS..:Indicate the output file from computeGCBias containing the observed and expected read frequencies per GC-content
# 1.4 BamCoverage
bamCoverage -b reads.bam -of bedgraph -o coverage.bedgraph
### bamCoverage
bamCoverage -b wt_H3K27me3_rep2.bam -o a.SeqDepthNorm.bw \
    --binSize 50
    --normalizeUsing RPKM
    --effectiveGenomeSize 2150570000
    --ignoreForNormalization chrX
    --extendReads
# Examples:
    bamCoverage -b wt_H3K27me3_rep2.bam -of bedgraph -o  a.SeqDepthNorm.bedgraph
# 1.5 BamCompare
bamCompare -p 4 -b1 treatment.bam -b2 control.bam -of bedgraph -o log2ratio.bedgraph
awk '{if($4!=0){print $1,$2,$3,$4}}' log2ratio.bedgraph|head -n20 |column -t
    # X  3050150  3050200  1
    # X  3050200  3050250  1.58496
    # X  3050250  3050300  1
    # X  3050550  3050650  -0.627081
    # X  3050900  3051000  -0.627081
    # X  3051700  3051750  -0.627081
    # X  3051750  3051800  -1.06273
    # X  3051800  3051850  -0.627081
    # X  3051950  3052050  -1.06273
    # X  3052650  3052750  1
    # X  3052950  3053050  1
    # X  3054800  3054900  -0.627081
    # X  3054900  3055000  -1.06273
    # X  3055050  3055150  -0.627081
    # X  3055350  3055550  -0.627081
    # X  3055700  3055900  1
    # X  3057650  3057750  -0.627081
    # X  3057850  3058000  -0.627081
    # X  3058000  3058100  -1.06273
    # X  3058100  3058150  -0.627081

# 1.6 bigWigCompare
bigwigCompare -p 4 -b1 h3k27ac.bw -b2 h3k4me3.bw -of bedgraph -o compare.bedgraph
# 1.7 computeMatrix 
# The out put can be used for plotheamap and plotprofiles tools, 
# Option1: draw from the first point to both side with -b "1000" bases
for i in *bw
do
computeMatrix scale-regions -S $i -R chr11_bed12.bed -b 1000 -o ${i}.gz
done
# Option2: Single input:
computeMatrix reference-point -S *bw -R chr11_bed12.bed -a 3000 -b 3000 -o all
# 1.8 alignSieve
# It will filter (length, mapping quality) to create new bam files with better quality
#Example:
alignmentSieve -p 4 -b wt_H3K27me3_rep2.bam --minMappingQuality 5 --samFlagInclude 16 --samFlagExclude 256 --ignoreDuplicates -o filtered.bam --filterMetrics metrics.txt



########################################################################################################################################
# 2.TOOL FOR QUALATY CONTROL:
# 2.1 PlotCorrelation:
plotCorrelation -in results.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot scatterplot --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.tab
plotCorrelation \
    -in readCounts.npz \
    --corMethod spearman --skipZeros \ #Method: Peason or Spearman
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \ # Chose type of Graph to plot # Color use for heatmap # If set, correlation will be set in the top
    -o heatmap_SpearmanCorr_readCounts.png   \
    --outFileCorMatrix SpearmanCorr_readCounts.tab
# 2.2 PlotPCA:
plotPCA -in readCounts.npz -o PCA_readCounts.png -T "PCA of read counts"
# 2.3 PlotFingerFrinting:
plotFingerprint -b *bam -l h3k4me3r1 h3k4me3r2 wt_H3K27ac_rep1 h3k27ac_rep2 --minMappingQuality 30 --skipZeros --region 19 --numberOfSamples 50000 -T "Fingerprints of different samples" --plotFile fingerprints.png --outRawCounts fingerprints.tab\ # To avoid overlapping bins

# 2.4 BamPEFRagmentSize
bamPEFragmentSize -hist fragmentSize.png -T "Fragment size of PE RNA-seq data" -b *.bam -samplesLabel h3r4me3r1 h3k4me3r2 h3k27me3r1 h3k27me3r2
        # BAM file : wt_H3K27me3_rep1.bam
        # Sample size: 3635

        # Fragment lengths:
        # Min.: 88.0
        # 1st Qu.: 155.0
        # Mean: 206.0346629986245
        # Median: 192.0
        # 3rd Qu.: 241.0
        # Max.: 979.0
        # Std: 68.97124100637404
        # MAD: 42.0
        # Len. 10%: 133.0
        # Len. 20%: 147.0
        # Len. 30%: 162.0
        # Len. 40%: 178.0
        # Len. 60%: 210.0
        # Len. 70%: 230.0
        # Len. 80%: 257.0
        # Len. 90%: 297.0
        # Len. 99%: 423.65999999999985


        # Read lengths:
        # Sample size: 3635

        # Min.: 27.0
        # 1st Qu.: 51.0
        # Mean: 50.76396148555708
        # Median: 51.0
        # 3rd Qu.: 51.0
        # Max.: 51.0
        # Std: 1.3640580845180024
        # MAD: 0.0
        # Len. 10%: 51.0
        # Len. 20%: 51.0
        # Len. 30%: 51.0
        # Len. 40%: 51.0
        # Len. 60%: 51.0
        # Len. 70%: 51.0
        # Len. 80%: 51.0
        # Len. 90%: 51.0
        # Len. 99%: 51.0

# 2.5 ComputeCGbias
computeGCBias -b H3K27Me3.bam --genome hg38.2bit -l 200 -freq freq_test.txt --region X --biasPlot test_gc.png
# 2.6 PlotCovergae
plotCoverage -b *.bam   --plotFile example_coverage -n 1000000 --plotTitle "example_coverage" --outRawCounts coverage.tab --region 19




#########################################################################################################################################
# 3. HEATMAPS AND SUMMARY PLOTS
# 3.1 plotHeatmap
# Use for computematrix scale
for i in *.gz
do
plotHeatmap -m $i -o ExampleHeatmap1.png 
done

# Use for computematrix reference point
# 3.2 plotProfiles
for i in *.gz
do 
plotProfile -m $i
     -out ${i}.png \
     --plotType=fill \ # add color between the x axis and the lines
     --perGroup \ # make one image per BED file instead of per bigWig file
     --colors red yellow blue \
     --plotTitle "Test data profile"
done

#Examples:
for i in *.gz
do
plotProfile -m $i -o ${i}.png --plotType=fill --colors red yellow blue --plotTitle "Test data profile" 
done
# 3.3 plotEnrichment
plotEnrichment -b Input.bam H3K4Me1.bam H3K4Me3.bam \
--BED up.bed down.bed \
--regionLabels "up regulated" "down regulated" \
-o enrichment.png



#########################################################################################################################################
# 4. MISSCELANOUS
# 4.1 computeMatrix Operation
# 4.2 eastimateReadingFiltering