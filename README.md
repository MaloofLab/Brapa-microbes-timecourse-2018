# Brapa-microbes-timecourse-2018

RNAseq analysis of microbe experiment run by Marc Brock.  

RNA libraries prepared by Amaryllis and sequenced at UCD

## download info

Download link: https://slims.bioinformatics.ucdavis.edu/solexa/view_run.php?id=3587

Julin Downloaded fastq files to Barbera /share/malooflab/Fastqs/Brapa_microbe/timecourse-2018/

Download command was:

    wget -r -nd http://slimsdata.genomecenter.ucdavis.edu/Data/mx7wuv39pp/Unaligned2/Project_JMMC_WYO004/
    
## concatenate reads

On the second download reads had not yet been concatenated across multiple lanes.  Concatenate with

    for f in $(ls *gz | sed s/_L.*// | uniq)
    do
        echo $f
        cat $f* > concatenated/${f}_R1_001.fastq.gz
    done
    
Then delete non-concatenated files and move the concatenated files up a directory level.
    
## sample info

See file `tube_no_legend_time_course_2018.xlsx`

Marc writes: 
>Here’s an excel file of tube_nos and the associated treatments, time points etc.  The only slightly confusing thing is that occasionally a pot assigned to a treatment didn’t have sufficient seedlings etc. and we had to substitute in a backup pot—hence the possible decoupling of pot number and tube number.   All that detail being said—the microbial treatments didn’t change—so this shouldn’t change your analysis IMO.

## JM Object name change on 2/19/20

twoafternoon.any.trtlive.samplingday.lrt.DEGs.all ___IS NOW__ twoafternoon.interaction.trtlive.samplingday.lrt.DEGs.all

dge.diurnal34.anytime.DEGs.all __IS NOW__ dge.diurnal34.interaction.trtlive.time.DEGs.all #interaction

dge.diurnal34.trt.DEGs.all __IS NOW__ dge.diurnal34.trtlive.DEGs.all #treatment

dge.diurnal34.trt.DEGs.all __IS NOW__ dge.diurnal34.any.trtlive.DEGs.all  # treatment and treatment:time interaction

dge.diurnal1314.anytime.DEGs.all __IS NOW__ dge.diurnal1314.interaction.trtlive.time.DEGs

dge.diurnal1314.trt.DEGs.all __IS NOW__ dge.diurnal1314.trtlive.DEGs.all

dge.diurnal1314.trt.DEGs.all __IS NOW__ dge.diurnal1314.any.trtlive.DEGs.all


## importing ASVs (email from Scott Klasek)

```
library("phyloseq")
phyloseq_object <- readRDS(file="path/to/rhizo.ps")
otu_table(phyloseq_object) # the read count data
tax_table(phyloseq_object) # taxonomy data
sample_data(phyloseq_object) # sample information
```