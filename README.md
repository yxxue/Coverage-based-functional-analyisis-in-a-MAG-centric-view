# Coverage-based-functional-analyisis-in-a-MAG-centric-view
A example of how to investigate functional potential potential based on coverage with a MAG-centric view.

![Workflow](./img/functional_analyisis_demo.jpg)

# Introduction

# Example
There are five samples in svalbard permafrost: sample_ids = ['2-1-2', '2-8-2', '2-9-2', '2-10-2', '2-12-2'].


* Calculate mapping statistics (coverage, # of mapped reads) for each sample

```
bbmap.sh ref=all_bins.ctg.fasta in1=2-1-2.R1.fa in2=2-1-2.R2.fa nzo=f ow=t threads=40 sortscafs=t scafstats=2-1-2.scafstats.txt covstats=2-1-2.covstats.txt statsfile=2-1-2.statsfile.txt 2> 2-1-2_map.log
```

* Merge coverage, and normalize coverage with TPM

```
python merge_and_norm_cov.py -i ./mapping_stats/ -s samples_mapnum.txt -m scaf2bin.txt
```

* Assign contigs into groups based on their normalized coverage distribution
```
Rscript group_ctg_by_norm_cov.R
```

* Calculate module abundance in each group with a MAG-centric view

```
python calculate_module_abundance_per_group.py -gk gene_ko_anno_ghostkoala.txt -ko ko00002.keg -cg NORM_group.csv
```
