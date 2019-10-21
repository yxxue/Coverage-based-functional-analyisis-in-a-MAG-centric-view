# Coverage-based-functional-analyisis-in-a-MAG-centric-view
A example of how to investigate functional potential potential based on coverage with a MAG-centric view.

![Workflow](./img/functional_analyisis_demo.jpg)

# Introduction
We developed a novel comparative strategy for investigating functional potential based on coverage in a MAG-centric view (See workflow above). 

Generally, metagenomic functional analysis was achieved by mapping short reads or assembled contigs with predicted genes against reference databases followed by parsing the result in gene or pathway level approaches. Gene-by-gene approaches utilizes most dominant gene products while overlooking the fact that biological functions rely on multiple genes while only part of them may be significantly abundant. For another, pathway-level analysis can miss nuanced differences in  functional variance as a core pathway could contain many shared sub-pathways or genes. 

Consequently we deployed comparative analysis strategy that utilizes KEGG Modules, a collection of manually defined functional units with several KO identifiers, which is commonly used to represent functional components. Comparing with pathway or gene enriched analysis, module-based analysis is able to provide a resolution for phenotypic features and link to specific metabolic capacity. Coverage is another important metagenomic characteristic that is currently not used beyond binning assembled contigs into MAGs. Taking all these into account, our novel comparative strategy for investigating functional potential based on coverage with a MAG-centric view captured several trends in Svalbard permafrosts. 

Although we have focused on permafrost metagenomics in this work, similar strategies applied here are applicable to other metagenomic studies, especially for well-characterized environments such as human gut with more accurate taxonomic classification and available MAGs.


# Example

Here we show a example of how to achieve this analysis with Svalbard permafrost metagnome. There are five samples in our dataset: sample_ids = ['2-1-2', '2-8-2', '2-9-2', '2-10-2', '2-12-2']. 


## Calculate mapping statistics (coverage, # of mapped reads, etc)    


Firstly, we merged contigs of all MAGs (all_bins.ctg.fasata), and then mapped using BBMAP to get mapping statistics(covstat, scafstat, statsfile) for each sample (see example below). Please add 'nzo=f' to make sure BBMAP will print contig with all coverage rather than only nonzero coverage.  

```
bbmap.sh ref=all_bins.ctg.fasta in1=2-1-2.R1.fa in2=2-1-2.R2.fa nzo=f ow=t threads=40 sortscafs=t scafstats=2-1-2.scafstats.txt covstats=2-1-2.covstats.txt statsfile=2-1-2.statsfile.txt 2> 2-1-2_map.log
```

## Merge coverage, and normalize coverage with TPM       

'./mapping_stats/' is the storage folder of your mapping statistics files.       
'samples_mapnum.txt' is a two columns tab-separated file: sample_id, total mapped reads (calculate from statsfile)   
'scaf2bin.txt' is a two columns tab-separated file: Contig_ID, Bin_ID 
```
> ls ./mapping_stats
2-1-2.covstats.txt   2-10-2.covstats.txt   2-12-2.covstats.txt   2-8-2.covstats.txt   2-9-2.covstats.txt
2-1-2.scafstats.txt  2-10-2.scafstats.txt  2-12-2.scafstats.txt  2-8-2.scafstats.txt  2-9-2.scafstats.txt
2-1-2.statsfile.txt  2-10-2.statsfile.txt  2-12-2.statsfile.txt  2-8-2.statsfile.txt  2-9-2.statsfile.txt
> head samples_mapnum.txt
2-1-2 9171534
2-8-2 11601008
2-9-2 16526122
2-10-2  18675925
2-12-2  17859538
> head scaf2bin.txt
k127_34665	metabat.113
k127_105473	metabat.113
k127_137477	metabat.113
k127_201593	metabat.115
```
This script will parse all coverage table of your input folder, add Bin_ID, merge into one coverage table, and normalize coverage by TPM caling factor(total mapped read counts divided with 1 million) for each sample.   
It will generate two output files:         
'ctg_bin_cov.csv' is the raw merged coverage table;         
'ctg_bin_cov.norm.csv' is the normalized coverage table.             

```
> python merge_and_norm_cov.py -i ./mapping_stats/ -s samples_mapnum.txt -m scaf2bin.txt

> head ctg_bin_cov.csv
Contig_ID	Bin_ID	2-1-2.raw_fold	2-8-2.raw_fold	2-9-2.raw_fold	2-10-2.raw_fold	2-12-2.raw_fold
k127_105473	metabat.113	3.5285	9.5581	48.0713	6.8607	10.2174
k127_137477	metabat.113	0.0493	19.4976	24.8368	8.5662	45.1643
k127_201593	metabat.113	0.0573	2.5481	52.6891	0.6004	0.72
> head ctg_bin_cov.norm.csv
Contig_ID	Bin_ID	2-1-2.norm_fold	2-8-2.norm_fold	2-9-2.norm_fold	2-10-2.norm_fold	2-12-2.norm_fold
k127_105473	metabat.113	0.38472299181358327	0.8239025436410353	2.9088070389411382	0.36735529833194336	0.5720976656842971
k127_137477	metabat.113	0.0053753276169504465	1.6806815407764566	1.5028813172261466	0.45867607628537815	2.528861608850128
k127_201593	metabat.113	0.00624759173329129	0.21964470673582845	3.1882313346107454	0.03214834071137039	0.040314592684312436

```

## Assign contigs into groups based on their normalized coverage distribution
We pre-defined several groups combining the coverage patterns with geographical significance.
![groups](./img/group.png)


```
Rscript group_ctg_by_norm_cov.R
```

* Calculate module abundance in each group with a MAG-centric view

```
python calculate_module_abundance_per_group.py -gk gene_ko_anno_ghostkoala.txt -ko ko00002.keg -cg NORM_group.csv
```
