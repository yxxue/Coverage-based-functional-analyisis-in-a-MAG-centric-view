# Coverage-based-functional-analyisis-in-a-MAG-centric-view
A example of how to investigate functional potential potential based on coverage with a MAG-centric view.

![Workflow](./functional_analyisis_demo.jpg)

# Introduction

# Example 
There are five samples in svalbard permafrost: sample_ids = ['2-1-2', '2-8-2', '2-9-2', '2-10-2', '2-12-2'].


* Calculation of mapping statistics

```
bbmap.sh ref=all_bins.ctg.fasta in1=2-1-2.R1.fa in2=2-1-2.R2.fa ow=t threads=40 sortscafs=t scafstats=2-1-2.bins.scafstats.txt covstats=2-1-2.bins.covstats.txt statsfile=2-1-2.bins.statsfile.txt 2> 2-1-2.bins.log
```



