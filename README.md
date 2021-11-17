# Illumina raw data transform to plink form
##### **This Rcode is used for transforming illumina raw gene-chip data to plink formate：made by Hanpeng Luo and email:luohanpeng@cau.edu.cn**

#### Raw data example (illumina):
[Header]
GSGT Version	1.9.4
Processing Date	10/19/2016 12:45 PM
Content		GGP_HDv3_C.bpm
Num SNPs	139376
Total SNPs	140668
Num Samples	307
Total Samples	307
[Data]

| SNP Name | Sample ID |  Allele1 - Forward| Allele2 - Forward |  Allele1 - Top |Allele2 - Top  | Allele1 - AB | Allele2 - AB | GC Score | X |Y  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|ARS-BFGL-BAC-10172|162120236|	G	|G	|G	|G	|B	|B	|0.9420	|0.042	|0.556|
|ARS-BFGL-BAC-1020	|162120236|	G	|G	|G	|G	|B	|B	|0.9489	|0.027	|0.472
|ARS-BFGL-BAC-10245	|162120236|	T	|C	|A	|G	|A	|B	|0.7277	|1.568	|1.245
|ARS-BFGL-BAC-10345	|162120236|	A	|C	|A	|C	|A	|B	|0.9411	|0.680	|0.617
##### Skip first 10 raws,get the data  below [Data] and title, 

#### The data will be transformed in plink form like:
##### First 6 col: family name 
```
TEST  11610  0  0  0  -9  A G A A G G A A A G A G G G C C A A A G 0 0 C G A ...
```

### WHAT SHOULU DO:

```R
#set path, final_report name, map name,out name
setwd("C:\\Users\\lhp\\Desktop\\ped\\")
source("https://raw.githubusercontent.com/hanpengluo/illumina-raw-data-transform/master/illumina_tran_plink.R")
final_report_name<-"Neogen_China_BOVUHDV03_20161018_FinalReport.txt"#finale report file name
map_name<-"SNP_Map.txt"#snp map file name
out_name<-"test" #output file name
#allele coding of snp "Forward","Top","AB"
coding_type<-"Top"
illumina_plink_tran(final_report_name,map_name,out_name,coding_type = "Top")
'''
