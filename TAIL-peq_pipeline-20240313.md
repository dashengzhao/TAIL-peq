---

**Sequencing data analyses pipeline of TAIL-peq genotyping method**

---

# 1.  Data demultiplexing

## 1.1  To split the raw sequencing data according to the barcodes located at the beginning of the read1 and read2 using fastq-multx tool

```shell
fastq-multx -B ./48_Barcodes.txt -x -m 1 -b Sample01-24-Raw.R1.fq.gz -o ./%.R1.fq
fastq-multx -B ./48_Barcodes.txt -x -m 1 -b Sample01-24-Raw.R2.fq.gz -o ./%.R2.fq
```

The format of 48_Barcodes.txt is listed as follows (the first column represents SP/AD's barcode ID corresponding to the specific individual, and the second column represents the barcode sequence. For example, TL-Hd1_501 and TL-AD2_701 may be corresponded to Sample01, and  TL-Hd1_524 and TL-AD2_724 may be corresponded to Sample24):

```# Document format interpretation
TL-Hd1_501	ACTGAT
TL-Hd1_502	ATGAGC
...
TL-Hd1_523	TCGAAG
TL-Hd1_524	TCGGCA
TL-AD2_701	CGATGT
TL-AD2_702	TGACCA
...
TL-AD2_723	GAGTGG
TL-AD2_724	GGTAGC
```

## 1.2 To get the read heads of fq file corresponding to one individual

```perl
perl head-get.pl  TL-Hd1_501.R1.fq	TL-Hd1_501.R1.fq.head.id
perl head-get.pl  TL-Hd1_501.R2.fq	TL-Hd1_501.R2.fq.head.id
perl head-get.pl  TL-AD2_701.R1.fq	TL-AD2_701.R1.fq.head.id
perl head-get.pl  TL-AD2_701.R2.fq	TL-AD2_701.R2.fq.head.id
```

The contents of head-get.pl are shown as follows:

``` perl
#!/usr/bin/perl -w
#usage: perl $0 input_file  output_file
my (@set);
open IN,$ARGV[0]; # AD1_R03.R1.fq
open OUT,">./$ARGV[1]"; #AD1_R03.R1.fq.head.id
while(<IN>){chomp;
        if($.%4==1){
                @set=split;
                print OUT "$set[0]\n";
                next;
        }else{
                next;
        }
}
close IN;
close OUT;
```

## 1.3 Since each individual corresponds to two barcodes which ligated with SP and AD primer (for example, Sample01 was amplified by TL-Hd1_501 and TL-AD2_701 primer pair), and each barcode could be existed in Read1 or Read2, thus each individual corresponds to four files which storing the heads of corresponding reads. However, some heads are redundant in four files, so we need to merge the four files and obtain the unique heads existing in four files

``` shell
perl R1R2_Head-Merge-Uniq_get.pl TL-Hd1_501.R1.fq.head.id  TL-Hd1_501.R2.fq.head.id TL-AD2_701.R1.fq.head.id  TL-AD2_701.R2.fq.head.id   Sample01-read.head
```

The contents of R1R2_Head-Merge-Uniq_get.pl are shown as follows:

``` perl
#!/usr/bin/perl -w
#usage: perl $0 input_file1  input_file2 input_file3 input_file4  output_file
use List::MoreUtils qw(uniq);
my (@set,@list,@new,$uniq);
open IN0,$ARGV[0];
while(<IN0>){chomp;
	push (@list,$_);
}
close IN0;

open IN1,$ARGV[1];
while(<IN1>){chomp;
	push (@list,$_);
}
close IN1;

open IN2,$ARGV[2];
while(<IN2>){chomp;
	push (@list,$_);
}
close IN2;

open IN3,$ARGV[3];
while(<IN3>){chomp;
	push (@list,$_);
}
close IN3;

@new=uniq @list;
open OUT,">./$ARGV[4]";
for(@new){
	@set=split /\@/;
	print OUT "$set[1]\n";
}
close OUT;
```

## 1.4 To extract the read1 and read2 sequences corresponding to one individual based on the unique heads of reads from raw mixed sequencing data

``` shell
seqkit grep -f Sample01-read.head  Sample01-24-Raw.R1.fq.gz -o Sample01-R1.fq.gz
seqkit grep -f Sample01-read.head  Sample01-24-Raw.R2.fq.gz -o Sample01-R2.fq.gz
```

# 2. Reads composition analysis (optional)

Given that each individual is amplified by a primer pair of SP and AD,  both ligating with distinct barcodes, the reads could be classified into the amplification products of SP+SP, AD+AD, SP+AD, single SP, single AD, and the remaining unknown reads. The proportion of each kind of reads is counted in this part.

## 2.1 To get basic statistics of your raw sequencing data

``` shell
seqkit stat Sample01-Raw.R1.fq.gz
```

You will get the number of heads/sequences. For an example shown below, there are a total of 3,333,333 heads/sequences in Sample01-Raw.R1.fq.gz.

```# output interpretation
file                    format  type   num_seqs      sum_len  min_len  avg_len  max_len
Sample01-Raw.R1.fq.gz  FASTQ   DNA   3,333,333  499,999,950      150      150      150
```

## 2.2 reads splitting

```shell
fastq-multx -B ./Sample01_Barcodes.txt -x -m 1 -b Sample01-Raw.R1.fq.gz -o ./%.R1.fq
fastq-multx -B ./Sample01_Barcodes.txt -x -m 1 -b Sample01-Raw.R2.fq.gz -o ./%.R2.fq
```

The contents of Sample01_Barcodes.txt are as follows (for example, Sample01 is amplified by primer pair of Hd1_B501 and AD2_B701, and the barcodes are B501 and B701, both with length of 6-nt):

```# Document format interpretation
Hd1_B501  ACTGAT
AD2_B701  CGATGT
```

You will get 4 files, named Hd1_B501.R1.fq, Hd1_B501.R2.fq, AD2_B701.R1.fq, and AD2_B701.R2.fq.

## 2.3 To get the read heads of each fq file

(The contents of head-get.pl could be checked in part of **1.2**)

```shell
perl head-get.pl  Hd1_B501.R1.fq	Hd1_B501.R1.fq.head.id
perl head-get.pl  Hd1_B501.R2.fq	Hd1_B501.R2.fq.head.id
perl head-get.pl  AD2_B701.R1.fq	AD2_B701.R1.fq.head.id
perl head-get.pl  AD2_B701.R2.fq	AD2_B701.R2.fq.head.id
```

## 2.4 To get the heads of paired reads amplified by SP+SP or AD+AD

```perl
perl common-get.pl  Hd1_B501.R1.fq.head.id  Hd1_B501.R2.fq.head.id  Hd1_B501.R1-R2.fq.head.id.common
perl common-get.pl  AD2_B701.R1.fq.head.id  AD2_B701.R2.fq.head.id  AD2_B701.R1-R2.fq.head.id.common
```

The contents of common-get.pl are shown as follows:

```perl
#!/usr/bin/perl -w
##usage: perl $0 input_file1  input_file2  output_file
my (%hash0);
open IN0,$ARGV[0];
while(<IN0>){chomp;
	$hash0{$_}=0;
}
close IN0;
open IN1,$ARGV[1];
open OUT,">./$ARGV[2]";
while(<IN1>){chomp;
	if(defined $hash0{$_}){
		print OUT "$_\n";
	}else{
		next;
	}
}
close IN1;
close OUT;
```

## 2.5 To get the heads of paired reads amplified by SP+AD

### 2.5.1 To get the remaining heads of reads with B501 (not include reads from SP+SP):

```shell
perl uniq-get.pl  Hd1_B501.R1-R2.fq.head.id.common  Hd1_B501.R1.fq.head.id    Hd1_B501.R1.fq.head.id.uniq
perl uniq-get.pl  Hd1_B501.R1-R2.fq.head.id.common  Hd1_B501.R2.fq.head.id  Hd1_B501.R2.fq.head.id.uniq
less Hd1_B501.R1.fq.head.id.uniq  Hd1_B501.R2.fq.head.id.uniq >> Hd1_B501.R1-R2.fq.head.id.uniq.merge
```

The contents of uniq-get.pl are shown as follows:

```perl
#!/usr/bin/perl -w
#usage: perl $0  input_file1  input_file2  output_file
my (%hash0);
open IN0,$ARGV[0]; # AD1_R03.R1-R2.fq.head.id.common
while(<IN0>){chomp;
	$hash0{$_}=0;
}
close IN0;
open IN1,$ARGV[1];
open OUT,">./$ARGV[2]";
while(<IN1>){chomp;
	if (defined $hash0{$_}){
		next;
	}else{
		print OUT "$_\n";
	}
}
close IN1;
close OUT;
```

### 2.5.2 To get the remaining heads of reads with B701 (not include reads from AD+AD):

(The contents of uniq-get.pl could be checked in part of **2.5.1**)

```shell
perl uniq-get.pl  AD2_B701.R1-R2.fq.head.id.common  AD2_B701.R1.fq.head.id    AD2_B701.R1.fq.head.id.uniq
perl uniq-get.pl  AD2_B701.R1-R2.fq.head.id.common  AD2_B701.R2.fq.head.id  AD2_B701.R2.fq.head.id.uniq
less AD2_B701.R1.fq.head.id.uniq  AD2_B701.R2.fq.head.id.uniq >> AD2_B701.R1-R2.fq.head.id.uniq.merge
```

### 2.5.3 heads of paired reads amplified by SP+AD:

(The contents of common-get.pl could be checked in part of **2.4**)

```shell
perl common-get.pl  Hd1_B501.R1-R2.fq.head.id.uniq.merge  AD2_B701.R1-R2.fq.head.id.uniq.merge  Hd1_B501-AD2_B701-R1-R2.fq.head.id.uniq.merge.common
```

## 2.6 To get the heads of reads amplified by single SP or single AD primer

(The contents of uniq-get.pl could be checked in part of **2.5.1**)

### 2.6.1 Single SP:

```shell
perl uniq-get.pl Hd1_B501-AD2_B701-R1-R2.fq.head.id.uniq.merge.common  Hd1_B501.R1-R2.fq.head.id.uniq.merge  Hd1_B501-singleStrandPCR.head.id
```

### 2.6.2 Single AD:

```shell
perl uniq-get.pl Hd1_B501-AD2_B701-R1-R2.fq.head.id.uniq.merge.common  AD2_B701.R1-R2.fq.head.id.uniq.merge  AD2_B701-singleStrandPCR.head.id
```

## 2.7 To calculate the proportion of each kind of reads

``` shell
wc -l Hd1_B501.R1-R2.fq.head.id.common  AD2_B701.R1-R2.fq.head.id.common  Hd1_B501-AD2_B701-R1-R2.fq.head.id.uniq.merge.common  Hd1_B501-singleStrandPCR.head.id  AD2_B701-singleStrandPCR.head.id
```

Here is an example of the output:

```# output interpretation
1030	Hd1_B501.R1-R2.fq.head.id.common
459260	AD2_B701.R1-R2.fq.head.id.common
2270008	Hd1_B501-AD2_B701-R1-R2.fq.head.id.uniq.merge.common
120188	Hd1_B501-singleStrandPCR.head.id
462815	AD2_B701-singleStrandPCR.head.id
3313301	total
```

Thus, the proportion could be calculated by :

SP+SP: 1030\*2/(3333333\*2)\*100%=0.03%;

AD+AD: 459260\*2/(3333333\*2)\*100%=13.78%;

SP+AD: 2270008\*2/(3333333\*2)\*100%=68.10%;

Single SP: 120188/(3333333\*2)\*100%=1.80%;

Single AD: 462815/(3333333\*2)\*100%=6.94%;

Unknown reads: 1-0.03%-13.78%-68.10%-1.80%-6.94%=9.35%.

# 3. Reads subsampling (optional)

To randomly subsample ~1.0 Gb sequencing reads from raw data of one individual if necessary, you can use seqtk tool, and the command is listed as follow:

``` shell
seqtk sample -s100 individual-1.R1.fq.gz 3333333|gzip >individual-1_1Gb.R1.fq.gz
seqtk sample -s100 individual-1.R2.fq.gz 3333333|gzip >individual-1_1Gb.R2.fq.gz
```

# 4. Barcode trimming

 Since the length of barcode seuqence is 6-nt, so we remove the first 6-nt bases at the beginning of Read1 and Read2:

```shell
fastp -i Sample01-R1.fq.gz -I Sample01-R2.fq.gz -f 6 -F 6  -A -L -Q -o Sample01-clean1-R1.fq.gz -O Sample01-clean1-R2.fq.gz 
```

# 5. Reads quality controlling

To remove the low-quality bases or reads in barcode-removed raw data by fastp tool (if you want to detect PTMA tags, we suggest you to skip this step). The corresponding command is listed as follow:

``` shell
fastp -i Sample01-clean1-R1.fq.gz -I Sample01-clean1-R2.fq.gz -q 20 -l 15 -o Sample01-clean2-R1.fq.gz -O Sample01-clean2-R2.fq.gz
```

# 6. Reads mapping and sorting

To align the clean reads to the corresponding reference genome. Here is an example:

```shell
#for downstream tag identification:
bwa mem -t 3 -M reference_genome.fa Sample01-clean1-R1.fq.gz Sample01-clean1-R2.fq.gz|samtools view -bS -q 30 --threads 3 >Sample01-clean1.bam; samtools sort -m 20000000000 Sample01-clean1.bam > Sample01-clean1.sort.bam; samtools index Sample01-clean1.sort.bam
#for downstream SNP identification:
bwa mem -t 3 -M reference_genome.fa Sample01-clean2-R1.fq.gz Sample01-clean2-R2.fq.gz|samtools view -bS -q 30 --threads 3 >Sample01-clean2.bam; samtools sort -m 20000000000 Sample01-clean2.bam > Sample01-clean2.sort.bam; gatk AddOrReplaceReadGroups -I Sample01-clean2.sort.bam -O Sample01-clean2.sort.RGadd.bam --RGLB Sample01 --RGPL ILLUMINA --RGPU Sample01 --RGSM Sample01; samtools index Sample01-clean2.sort.RGadd.bam 
```

# 7. Tag identification

## 7.1 To calculate the coverage depth of each base across the genome

``` shell
samtools depth Sample01-clean1.sort.bam > Sample01-clean1.sort.bam.dep
```

## 7.2 To find potential tags (coverage depth ≥ 3×) across the whole genome

```shell
perl find_tags.pl -i Sample01-clean1.sort.bam.dep -d 3 -o Sample01-clean1.sort.bam.dep.tags1
```

The contents of find_tags.pl are shown as follows:

```perl
#!/usr/bin/perl -w
#usage: perl $0 -i input_file  -d 3 -o output_file
#-i: input
#-d: the depth of tags site should bigger than this value.)
#-o: output
my(@list,%hash1,%hash2,%hash3,@row,$name,$hou1,$all_line,$max,@hang,$qian1);
use Getopt::Std;
use vars qw($opt_i  $opt_d $opt_o);
getopts ('i:d:o:');
open IN,$opt_i;
while(<IN>){chomp;
	@list=split;
	$pos=$list[0]."-".$list[1]; 
	$hash1{$.}=$list[2];
	$hash2{$.}=$_;      
	$hash3{$pos}=$_;    
}
close IN;
open OUT,">./$opt_o";
print OUT "#Chr\tPos(1-based)\tdepth(x)\n";
my $line=1;
$all_line=keys %hash2;
while($line<$all_line){ 
	if($hash1{$line}>=$opt_d){ 
		@list=split /\s+/,$hash2{$line}; 
		$hou1=$list[1]-1;                
		$name=$list[0]."-".$hou1;        
		$qian1=$list[1]+1;
		$id=$list[0]."-".$qian1;        
		if(defined $hash3{$name}){       
			@row=split /\s+/,$hash3{$name};
		}else{
			$row[2]=0; 
			$row[0]=$list[0]; 
		}

		if(defined $hash3{$id}){    
			@hang=split /\s+/,$hash3{$id};
		}else{
			$hang[2]=0; 
			$hang[0]=$list[0]; 
		}

		if(($list[0] eq $row[0]) and  ($list[2]>=$row[2]+$opt_d) ){
			print OUT "$hash2{$line}\n";
			$line++;
		}elsif(($list[0] eq $hang[0]) and ($list[2]>=$hang[2]+$opt_d)){
			print OUT "$hash2{$line}\n";
			$line++;
		}else{
			$line++;
		}
	}else{
		$line++;
	}
}
close OUT;
```

## 7.3 To remove false tags (mainly caused by InDels) detected in the previous step and also listed the CIGAR values of each mapped reads beginning at the corresponding tags

```shell
perl find_readsOnTags.pl Sample01-clean1.sort.bam.dep.tags1  Sample01-clean1.sort.bam  Sample01-clean1.sort.bam.dep.tags2
less Sample01-clean1.sort.bam.dep.tags2|awk '{print $1"-"$2}' > Sample01-clean1.sort.bam.dep.tags2.cord
```

The contents of find_readsOnTags.pl are shown as follows:

```perl
#!/usr/bin/perl -w
#usage: perl $0 input_file1  input_file2  output_file
#
my (@list,$id,%hash,@row,$name,$pos,$samp);
open IN0,"$ARGV[0]";
while(<IN0>){chomp;
	next if /^#/;
	@list=split;
	$id=$list[0]."-".$list[1];
	$hash{$id}=0;
	
}
close IN0;
open IN1,"samtools view $ARGV[1]|";
while(<IN1>){chomp;
	next if /^#/;
	@row=split;
	if($row[8]>0 and $row[8]<2000){
		$name=$row[2]."-".$row[3];
		if(defined $hash{$name}){
			push(@{$name},$row[5]);
		}else{
			next;
		}
	}elsif($row[8]<0){
		$pos=$row[7]-$row[8]-1;
		$samp=$row[2]."-".$pos;
		if(defined $hash{$samp}){
			push(@{$samp},$row[5]);
		}else{
			next;
		}
	}else{
		next;
	}
}
close IN1;
open IN0,"$ARGV[0]";
open OUT,">./$ARGV[2]";
while(<IN0>){chomp;
	next if /^#/;
	@list=split;
	$id=$list[0]."-".$list[1];
	if(@{$id} !=()){
		$read=join(",",@{$id});
		print OUT "$list[0]\t$list[1]\t$read\n" if (@{$id} >=3);
		$read="";
	}else{
		next;
	}
}
close IN0;
close OUT;

```

## 7.4 To count the number of tags:

``` shell
wc -l Sample01-clean1.sort.bam.dep.tags2.cord
```

The format of Sample01-clean1.sort.bam.dep.tags2.cord is shown as follows:

``` # Document format interpretation
Chr1-1326
Chr1-13043
Chr1-14587
...
```

## 7.5 To count the number of overlapped tags among three replicates:

```shell
perl shared_tags_analysis.pl -a Sample01-clean1-rep1.sort.bam.dep.tags2.cord  -b Sample01-clean1-rep2.sort.bam.dep.tags2.cord  -c Sample01-clean1-rep3.sort.bam.dep.tags2.cord  -o  Sample01-clean1-Rep123-dep3-tags-sharedAnalysis.out
```

The contents of shared_tags_analysis.pl are shown as follows:

```perl
#!/usr/bin/perl -w
#usage: perl $0  -a input_file1 -b input_file2 -c input_file3 -o output_file
#计算三个输入文件之间的相互重叠关系，以便venn图制作。
use Getopt::Std;
use vars qw($opt_a $opt_b $opt_c $opt_o);
getopts ('a:b:c:o:');
use List::Util qw /uniq/;
my (%hash0,%hash1,%hash2,@all,@uniq,$a,$b,$c,$ab,$ac,$bc,$abc,$uniq_a,$uniq_b,$uniq_c);
open IN0,$opt_a;
while(<IN0>){chomp;
	$hash0{$_}=1;
	push(@all,$_);
}
close IN0;
open IN1,$opt_b;
while(<IN1>){chomp;
	$hash1{$_}=2;
	push(@all,$_);
}
close IN1;
open IN2,$opt_c;
while(<IN2>){chomp;
	$hash2{$_}=3;
	push(@all,$_);
}
close IN2;

$a=keys %hash0;
$b=keys %hash1;
$c=keys %hash2;
$abc=0;
$ab=0;
$ac=0;
$bc=0;
@uniq=uniq @all;
my $all_num=@uniq;
for (@uniq){
	if( (defined $hash0{$_}) and (defined $hash1{$_}) and (defined $hash2{$_})  ){
		$abc ++;
	}
	if((defined $hash0{$_}) and (defined $hash1{$_}) and ($hash2{$_} !=3)){
		$ab ++;
	}
	if(($hash0{$_} !=1) and (defined $hash1{$_}) and (defined $hash2{$_})){
		$bc ++;
	}
	if((defined $hash0{$_}) and ($hash1{$_} !=2) and (defined $hash2{$_})){
		$ac ++;
	}
}
$uniq_a=$a-$abc-$ab-$ac;
$uniq_b=$b-$abc-$bc-$ab;
$uniq_c=$c-$abc-$ac-$bc;
open OUT,">./$opt_o";
print OUT "The number of total tags in a file is: $a.\n";
print OUT "The number of total tags in b file is: $b.\n";
print OUT "The number of total tags in c file is: $c.\n";
print OUT "The number of shared tags num in abc three files is : $abc.\n";
print OUT "The number of uniq tags in a file is : $uniq_a.\n";
print OUT "The number of uniq tags in b file is : $uniq_b.\n";
print OUT "The number of uniq tags in c file is : $uniq_c.\n";
print OUT "The number of shared tags in ab two files is : $ab.\n";
print OUT "The number of shared tags in bc two files is : $bc.\n";
print OUT "The number of shared tags in ac two files is : $ac.\n";
print OUT "The number of uniq tags in abc three files num is: $all_num.\n\n";
close OUT;

```

# 8. SNP calling and filtering

## 8.1 To get SNPs on individual-level

### 8.1.1 To get all SNPs across the whole genome

```shell
# The filtering parameters for SNPs could be adjusted according to your specific research. Here is an example:
gatk --java-options "-Xmx6G" HaplotypeCaller --minimum-mapping-quality 30  --max-reads-per-alignment-start 1000 --min-base-quality-score 20 -R  reference_genome.fa -I Sample01-clean2.sort.RGadd.bam --emit-ref-confidence GVCF  -O Sample01.gvcf
gatk GenotypeGVCFs --call-genotypes true -R reference_genome.fa -V Sample01.gvcf -O Sample01.vcf
gatk SelectVariants -R reference_genome.fa --select-type-to-include SNP -V Sample01.vcf -O Sample01-snp.vcf
gatk VariantFiltration -R reference_genome.fa  -V Sample01-snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O ./Sample01-snp.vcf.filtering
less Sample01-snp.vcf.filtering |grep -v "Filter" > Sample01-snp.vcf.filtered
```

### 8.1.2 To get all SNPs based on specific genomic coordinates

```shell
#Method_1 (The GATK engine recognizes the .bed extension and interprets the coordinate system accordingly.However, you should be aware that this file format is 0-based for the start coordinates, so coordinates taken from 1-based formats [e.g. if you're cooking up a custom interval list derived from a file in a 1-based format] should be offset by 1):
gatk  --java-options "-Xmx6G" HaplotypeCaller -R reference_genome.fa -I Sample01-clean2.sort.RGadd.bam  -O test2.vcf -L pos.bed

#Method_2:
samtools mpileup Sample01-clean2.sort.RGadd.bam -g -Q 20 -q 30 -l Specific_pos.list -f reference_genome.fa --VCF -o Sample01-clean2.vcf
```

The format of pos.bed (0-based):

```# Document format interpretation
Chr1	1248	1249
Chr1	1276	1277
...
```

The format of Specific_pos.list (1-based):

```# Document format interpretation
Chr1	1249
Chr1	1277
Chr1	1278
Chr1	1297
Chr1	1299
...
```

## 8.2 To get SNPs on population-level

### 8.2.1 To get the GVCF file of each individual

```shell
gatk  --java-options "-Xmx6G" HaplotypeCaller -R reference_genome.fa --emit-ref-confidence GVCF -I Sample01-clean2.sort.RGadd.bam -O Sample01.sort.RGadd.bam.gvcf
...
gatk  --java-options "-Xmx6G" HaplotypeCaller -R reference_genome.fa --emit-ref-confidence GVCF -I Sample100-clean2.sort.RGadd.bam -O Sample100.sort.RGadd.bam.gvcf
```

### 8.2.2 To merge all the GVCF files into one

```shell
gatk CombineGVCFs -R reference_genome.fa -V Sample01.sort.RGadd.bam.gvcf -V Sample02.sort.RGadd.bam.gvcf -V Sample03.sort.RGadd.bam.gvcf ... -V Sample100.sort.RGadd.bam.gvcf  -O Sample01-100-CombineGVCFs.gvcf
```

### 8.2.3 To call variations from GVCF file:

```shell
gatk GenotypeGVCFs --call-genotypes true -R reference_genome.fa  -V Sample01-100-CombineGVCFs.gvcf  -O Sample01-100-CombineGVCFs.vcf
```

### 8.2.4 To select raw SNPs:

```shell
gatk SelectVariants -R reference_genome.fa -select-type-to-include SNP -V Sample01-100-CombineGVCFs.vcf -O Sample01-100-CombineGVCFs.vcf.snp
```

### 8.2.5 To get high-quality SNPs:

```shell
# The filtering parameters for SNPs could be adjusted according to your specific research. Here is an example:
gatk VariantFiltration -R reference_genome.fa -V ./Sample01-100-CombineGVCFs.vcf.snp --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O ./Sample01-100-CombineGVCFs.vcf.snp.filtering
less Sample01-100-CombineGVCFs.vcf.snp.filtering |grep -v "Filter" > Sample01-100-CombineGVCFs.vcf.snp.filtered
vcftools --vcf Sample01-100-CombineGVCFs.vcf.snp.filtered --max-missing 0.8 --minQ 30 --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --out Sample01-100-HighQuality-snp.vcf
```

After getting the VCF file including the high-quality SNPs, you can perform the subsequent population genetics analyses, such as PCA, Phylogenetic tree construction, Population admixture inference, LD decay analysis, GWAS, etc.