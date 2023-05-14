
use strict;
#use warnings;

my $file="input.maf";
my %hash=();
my @sampleArr=();
my %sampleHash=();
my $gene="all";
my %fieldHash=();
my %countHash=();

my $lineCount=0;
my @samp1e=(localtime(time));
open(RF,"$file") or die $!;
while(my $line=<RF>){
	next if($line=~/^\n/);
	next if($line=~/^\#/);
	$lineCount++;
	chomp($line);
	my @arr=split(/\t/,$line);
	if($lineCount==1){
		for(my $i=0;$i<=$#arr;$i++){
			$fieldHash{$arr[$i]}=$i;
		}
		next;
	}
	
	#去除氨基酸没有改变的
	if($arr[$fieldHash{"Amino_acids"}] eq ""){
		next;
	}
	if($arr[$fieldHash{"Variant_Classification"}] eq "Silent"){
		next;
	}
	if($arr[$fieldHash{"Variant_Classification"}] eq "Splice_Region"){
		next;
	}
	
	my $sampleName=$arr[$fieldHash{"Tumor_Sample_Barcode"}];
	my @sampleNameArr=split(/\-/, $sampleName);
	my $sample="$sampleNameArr[0]-$sampleNameArr[1]-$sampleNameArr[2]";
	#统计每个样品突变位点数目
	$countHash{$sample}++;
	
	#构建突变hash
	unless($sampleHash{$sample}){
		push(@sampleArr,$sample);
		$sampleHash{$sample}=1;
	}
	my $geneField=$fieldHash{"Hugo_Symbol"};
	if( ($gene eq "all") || ($gene eq $arr[$geneField]) ){
		$hash{$arr[$geneField]}{$sample}=1;
	}
}
close(RF);

#输出突变矩阵
my %geneCountHash=();
open(WF,">mutMatrix.txt") or die $!;
print WF "Gene\t" . join("\t",@sampleArr) . "\n";
foreach my $key(keys %hash){
	print WF $key;
	foreach my $sample(@sampleArr){
		if(exists ${$hash{$key}}{$sample}){
			print WF "\t1";
			$geneCountHash{$key}++;
		}
		else{
			print WF "\t0";
		}
	}
	print WF "\n";
}
close(WF);

#输出每个样品突变count
open(WF,">TMB.txt") or die $!;
print WF "id\tTMB\n";
foreach my $key (keys %countHash){
	my $tmb=$countHash{$key}/38;
	print WF "$key\t$tmb\n";
}
close(WF);

#输出每个基因突变count
open(COUNT,">geneMut.txt") or die $!;
print COUNT "Gene\tNum\n";
foreach my $countKey(sort{$geneCountHash{$b}<=>$geneCountHash{$a}} keys(%geneCountHash)){
	print COUNT "$countKey\t$geneCountHash{$countKey}\n";
}
close(COUNT);



