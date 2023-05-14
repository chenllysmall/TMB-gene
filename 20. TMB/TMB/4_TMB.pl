
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
	
	#ȥ��������û�иı��
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
	#ͳ��ÿ����Ʒͻ��λ����Ŀ
	$countHash{$sample}++;
	
	#����ͻ��hash
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

#���ͻ�����
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

#���ÿ����Ʒͻ��count
open(WF,">TMB.txt") or die $!;
print WF "id\tTMB\n";
foreach my $key (keys %countHash){
	my $tmb=$countHash{$key}/38;
	print WF "$key\t$tmb\n";
}
close(WF);

#���ÿ������ͻ��count
open(COUNT,">geneMut.txt") or die $!;
print COUNT "Gene\tNum\n";
foreach my $countKey(sort{$geneCountHash{$b}<=>$geneCountHash{$a}} keys(%geneCountHash)){
	print COUNT "$countKey\t$geneCountHash{$countKey}\n";
}
close(COUNT);



