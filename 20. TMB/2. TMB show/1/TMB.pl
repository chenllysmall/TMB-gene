use strict;
#use warnings;

my %hash=();
open(RF, "Trainrisk.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/, $line);
	$hash{$arr[0]}=$arr[$#arr];
}
close(RF);

my %fieldHash=();
my @samp1e=(localtime(time));
my $lineCount=0;
open(RF, "input.maf") or die $!;
open(LOW, ">low.maf") or die $!;
open(HIGH, ">high.maf") or die $!;
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
		print LOW $line . "\n";
		print HIGH $line . "\n";
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
	my @sampleArr=split(/\-/, $sampleName);
	my $subSampleName="$sampleArr[0]-$sampleArr[1]-$sampleArr[2]";
	if(exists $hash{$subSampleName}){
		if($hash{$subSampleName} eq "high"){
			$arr[$fieldHash{"Tumor_Sample_Barcode"}]=$subSampleName;
			print HIGH join("\t", @arr) . "\n";
		}
		if($hash{$subSampleName} eq "low"){
			$arr[$fieldHash{"Tumor_Sample_Barcode"}]=$subSampleName;
			print LOW join("\t", @arr) . "\n";
		}
	}
}
close(HIGH);
close(LOW);
close(RF);

