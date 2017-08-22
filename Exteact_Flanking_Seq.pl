use strict;
use warnings;
use POSIX qw(tmpnam);
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);
use Time::localtime;

## ======================================
## Usage: see -h
## ======================================

sub usage{
  warn <<END;
  Usage:
  Run by typing: perl Exteact_Flanking_Seq.pl -genomefile [run directory] -vcffile [vcf file] -len [length (bp)] -outfile [output file]
    Required params:
	-i|vcffile							[s]	VCF-type file in GATK Format
	-o|outfile							[s]	Output file (.fa)
	-g|genomefile						[s]	Genome file (.fa)
	-l|len								[i]	Flanking length on one side (bp)
   
END
  exit;
}
## ======================================
## Get options
## ======================================

my %opt;
%opt = (
    'help'				=> undef,
    'debug'				=> undef,
    'vcffile'		    => undef,
    'outfile'			=> undef,
	'genomefile'		=> undef,
	'len'		=> undef
);

die usage() if @ARGV == 0;
GetOptions (
  'h|help'				=> \$opt{help},
  'debug'				=> \$opt{debug},
  'i|vcffile=s'			=> \$opt{vcffile},
  'o|outfile=s'			=> \$opt{outfile},
  'g|genomefile=s'		=> \$opt{genomefile},
  'l|len=i'				=> \$opt{len}
) or die usage();

#check input paramaters
die usage() if $opt{help};
die usage() unless ( $opt{vcffile} );
die usage() unless ( $opt{outfile} );
die usage() unless ( $opt{genomefile} );
die usage() unless ( $opt{len} );

## ======================================
## Input Genome fasta
## ======================================

print "The input genome file is $opt{genomefile} \n";
my %seq;
open(FA,$opt{genomefile}) or die "No Genome fasta file input \n";
$/=">";
my $i=0;
while(<FA>){
	if (/^(\S+).*?\n(.*)/ms){
		my $key = $1; my $value = $2;
		$value =~ s/>//gm;
		$value =~ s/\*$//gm;
		$value =~ s/\n//gm;
		$key =~ s/\n//mg;
		$seq{$key} = $value;
		$i++;
	}
}
close FA;
$/="\n";
print "The Contigs in $opt{genomefile} is $i \n";

## ======================================
## Input VCF File and extract flanking sequence
## ======================================

print "The input VCF-type file is $opt{vcffile} \n";
open (VCF,$opt{vcffile}) or die "No VCF file input! \n";
open (OUT,">$opt{outfile}");
my $j=0;
while(<VCF>){
	chomp;s/\r//;
	next if /^#/;
	my @g = split /\t/,$_;
	my $chr = $g[0];
	my $snp_pos = $g[1];
	my $genotype  = "\[".$g[3]."\/".$g[4]."\]";
	if (exists $seq{$chr}){
		my $chr_len = length($seq{$chr})-1;
		my $left_start = $snp_pos - $opt{len} - 1;
		my $left_end = $snp_pos - 1 - 1;
		if ($left_start < 0 ) {$left_start = 0};
		my $right_end = $snp_pos + $opt{len} - 1;
		my $right_start = $snp_pos + 1 - 1;
		if ($right_end > $chr_len){$right_end = $chr_len};
		my $left_seq = substr($seq{$chr},$left_start,($left_end-$left_start+1));
		my $right_seq = substr($seq{$chr},$right_start,($right_end-$right_start+1));
		my $snpid = $chr."_".$snp_pos;
		print OUT $snpid."\t".$chr."\t".$snp_pos."\t".$left_seq.$genotype.$right_seq."\n";
		$j++;
	}
}
print "The number of flanking sequences for snp in $opt{vcffile} is $j \n";

