#!/usr/bin/perl -w
#
# Author: Lavanya Rishishwar
# Creation Date  : 17th Sep 2014
# Modified Date  : 19th Sep 2016
#
#############################################################
use strict;
use Getopt::Long;
use threads;
use threads::shared;
#############################################################
my $usage = "Parallel script for computing pairwise ANIm values.\nRequires mummer to be installed.\n\nUsage instructions:\n$0 [-out <output ANI matrix. Default:ani.txt>] [-ext <extension of the fasta files. Default:fasta>] [-threads <number of threads. Default: 10>] [-folder <Folder in which to look for file. Default: ./>] [-distance <FLAG. Prints distances (100-ani) instead of similarity>] [-help <FLAG. Prints the help>]\n";
my $ext="fasta";          # Input Queries Folder
my $out = "ani.txt"; # Output File
my $threads = 10;
my $help = 0;
my $dir = "./";
my $getDistance = 0;
my %ani;
share(%ani);
# Get the arguments
my $args  = GetOptions ("ext=s"    	=> \$ext,
                        "threads=s"	=> \$threads,
                        "folder=s"	=> \$dir,
                        "distance"	=> \$getDistance,
						"help"		=> \$help,
                        "out=s"     => \$out);
my @fastas = <$dir/*.$ext>;
my @names;
foreach (@fastas){
	my $this = `basename $_`;
	$this =~ s/.$ext$//;
	chomp $this;
	push(@names, $this);
}
#@fastas = splice(@fastas, 0, 4);

if($help > 0){
	print STDERR $usage;
	exit 0;
}

if(@fastas == 0){
	print STDERR "ERROR: No FASTA file found.  Looking for *.$ext in the directory: $dir - Exiting!\n\n$usage";
	exit 1;
}

#############################################################
sub runDnaDiff{
	my ($ref, $query, $num) = @_;
	my $this_start = time();
	print "\n  [LOGMSG] : Comparing reference file: $ref with query file: $query\n";
	#my $rand = int(rand(1000));
	#print "dnadiff -p temp-$$ $ref $query\n";
	system("dnadiff -p temp-$$-$num $ref $query");
	my $ani = `head -19 temp-$$-$num.report | tail -1 | awk '{print \$2}'`;
	`rm temp-$$-$num.*`;
	chomp $ani;
	unless(exists $ani{$ref}){
		$ani{$ref} = &share({});
		$ani{$query} = &share({});
	}
	$ani = 100 - $ani if($getDistance > 0);
	$ani{$ref}{$query} = $ani;
	$ani{$query}{$ref} = $ani;
	print STDERR "  [LOGMSG] : $ref-$query ANI = $ani{$ref}{$query}\n";
	my $this_end = time();
	print STDERR "  [LOGMSG] : This took: ".($this_end-$this_start)." sec\n";
}
sub printToFile{
	open OUT, ">$out" or die "Cannot create $out: $!\n";
	for(my $j=0; $j<@names; $j++){
		print OUT "\t$names[$j]";
	}
	print OUT "\n";
	for(my $i=0; $i<@fastas; $i++){
		print OUT "$names[$i]";
		for(my $j=0; $j<@fastas; $j++){
			#$ani{$fastas[$i]}{$fastas[$j]} //= ".";
			print "Checking for $fastas[$i]-$fastas[$j]\n";
			print OUT "\t".$ani{$fastas[$i]}{$fastas[$j]};
		}
		print OUT "\n";
	}
	close OUT;
}
#############################################################
#my @ref = <../ref/*.$ext>;
#print STDERR @ref;
# Populate threads
my @queue;
for(my $i=0; $i<@fastas; $i++){
	$ani{$fastas[$i]} = &share({});
	$ani{$fastas[$i]}{$fastas[$i]} = 100;
	for(my $j=$i+1; $j<@fastas; $j++){
		push(@queue, [$fastas[$i], $fastas[$j]]);
	}
}
# start ani runs
my $start_run = time();
for(my $i=0; $i < @queue; $i++){
	my $this_start = time();
	my $t = async{runDnaDiff($queue[$i][0],$queue[$i][1], $i)};
	my @joinable = threads->list(threads::joinable);
	$_->join() for(@joinable);
	my @running = threads->list(threads::running);
	while(@running == $threads){
		@running = threads->list(threads::running);
		sleep(1);
	}
}
my $end_run = time();
print STDERR "Complete job took: ".($end_run - $start_run)." sec\n";
my @running = threads->list(threads::running);
while(@running > 0){
	my @joinable = threads->list(threads::joinable);
	$_->join() for(@joinable);
	@running = threads->list(threads::running);
	sleep(1);
}
threads->create( \&printToFile)->join;

