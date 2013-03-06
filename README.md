## OVERVIEW
This is the "lean" development branch, optimized for pan-genomic analyses of large numbers of sequences.

Panseq determines the core and accessory regions among a collection of genomic sequences based on user-defined parameters. It readily extracts regions unique to a genome or group of genomes, identifies SNPs within shared core genomic regions, constructs files for use in phylogeny programs based on both the presence/absence of accessory regions and SNPs within core regions

If you find Panseq useful please cite:

> Pan-genome sequence analysis using Panseq: an online tool for the rapid analysis of core and accessory genomic regions.
> Laing C, Buchanan C, Taboada EN, Zhang Y, Kropinski A, Villegas A, Thomas JE, Gannon VP.
> BMC Bioinformatics. 2010 Sep 15;11:461.

> Identification of Salmonella enterica species- and subgroup-specific genomic regions using Panseq 2.0.
> Laing C, Villegas A, Taboada EN, Kropinski A, Thomas JE, Gannon VP.
> Infect Genet Evol. 2011 Dec;11(8):2151-61.


## USAGE

The Panseq standalone script can be accessed from:

> lib/panseq.pl


## Setup

Panseq requires Perl 5.10 or greater, the provided Perl modules, and the following CPAN packages:

* Parallel::ForkManager
* Log::Log4perl
* Role::Tiny
* Bio::SeqIO
* Bio::DB::Fasta
* Bio::Seq
* Tie::Log4perl
* Statistics::R

To install, and automatically retrieve the CPAN packages, do the following:

	perl Build.PL
	./Build installdeps


The following free, external programs must also be installed:

* [MUMmer:](http://mummer.sourceforge.net)
* [BLAST+:](http://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [Muscle:](www.drive5.com/muscle)

## Running Panseq

All the adjustments to Panseq are made by modifying a tab-delimited configuration file, which is specified as the only argument to the script.

	perl panseq.pl settings.txt


Below is an example configuration file:

	queryDirectory	/home/phac/panseq/queryLarge/
	referenceDirectory	/home/phac/panseq/referenceLarge/
	baseDirectory	/home/phac/panseq/output/panseq2test/
	numberOfCores	22
	mummerDirectory	/home/phac/MUMmer3.23/
	blastDirectory	/home/phac/ncbi-blast-2.2.25+/bin/
	minimumNovelRegionSize	500
	novelRegionFinderMode	no_duplicates
	muscleExecutable	/usr/bin/muscle3.8.31_i86linux64
	fragmentationSize	500
	percentIdentityCutoff	85
	coreGenomeThreshold	2
	mode 	pan


* ‘queryDirectory’ should contain the full directory path of the folder where all of the query sequences you are interested in comparing reside. Panseq will use the entire contents of this folder. 

* ‘referenceDirectory’ should contain the full directory path of the folder where all of the reference sequences you are interested in comparing reside. Panseq only uses this folder for Novel Region Comparisons of type 'common_to_all' and 'no_duplicates'. All other analyses use the contents of the Query folder only.

* ‘baseDirectory’ is the directory where all the output from Panseq is placed, and should be the full directory path. WARNING: This directory is automatically cleared of all contents at the beginning of each run.

* 'numberOfCores' sets the number of processors available to Panseq. Increasing this can reduce run times.

* 'mummerDirectory' specifies the full path to the folder containing the nucmer program.

* 'blastDirectory' specifies the full path to the BLAST+ bin directory.

* 'minimumNovelRegionSize' sets the size in bp of the smallest region that will be kept by the Novel Region Finder; all regions found below this value will not be kept.

* 'novelRegionFinderMode'  sets the type of novel region analysis that will be performed. ‘no_duplicates’ finds the novel regions among one or more query strains with respect to the reference strains selected. ‘unique’ finds sequence regions that are unique to each sequence among all of the strains selected; only strains in the 'queryDirectory' will be considered for analysis. ‘common_to_all’ finds sequence regions that are shared among all query strains but absent among all reference strains.

* 'muscleExecutable' specifies the full path to the muscle executable file.

* 'fragmentationSize' when running a ‘core’ analysis with 'coreInputType' set to 'panGenome', it determines the size of the fragments that the genomic sequences are segmented into.

* 'percentIdentityCutoff' when running a ‘core’ analysis sets the threshold of sequence identity for determining whether a fragment is part of the ‘core’ or ‘accessory’ genome.

* 'coreGenomeThreshold' defines the number of the initial sequences that a segment must be found in to be considered part of the 'core' genome; multi-fasta files of a single sequence are treated as a single sequence.

* 'mode' is currently required to be 'pan', for pan-genome analyses.