## OVERVIEW

Panseq determines the core and accessory regions among a collection of genomic sequences based on user-defined parameters. It readily extracts regions unique to a genome or group of genomes, identifies SNPs within shared core genomic regions, constructs files for use in phylogeny programs based on both the presence/absence of accessory regions and SNPs within core regions.

It also provides a loci selector that efficiently computes the most discriminatory loci from a tab-delimited dataset.

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

The loci finder can be accessed from:

> lib/loci_selector.pl


## SETUP

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

## Testing your installation

	perl t/output.t

This will run a test suite against the included test data to ensure that Panseq is configured and working correctly. All tests should pass. The default location for the external programs used by the test suite are as follows:

	#program locations
	my $blastDirectory = '/usr/bin/';
	my $mummerDirectory = '/usr/bin/';
	my $muscleExecutable = '/usr/bin/muscle';

If your system setup is different, edit the t/output.t file, beginning at line 16 to change the directories where Blast+ and MUMmer reside, as well as the executable location for Muscle.

## Running Panseq

All the adjustments to Panseq are made by modifying a tab-delimited configuration file, which is specified as the only argument to the script.

	perl lib/panseq.pl settings.txt

Below is an example configuration file for panseq.pl:

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
	runMode 	pan


* ‘queryDirectory’ should contain the full directory path of the folder where all of the query sequences you are interested in comparing reside. Panseq will use the entire contents of this folder. 

* ‘referenceDirectory’ should contain the full directory path of the folder where all of the reference sequences you are interested in comparing reside. Panseq only uses this folder for Novel Region Comparisons of type 'common_to_all' and 'no_duplicates'. All other analyses use the contents of the Query folder only.

* ‘baseDirectory’ is the directory where all the output from Panseq is placed, and should be the full directory path. 

* 'numberOfCores' sets the number of processors available to Panseq. Increasing this can reduce run times.

* 'mummerDirectory' specifies the full path to the folder containing the nucmer program.

* 'blastDirectory' specifies the full path to the BLAST+ bin directory.

* 'minimumNovelRegionSize' sets the size in bp of the smallest region that will be kept by the Novel Region Finder; all regions found below this value will not be kept.

* 'novelRegionFinderMode'  sets the type of novel region analysis that will be performed. ‘no_duplicates’ finds the novel regions among one or more query strains with respect to the reference strains selected. ‘unique’ finds sequence regions that are unique to each sequence among all of the strains selected; only strains in the 'queryDirectory' will be considered for analysis. 

* 'muscleExecutable' specifies the full path to the muscle executable file.

* 'fragmentationSize' when running in mode `pan' determines the size of the fragments that the genomic sequences are segmented into.

* 'percentIdentityCutoff' when running in mode `pan' sets the threshold of sequence identity for determining whether a fragment is part of the ‘core’ or ‘accessory’ genome.

* 'coreGenomeThreshold' defines the number of the initial sequences that a segment must be found in to be considered part of the 'core' genome; multi-fasta files of a single sequence are treated as a single sequence.

* 'runMode' can be either 'novel' or 'pan', for novel-region finding and pan-genome analyses respectively.

##Detailed explanation of Panseq

###Novel Region Finder
 The Novel Region Finder currently has two modes implemented: "no_duplicates" and "unique". The no_duplicates mode identifies any genomic regions present in any of the query sequences that are not present in any of the reference sequences, and returns these regions in multi-fasta format. The "unique" mode finds genomic regions that are unique to each of the query sequences and not present in any of the reference sequences.

 These comparisons are done using the nucmer program from MUMmer 3, the parameters of which can be adjusted by the user prior to submitting an analysis. A sample of the output from novelRegions.fasta looks as follows:

	>lcl|strain1|contig000001_(505..54222)
	CCGTACGGGATTA...

Where the name of the contig containing the novel region is listed, followed by the nucleotide positions that were determined to be "novel" based on the comparison run. Lastly, the corresponding novel nucleotide sequence is included. 

###Core / Accessory Analysis

 All of the query strains are used to determine a non-redundant pan-genome. This is done by choosing a seed sequence for the pan-genome and iteratively building the pan-genome by comparing non-seed sequences to the "pan-genome", using the Novel Region Finder described above. For each comparison, sequences not present in the "pan-genome" are added, and the expanded "pan-genome" is used for the comparison against the next sequence. This iteration continues until a non-redundant pan-genome has been constructed.

Following the creation of this pan-genome for the selected sequences, the pan-genome is segmented into fragments of user-defined size. These fragments are subsequently queried against all of the sequences in the query list using blastn. The presence or absence of each pan-genome fragment is determined for each query, based on the Sequence Identity threshold set by the user. Pan-genome fragments present in a minimum number of genomes (determined by the core genome threshold) are aligned using Muscle.

Single nucleotide polymorphisms (SNPs) in these alignments are determined and used to generate a Phylip formatted file of all SNPs for use in downstream phylogenetic analyses. A phylip formatted file of pan-genome fragment presence / absence is also created, as are tab-delimited tables for both SNP and pan-genome fragment presence / absence (snp_table.txt and binary_table.txt, respectively), and a detailed result file listing the names, positions and values of the SNP and binary data (core_snps.txt and pan_genome.txt).

The detailed results file provides data in six columns: Locus Id, Locus Name, Genome, Allele, Start bp and Contig. Locus Id will be a 10-digit number that is a unique identifier for the locus; this number is used in the tabular output for both the SNP and binary data, and can be used for cross-referencing. Locus Name and Genome provide the human-readable names for the locus and the Genome. The Allele columns lists the actual data for the comparison. For pan-genome fragment presence / absence this is binary "0" or "1" data. For the SNP table, this is "A", "C", "T", "G" or "-". Start bp refers to the nucleotide position of the locus in base pairs, for example "45933"; the nucleotide position information is for the start of the fragment for the binary data. Lastly, the Contig column lists the name of the contig the locus is found on, which may differ than the Genome column, for genomes comprised of multi-fasta files. 

##Running the loci selector

The loci selector takes two command line arguments, with an optional third:

	perl loci_selector.pl input_file number_of_loci maximize_pod > output_file

The number_of_loci can be any positive integer, or 'best', which will find the minimum number of loci to provide a unique fingerprint for each data column, if possible. The default value of maximize_pod is 0, but can be set to 1 to disable masking of previously used locus pairs when calculating the points of discrimination. The default is recommended. Output is to STDOUT, and can be redirected to output_file.

## Detailed explanation of loci selector

The loci selector constructs loci sets that are maximized with respect to the unique number of fingerprints produced among the input sequences as well as the discriminatory power of the loci among the input sequences. The final loci set is iteratively built, in the following steps, given a tab-delimited table with loci names in the first column, sequence names in the first row, and single character data filling the matrix. Missing data is denoted by the characters '?', '-', or '.' :

* (1) Each potential available locus is evaluated for the number of unique fingerprints that would result from its addition to the final loci set. All loci that would generate the maximum number of unique fingerprints in this respect are evaluated in step (2).

* (2) All loci from step (1) are evaluated for their discriminatory power among the sequences, which is given as points of discrimination (POD). The POD for a locus is calculated as follows. A listing of all possible pair-wise comparisons is constructed; for example, if the input table consisted of three sequences, A, B and C, the list would consist of A-B, A-C and B-C. Next, it is determined whether or not the sequences in each pair-wise comparison contain the same single character denoting the locus state. If they do, a value of 0 is assigned; if they differ a value of 1 is assigned. The POD is then the summation of all pair-wise comparisons that differ for that locus. With our previous example, if A-B = 1, A-C = 1 and B-C = 0, the POD for that locus would be 2.

* (3) The locus with the highest value from step (2) is selected for addition to the final loci set and removed from the pool of candidate loci. If two or more loci tie in value, one is randomly selected. If all possible unique fingerprints have been found, the algorithm continues with (4); if additional unique fingerprints are possible, the algorithm continues with (5).

* (4) Sequence pairs for which the allele of the locus chosen in (3) differ are temporarily excluded from the analysis ("masked"). This ensures loci that differ between other pairs of strains are preferentially considered. Consider our A, B and C example with pair-wise comparisons of A-B = 1, A-C = 1 and B-C = 0. In the case of this locus being chosen, the sequence pairs A-B and A-C would be temporarily removed from the analysis ("masked"), leaving only loci that differed between B-C as viable options. Setting maximumPod = 1 prevents this masking step, which can be useful if one is only interested in loci that offer the most discrimination, regardless of what locus pairs offer that discrimination.

* (5) Once a locus has been chosen:

  * a) the specified number of loci has been reached (all unique fingerprints in the case of 'best') and the algorithm terminates; or

  * b) the specified number of loci has not been reached and there are remaining fingerprints possible, or sequence pairs for which differences exist. The algorithm returns to (1); or

  * c) there are no remaining fingerprints possible and no sequence pairs for which differences exist. At such time, all sequence pairs are again considered part of the analysis ("unmasked"). If no differences among any sequence pairs exist at this point, the algorithm terminates; if differences remain, the algorithm returns to (1). 
