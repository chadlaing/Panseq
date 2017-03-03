[![Master branch build status](https://travis-ci.org/chadlaing/Panseq.svg?branch=master "Master Build Status")](https://travis-ci.org/chadlaing/Panseq)

## OVERVIEW

_**Panseq**_ determines the core and accessory regions among a collection of genomic sequences based on user-defined parameters. It readily extracts regions unique to a genome or group of genomes, identifies SNPs within shared core genomic regions, constructs files for use in phylogeny programs based on both the presence/absence of accessory regions and SNPs within core regions.

It also provides a loci selector that efficiently computes the most discriminatory loci from a tab-delimited dataset.

If you find _**Panseq**_ useful please cite:

> Pan-genome sequence analysis using Panseq: an online tool for the rapid analysis of core and accessory genomic regions.
> Laing C, Buchanan C, Taboada EN, Zhang Y, Kropinski A, Villegas A, Thomas JE, Gannon VP.
> BMC Bioinformatics. 2010 Sep 15;11:461.

> Identification of Salmonella enterica species- and subgroup-specific genomic regions using Panseq 2.0.
> Laing C, Villegas A, Taboada EN, Kropinski A, Thomas JE, Gannon VP.
> Infect Genet Evol. 2011 Dec;11(8):2151-61.


## USAGE

The _**Panseq**_ standalone script can be accessed from:

> lib/panseq.pl

The loci finder can be accessed from:

> lib/loci_selector.pl


## SETUP

_**Panseq**_ requires Perl 5.10 or greater, and the following CPAN package to be installed:

    Module::Build
    
This package can be installed from the command-line using `cpan -i Module::Build`.
Following this, to install, and automatically retrieve the required CPAN packages for _**Panseq**_, do the following:

	perl Build.PL
	./Build installdeps


The following free, external programs must also be installed:

* [MUMmer:](http://mummer.sourceforge.net)
* [BLAST+:](http://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [Muscle:](http://www.drive5.com/muscle)

and optionally

* [cd-hit](http://weizhongli-lab.org/cd-hit/)

## Testing your installation

	perl t/output.t

This will run a test suite against the included test data to ensure that _**Panseq**_ is configured and working correctly. All tests should pass. The `cd-hit` tests will only be run if `cd-hit` is found. _**Panseq**_ checks for the installed programs on the local path, but they can be optionally specified as follows:

	perl t/output.t --blastDirectory '/home/user/local_blast/' --mummerDirectory '/home/user/local_mummer/' --cdhitDirectory '/home/user/cdhit' --muscleExecutable '/home/dir/muscle_executable'

## Running _**Panseq**_

All the adjustments to _**Panseq**_ are made by modifying a tab-delimited configuration file, which is specified as the only argument to the script.

	perl lib/panseq.pl settings.txt

Below is an example configuration file for `panseq.pl`:

	queryDirectory	/home/phac/panseq/queryLarge/
	referenceDirectory	/home/phac/panseq/referenceLarge/
	baseDirectory	/home/phac/panseq/output/panseq2test/
	numberOfCores	22
	mummerDirectory	/home/phac/MUMmer3.23/
	blastDirectory	/home/phac/ncbi-blast-2.2.29+/bin/
	minimumNovelRegionSize	500
	novelRegionFinderMode	no_duplicates
	muscleExecutable	/usr/bin/muscle3.8.31_i86linux64
	fragmentationSize	500
	percentIdentityCutoff	85
	coreGenomeThreshold	2
	runMode 	pan
	
Advanced Options

	queryFile	/home/phac/fileOfQuerySequence.fasta
	cdhitDirectory  /home/phac/cd-hit/
	storeAlleles	1
	allelesToKeep	2
	nameOrId	name
	frameshift	1
	overwrite	1
	maxNumberResultsInMemory 	500
	blastWordSize	11
	nucB    200
	nucC    65
	nucD    0.12
	nucG    90
	nucL    20
	cdhit   1

### Settings and their [DEFAULTS]

* `queryDirectory` [REQUIRED] should contain the full directory path of the folder where all of the query sequences you are interested in comparing reside. _**Panseq**_ will use the entire contents of this folder. 

* `baseDirectory` [REQUIRED] is the directory where all the output from _**Panseq**_ is placed, and should be the full directory path. 

* `runMode` [REQUIRED] can be either `novel` or `pan`, for novel-region finding and pan-genome analyses respectively.

* `referenceDirectory` [OPTIONAL] should contain the full directory path of the folder where all of the reference sequences you are interested in comparing reside. During the identification of novel regions step, these reference sequences will be screened out.

* `queryFile` [OPTIONAL] is an input file of fasta-formatted sequences. If the mode is set to `pan`, the sequences in the `queryFile` will be used instead of generation a pan-genome. Thus, the distribution and SNPs of the input sequences will be determined for all genomes in the queryDirectory. This can be useful, for example, for quickly generating a table of + / - values for a set of input genes against a set of genomes in the `queryDirectory`.

* `numberOfCores` [1] sets the number of processors available to _**Panseq**_. Increasing this can reduce run times.

* `mummerDirectory` [$PATH] specifies the full path to the folder containing the `nucmer` program.

* `blastDirectory` [$PATH] specifies the full path to the `blast+` bin directory.

* `cdhitDirectory` [$PATH] specifies the full path to the `cd-hit` program directory. 

* `muscleExecutable` [$PATH] specifies the full path to the `muscle` executable file.

* `minimumNovelRegionSize` [0] sets the size in bp of the smallest region that will be kept by the Novel Region Finder; all regions found below this value will not be kept.

* `fragmentationSize` [0] when running in mode `pan` determines the size of the fragments that the genomic sequences are segmented into. When set to `0`, no fragmentation of the input is done, which can be useful if specifying input via the `queryFile` option.

* `percentIdentityCutoff` [85] when running in mode `pan` sets the threshold of sequence identity for determining whether a fragment is part of the `core` or `accessory` genome.

* `coreGenomeThreshold` [3] defines the number of input sequences that a segment must be found in to be considered part of the `core` genome; multi-fasta files of a single genome are treated as a single sequence.


* `storeAlleles` [0] if set to 1, will store the allele matching the query sequence for each of the genomes and output them to `locus_alleles.txt`.

* `allelesToKeep` [1] if set, and if `storeAlleles` is set, determines the number of alleles per genome to keep, if multiple exist. They will be output to the `locus_alleles.txt` file, and every allele after the first will be appended with a `_a#` tag, where `#` is the allele number (eg. `_a2`).

* `nameOrId` [`id`] determines whether the individual locus ID string of numbers is output, or the name based on the query sequence in the files `binary_table.txt` and `snp_table.txt`.

* `frameshift` [0] includes frameshift only differences in SNP counts. Default behavior is to include only positions where there are also nucleotide differences. If gap-only differences are required, set this option to 1.

* `overwrite` [0] determines whether or not the specified `baseDirectory` will be overwritten if it already exists. This will cause all data in the existing directory to be lost. 

* `maxNumberResultsInMemory` [500] sets the number of pan-genome results to process before emptying the memory buffers and printing to file. Set this number higher if you want to limit the number of I/O operations. If you run into memory issues, lower this number.

* `blastWordSize` [20] sets the word size for the blastn portion of _**Panseq**_. For small values of `fragmentationSize` or `percentIdentityCutoff`, hits may be missed unless this value is lowered. (The default value for the `blastn` program is 11; _**Panseq**_ sets this to 20 as the default).

* `cdhit` [0] determines whether or not `cd-hit-est` is run on the pan-genome before identifying the distribution of the pan-genome (and SNPs among core regions) among the input sequences. Percent identity cutoff for `cd-hit-est` is taken from `percentIdentityCutoff`.

## Format of multi-fasta files ##

_**Panseq**_ currently only accepts fasta or multi-fasta formatted files. More than one genome may be in a single file, but for all genomes consisting of more than one contig, a distinct identifier must be present in the fasta header of each contig belonging to the same genome. For example, you have just assembled a new genome and are eager to analyze it. Your file consists of a number of contigs, similar to:

	>contig000001
	ACTGTTT...

	>contig000002
	CGGGATT...

The unique identifier could be the strain name or anything else of your choosing, but it must be included using the "local" designation: lcl|unique_identifer|. To reformat the above contigs, find and replace all ">" characters in your multi-fasta file with >lcl|unique_identifer|. Thus, if the unique identifier were "strain1", the reformatted contigs would look as follows:

	>lcl|strain1|contig000001
	ACTGTTT...

	>lcl|strain1|contig000002
	CGGGATT...

Common database file formats are supported by default, such as ref|, gb|, emb|, dbj|, and gi| and do not need to be modified as described above. For legacy purposes, the name=|unique_identifier| is supported in addition to lcl|unique_identifier|. Please note that spaces are not permitted in the unique identifier. Only letters (A-Z, a-z), digits (0-9) and the underscore "_" are valid characters. 

##Description of output files
- `accessoryGenomeFragments.fasta`: based on the run settings, all pan-genome fragments that are considered "accessory".
- `binary.phylip`: the presence / absence of the pan-genome among all genomes in the `queryDirectory` in phylip format.
- `binary_table.txt`: the presence / absence of the pan-genome among all genomes in the `queryDirectory` in tab-delimited table format.
- `core_snps.txt`: based on the run settings, a tab-delimited, detailed results file of all SNPs found. Includes genome name, contig name, nucleotide variant, and base-pair position.
- `coreGenomeFragments.fasta`: based on the run settings, all pan-genome fragments that are considered "core".
- `Master.log`: the log detailing program execution.
- `pan_genome.txt`: based on the run settings, a tab-delimited, detailed results file of all pan-genome regions. Includes genome name, contig name, presence / absence, and base-pair position for the pan-genome regions.
- `panGenome.fasta`: the non-fragmented pan-genome for the genomes in `queryDirectory`.
- `panGenomeFragments.fasta`: the fragmented pan-genome based on the `fragmentationSize` parameter.
- `phylip_name_conversion.txt`: the genomes in the phylip file are labeled as sequential numbers. This file maps the numbers back to the original names given in the input fasta files. Can be used by the `lib/treeNumberToName.pl` script to automatically convert a newick file labeled with numbers to the appropriate genome names.
- `snp.phylip`: a concatenated alignment of all SNPs found in the "core" genome regions for all genomes in the `queryDirectory`l
- `snp_table.txt`: the nucleotide values for all SNPs found in the "core" genome regions in tab-delimited table format.


##Detailed explanation of _**Panseq**_

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
