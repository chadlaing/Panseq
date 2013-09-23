#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use DBI;

my $baseDirectory = '/home/chad/panseq/output/genodo_sep6/';

_createTree($ARGV[0]);

sub _createTree{
	my $type = shift;

	#define SQLite db
	my $dbh = (DBI->connect("dbi:SQLite:dbname=" . $baseDirectory . "temp_sql.db","","")) or die("Could not connect to SQLite DB");
	
	my $sql = qq{SELECT strain.name,results.value, locus.id
		FROM results
		JOIN contig ON results.contig_id = contig.id
		JOIN strain ON contig.strain_id = strain.id
		JOIN locus ON results.locus_id = locus.id
		WHERE results.type = '$type'
		ORDER BY locus.id ASC
	};

	my $sth = $dbh->prepare($sql);
	my $didSucceed = $sth->execute();

	my $tableFH = IO::File->new('>' . $baseDirectory . $type . '_table.txt') or die("$!");
	my %results;
	my %loci;
	my $locus;
	my @genomeOrder;
	my $counter=0;
	while(my $row = $sth->fetchrow_arrayref){
		$counter++;
		if($counter % 100){
			
		}
		else{
			print '.';
		}	
		
	    if(defined $results{$row->[0]}){
	    	push @{$results{$row->[0]}},$row->[1];
	    }
	    else{
	    	$results{$row->[0]}=[$row->[1]];
	    }

	    if(defined $locus && ($locus ne $row->[2])){
	    	unless(defined $genomeOrder[0]){
	    		@genomeOrder = sort keys %loci;
	    		$tableFH->print("\t" . join("\t",@genomeOrder) . "\n");
	    	}
	    	$tableFH->print($locus);
	    	foreach my $genome(@genomeOrder){
	    		$tableFH->print("\t" . $loci{$genome});
	    	}
	    	$tableFH->print("\n");
	    }
	    $locus = $row->[2];
	    $loci{$row->[0]}=$row->[1];
	}
	
	if(defined $locus){
		$tableFH->print($locus);
		foreach my $genome(@genomeOrder){
			$tableFH->print("\t" . $loci{$genome});
		}
		$tableFH->print("\n");
	}
	
	$dbh->disconnect();
	$tableFH->close();	

	my $nameConversion = _printPhylipFile($type,\%results);


	if($type eq 'binary'){
		_printConversionInformation($nameConversion);
	}	
}


sub _printConversionInformation{
	my $hashRef =shift;

	my $conversionFH = IO::File->new('>' . $baseDirectory . 'phylip_name_conversion.txt') or die "$!";
	$conversionFH->print(
		'Number' . "\t" . 'Name' . "\n"
	);

	foreach my $number(sort keys %{$hashRef}){
		$conversionFH->print($number . "\t" . $hashRef->{$number} . "\n");
	}

	$conversionFH->close();	
}

sub _printPhylipFile{
	my $table = shift;
	my $results = shift;

	my $outFH = IO::File->new('>' . $baseDirectory . $table . '.phylip') or die "$!";

	my $counter=1;
	my %nameConversion;
	foreach my $genome(sort keys %{$results}){
		$nameConversion{$counter}=$genome;

		if($counter==1){
			$outFH->print(scalar(keys %{$results}) . "\t" . scalar(@{$results->{$genome}}) . "\n");
		}

		$outFH->print($counter . "\t" . join('',@{$results->{$genome}}) . "\n");
		$counter++;
	}
	$outFH->close();
	return \%nameConversion;
}


