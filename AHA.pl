#!usr/bin/perl -w
use strict;

# Create an empty array to hold primers
my (@primers);
# Read in flags from the command line
my ($insSeqFile, $vectSeqFile, $printFile) = checkFlags (\@ARGV);

# Read in data from the insert sequence file [$insSeqFile] and the name of the vector sequence file [$vectSeqFile] and assign these to the
# sequence array [@sequences] and the vector file identifier [$vector]
my (@sequences) = readInsertSeqs ($insSeqFile);
my ($vector) = getVector($vectSeqFile);

# Checks to see how many inserts are intended to be inserted into the vector
if ($#sequences < 3) {
	##my (@primers) = makePrimers (\@sequences, $vector, 1);
	die ("Sorry, this program does not support less than two insert sequences at the moment. Check back later.\n");
} elsif ($#sequences == 3) {
	@primers = makePrimers (\@sequences, $vector, 2);
} else {
	die ("Sorry, this program cannot handle more than two insert sequences at the moment. Check back later.\n");
}

# Create variables for printing the DELPRIMERS format
my ($USName, $DSName, $VectName) = ($sequences[0], $sequences[2], $vector);
my (@dateTime) = localtime ();
$dateTime[5] += 1900;
my (@days) = ("SUN", "MON", "TUE", "WED", "THU", "FRI", "SAT");
my (@months) = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOv", "DEC");

# Print the primers to the designated output
open (DELPRIMERS, ">" . $printFile) || die ("Cannot open $printFile to write: $!\n");
write (DELPRIMERS);
close (DELPRIMERS) || die ("Cannot close $printFile: $!\n");
print ("Primer sequences were successfully generated and printed to $printFile\n");
#Subroutines ------->
sub checkFlags {
	my ($ARGVRef, $insFile, $vectFile, $printFile, $i, $response) = @_;
	if (defined ($ARGVRef)) {
		for ($i = 0; $i <= $#{$ARGVRef}; $i++ ) {
			if (${$ARGVRef}[$i] eq "-h" || ${$ARGVRef}[$i] eq "--help") {
				open (HELP, "<help.txt") || die ("Cannot open 'help.txt' for reading: $!\n");
				while (<HELP>) {
					print $_;
				}
				close (HELP) || die ("Cannot close 'help.txt' after reading: $!\n");
				exit (0);
			} elsif (${$ARGVRef}[$i] eq "-i" || ${$ARGVRef}[$i] eq "--insert" && ($i + 1) <= $#{$ARGVRef}) {
				chomp ($insFile = ${$ARGVRef}[$i + 1]);
			} elsif (${$ARGVRef}[$i] eq "-v" || ${$ARGVRef}[$i] eq "--vector" && ($i + 1) <= $#{$ARGVRef}) {
				chomp ($vectFile = ${$ARGVRef}[$i + 1]);
			} elsif (${$ARGVRef}[$i] eq "-p" || ${$ARGVRef}[$i] eq "--printFile" && ($i + 1) <= $#{$ARGVRef}) {
				chomp ($printFile = ${$ARGVRef}[$i + 1]);
				if (-e $printFile) {
					print ("$printFile already exists. Do you really wish to overwrite this file rather than append it? [y/N]: ");
					chomp ($response = <STDIN>);
					print ("$response\n");
					($response eq "y") || ($response eq "Y") || ($printFile = ">" . $printFile);
				}
			} elsif (${$ARGVRef}[$i] =~ /^-+.*/) {
				die ("${$ARGVRef}[$i] is an invalid flag option. Program terminated.\nFor more information use the -h/--help flag\n");
			}
		}
	}
	return ($insFile, $vectFile, $printFile);
}
sub readInsertSeqs {
	my ($fileName, $i) = @_;
	open (INSERTSEQS, $fileName) || die ("Cannot open $fileName for reading: $!\n");
	chomp (my (@lines) = <INSERTSEQS>);
	close (INSERTSEQS) || die ("Cannot close $fileName: $!\n");
	for ($i = 0; $i <= $#lines; $i+=2) {
		$lines[$i] =~ s/^>(.*)/$1/;
	}
	return (@lines);
}
sub getVector {
	my ($fileName) = @_;
	if ($fileName =~ /pjq200sk/i) {
		return ("pJQ200sk");
	} elsif ($fileName =~ /pbbpgdh/i) {
		die ("Sorry, this progrma cannot handle pBBPgdh primer sequence construction at the moment. Check back later.\n");
		#return ("pBBPgdh");
	} elsif (-e $fileName) {
		return ($fileName);
	} else {
		die ("$fileName does not exist. Program was terminated because no valid vector was loaded in.\n");
	}
}
sub revComp {
my ($seq, $revOnly, $compOnly) = @_;
	$seq =~ tr/a-z/A-Z/;
	if (not ($revOnly)) {
		$seq =~ tr/[ATCGN]/[tagcn]/;
		$seq =~ s/^(.*)/\U$1/;
	}
	if (not ($compOnly)) {
		$seq = reverse ($seq);
	}

	return ($seq);

}
sub findHomologies {
	my ($insertRef, $num, $vect, $const, $i, $seq, @tempSeq, @homologies) = @_;
	$const++;
	for ($i = 1; $i < ($num * 2); ($i += 2)) {
		$seq = ${$insertRef}[$i];
		{
			if (not ($vect)) {
				$seq =~ s/.*?\[(.+?)\](.{25})(.*)/$1,$2,$3/;
				@tempSeq = split (",", $seq);
				push (@homologies, $tempSeq[0] . "-" . $tempSeq[1]);
				if ($homologies[$i - $const] =~ /.{25,}.*-.*/) {
					$homologies[$i - $const] = $tempSeq[0];
				} else {
					$homologies[$i - $const] =~ s/(.{26}).*/$1/;
				}
				$seq = join ("", $tempSeq[1], $tempSeq[2]);
				if ($const == 0) {
					last;
				}
			}
			$seq =~ s/.*?(.{25})\[(.+?)\](.*)/$1,$2,$3/;
			@tempSeq = split (",", $seq);
			push (@homologies, $tempSeq[0] . "-" . $tempSeq[1]);
			if ($homologies[$i - $vect] =~ /.*-.{25,}/) {
				$homologies[$i - $vect] = $tempSeq[1];
			} else {
				$homologies[$i - $vect] =~ s/.*(.{26})/$1/;
			}
			if ($vect) {
				$seq = $tempSeq[2];
				$vect--;
				$const--;
				redo;
			}
		}

	}
	return (@homologies);
}
sub hom2Prim {
	my ($homologRef, $num, $i, $primer, @primers) = @_;
	if ($num == 1) {
		die ("Sorry, this program does not support single insert sequences at the moment. Check back later.\n")
	} elsif ($num == 2) {
		for ($i = 0; $i < (($num * 2) + 2); $i+=2) {
			$primer = ${$homologRef}[$i];
			$primer =~ /[-]/ && $primer =~ s/(.+)-.+/$1/;
			push (@primers, $primer);
			$primer = ${$homologRef}[$i + 1];
			$primer =~ /[-]/ && $primer =~ s/.+-(.+)/$1/;
			$primer = revComp($primer, 0, 0);
			push (@primers, $primer);
		}
		$primers[0] = ${$homologRef}[5] . "=" . $primers[0];
		$primers[1] = revComp (${$homologRef}[2], 0, 0) . "=" . $primers[1];
		$primers[2] = ${$homologRef}[1] . "=" . $primers[2];
		$primers[3] = revComp (${$homologRef}[4], 0, 0) . "=" . $primers[3];
		$primers[4] = ${$homologRef}[3] . "=" . $primers[4];
		$primers[5] = revComp (${$homologRef}[0], 0, 0) . "=" . $primers[5];
		for ($i = 0; $i <= $#primers; $i++) {
			$primers[$i] =~ tr/-//d;
			$primers[$i] =~ tr/=/-/;
			$primers[$i] =~ s/.*?(.{25}-)(.*)/$1$2/;
		}
	} else {
		die ("Sorry, this program cannot handle more than two insert sequences at the moment. Cehck back later.\n");
	}
	return (@primers);
}
sub makePrimers {
	my ($insertRef, $vector, $num, @homologies) = @_;
	@homologies = findHomologies ($insertRef, $num, 0);
	if ($vector eq "pJQ200sk") {
	#	push (@homologies, ("GCTCTAGAACTAGTGGATCCCCCGG","CCCTCACTAAAGGGAACAAAAGCTGGA"));
		push (@homologies, ("GCTCTAGAACTAGTGGATCCCCCGG","CCCTCACTAAAGGGAACAAAAGCTGGAG"));
#	} elsif ($vector eq "pBBPgdh") {
#		push (@homologies, ("GCTTGATATCGAATTCCTGCAGCCCG", "GCGCTCACTGGCCGTCGTTTTACAA"));
	} else {
		my (@vector) = readInsertSeqs ($vector);
		push (@homologies, reverse (findHomologies(\@vector, 1, 1)));
	}
	my (@primers) = hom2Prim (\@homologies, $num);
	return (@primers);
}

#Formats ->
format DELPRIMERS =
Created on @<< @<< @< @<<< at @<:@<:@<
$days[$dateTime[6]], $months[$dateTime[4]], $dateTime[3], $dateTime[5], $dateTime[2], $dateTime[1], $dateTime[0]
------------------------------------------------------------------------------------------
Upstream:   @<<<<<<<<<<<<<<<<
$USName
Downstream: @<<<<<<<<<<<<<<<<
$DSName
Vector:     @<<<<<<<<<<<<<<<<
$VectName


Primers:            (5')     Overhang homology    -    Primer                         (3')
==========================================================================================
Upstream Forward:   (5') @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< (3')
$primers[0]
Upstream Reverse:   (5') @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< (3')
$primers[1]
Downstream Forward: (5') @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< (3')
$primers[2]
Downstream Reverse: (5') @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< (3')
$primers[3]
Vector Forward:     (5') @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< (3')
$primers[4]
Vector Reverse:     (5') @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< (3')
$primers[5]
..........................................................................................

.
