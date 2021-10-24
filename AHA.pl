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
	my ($ARGVRef, $i, $response, @flag, @files, %ARGVHash) = @_;


	# Check if a help flag was specified (matches -h and --help)
	if (grep /^-+h.*\b/, @{$ARGVRef}) {

		# Opens the file 'help.txt' for reading, prints the contents to STDOUT for the user, closes the file and exits the script
		open (HELP, "<help.txt") || die ("Cannot open 'help.txt' for reading: $!\n");
		while (<HELP>) {
			print $_;
		}
		close (HELP) || die ("Cannot close 'help.txt' after reading: $!\n");
		exit (0);
	}

	# Checks whether the number of flags and values are 6
	# This will catch if there is an unexpeected flag:value ratio or the wrong number of flags, and ensures the @{$ARGVRef} array can be mapped to a hash
	($#{$ARGVRef} + 1) == 6 || die ("Looks like you have the wrong number of flags, flag values, or both. Please check the help documentation with '-h' or '--help' and the README.md document.\n");

	# Convert the @ARGV values to a hash for simpler flag inquiries
	%ARGVHash = @{$ARGVRef};

	# Check for invalid flags since all valid flags should start with -{1,2}[ivp].
	# This should catch when wrong flags are called or if a flag has more than one value.
	grep (/^-{1,2}[^ivp\W]{1}.*\b/, keys (%ARGVHash)) && die ("You have used an invalid flag. Please check the help documentation with '-h' or '--help' and the README.md document.\n");

	# Loops through each expected flag and stores the associated values in the @flag array
	foreach $i (qw/i v p/) {
		@flag = grep (/^-+$i.*\b/, keys (%ARGVHash));
		push (@files, $ARGVHash{$flag[0]});
	}

	# Check if the file exists, and provide a warning and option to append rather than rewrite
	if (-e $files[2]) {
		print ("$files[2] already exists. Do you really wish to overwrite this file rather than append it? [y/N]: ");
		chomp ($response = <STDIN>);

		# prepend '>' to append the specified output file if the user requests
		($response =~ /^y.*\b/i) || ($files[2] = ">" . $files[2]);
	}

	# Return the user specified values for the insert sequence file, vector/vector sequence file, and the output file
	return (@files);
}

sub readInsertSeqs {
	my ($fileName, $i) = @_;

	# Open the user specified sequence file and read in all lines to the @lines array, removing newline characters, and then close the file
	open (INSERTSEQS, "<" .$fileName) || die ("Cannot open $fileName for reading: $!\n");
	chomp (my (@lines) = <INSERTSEQS>);
	close (INSERTSEQS) || die ("Cannot close $fileName: $!\n");

	# Removes the leading '>' from each sequence name line
	for ($i = 0; $i <= $#lines; $i+=2) {
		$lines[$i] =~ s/^>(.*)/$1/;
	}

	# Returns the @lines array containing the renamed sequences
	return (@lines);
}

sub getVector {
	my ($fileName) = @_;

	# Check if the user-specified vector matches a preloaded sequence and either return the sequence or a message indicating the status of the sequence then kill the program
	if ($fileName =~ /pjq200sk/i) {
		return ("pJQ200sk");
	} elsif ($fileName =~ /pbbpgdh/i) {
		die ("Sorry, this program cannot handle pBBPgdh primer sequence construction at the moment. Check back later.\n");
		#return ("pBBPgdh");
	} elsif (-e $fileName) {
		# Return the filename for a vector sequence if it exists
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
