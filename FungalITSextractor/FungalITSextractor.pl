#!/usr/bin/perl -w
use strict;

# Program name: ITS1/ITS2 extractor for fungal ITS sequences
# Program purpose: Extraction of ITS1/ITS2 from fungal ITS sequences in the FASTA format
# Copyright (C) 2009 R. Henrik Nilsson (henrik.nilsson@dpes.gu.se) et al.
# Henrik Nilsson, University of Gothenburg, Department of Plant and Environmental Sciences, Box 461, 405 30 Göteborg,  Sweden
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program (as gpl.txt); if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

my $infile="indata/indata.fasta";                                                                   # location of indata file

my $i=0; my $ITS2; my $ITS1; my $n58S;                                 # loop counter and variables to hold the sequence data
my $starttime=localtime;                                                                                     # get start time

my @Accno;                                                                     # where the query accession numbers are stored
my @Sequence;                                                                       # where the query sequence data is stored
my @date;                                                                              # where to store humanly readable time
my $outfile;                                                                                                # name of outfile

my $n58S_End;                                                    # to hold the absolute position of the  last (3') bp of 5.8S
my $n58S_Start;                                                  # to hold the absolute position of the first (5') bp of 5.8S
my $nLSU_Start;                                                  # to hold the absolute position of the first (5') bp of nLSU
my $nSSU_End;                                                    # to hold the absolute position of the  last (3') bp of nSSU 

my ($ITS1extracted, $ITS2extracted, $ITS1andITS2extracted);    # variables to keep track of what has been extracted for query
my $ITS1extractedcounter=0;                                          # to keep a global count of the number of extracted ITS1
my $ITS2extractedcounter=0;                                          # to keep a global count of the number of extracted ITS2
my $ITS1andITS2extractedcounter=0;                                                # global count: number of joint extractions
my $ITS1notITS2extractedcounter=0;                                                   # global count: ITS1 extracted, not ITS2
my $notITS1butITS2extractedcounter=0;                                                # global count: ITS2 extracted, not ITS1
my $noneextractedcounter=0;                                                            # global count: no subregion extracted
my $falsepositivescounter=0;

my $ITS1Starts=-1; my $ITS1Ends=-1; my $ITS2Starts=-1; my $ITS2Ends=-1; # for output to screen: absolute positions of regions
my $locallength=0;                                                                       # for computation of sequence length
my $warning;                                           # keep track of instances of palse positives; don't extract from these
my $warningreason;                                                                        # stores the reason for the problem

my $revcomp; my $iplusone;                            # flag for reverse complementary entries; humanly readable loop counter

&CheckForCriticalFiles;  # subroutine to check that everything looks OK prior to start (and to kill program if things=not OK)
&LoadIndata;                                                     # subroutine to load indata (in the FASTA format) into array

my $length = scalar @Accno -1;                                                  # $length holds the number of query sequences

my $sequencefile='UPsequence.txt';                                                  # name of file containing target sequence
my $number=0; my $numberofquerysequences = scalar @Accno;        # how many query sequences are there (not for loop purposes)
my $ITS2sequence; my $ITS1sequence;                                            # temporary storage of extracted sequence data
my $brief; my $temp;                                     # shorten any long accession numbers for structured output to screen

for ($i=0; $i<=$length ; $i++) {                                                              # loop over all query sequences

    $ITS1Starts=-1; $ITS1Ends=-1; $ITS2Starts=-1; $ITS2Ends=-1;   # set to -1 to handle cases where gene starts at position 0
    $locallength=length($Sequence[$i]);                                          # length of query sequence for screen output

    $number=$i;                                                        # use number to avoid handling the loop counter itself
    $brief=''; $n58S="no"; $temp="no";           # reset variables between each iteration; here brief name for query sequence
    $ITS1extracted="no";          # "no" signals that no ITS1 was extracted from the query sequence; reset between iterations
    $ITS2extracted="no";          # "no" signals that no ITS2 was extracted from the query sequence; reset between iterations
    $ITS1andITS2extracted="no";                 # "no" signals that the union was not extracted; reset between the iterations
    $warning="no";                            # set warning to no; it will be set to "yes" below if there are false positives
    $warningreason="none";                                                       # reset reason for any warning for each loop    

    $iplusone=$i+1;                                       # for humanly readable loop counter; not used for looping obviously
    if (length $Accno[$number]>20) {               # if the name of the accession name is longer than 20 chars, then truncate
	$brief=substr($Accno[$number],0,20);                  # get the first 20 chars of the query sequence by 3' truncation
	$brief = $brief . "..."; }                                # and append "..." to indicate that truncation has occurred
    else { 
	$brief = $Accno[$number];            # else no need for truncation; let's instead add spaces to create uniform output
	while (length $brief < 23) { $brief = $brief . " "; } }                            # adds the number of spaces needed
    
    $brief = $brief . " $locallength bp.";                              # append the length of the query to the output string
    
    print "($iplusone of $numberofquerysequences) $brief\t";                    # print progress report, counters, query name
    print OUTPUT "($iplusone of $numberofquerysequences) $brief\t";                                    # and echo to tsv file
    $n58S_End=$nLSU_Start="no";                                                  # reset trackkeepers between loop iterations
    $n58S_Start=$nSSU_End="no";                                                  # reset trackkeepers between loop iterations
    $revcomp="no";                                                               # reset trackkeepers between loop iterations
    $ITS2sequence="no";                                                          # reset trackkeepers between loop iterations

    open(FILE,">$sequencefile") or &err;                                                     # create FASTA file of the entry
    print FILE ">$Accno[$number]\n$Sequence[$number]\n";                      # print accno, sequence to temporary query file
    close(FILE) or &err;                                                                                  # close file or die

    $ITS2sequence=&GetITS2;                                      # call subroutine to try to extract ITS2 from query sequence
    $ITS1sequence=&GetITS1;                                      # call subroutine to try to extract ITS1 from query sequence


      # below: a large block to examine whether the internal order of nSSU, 5.8S, and nLSU makes sense. if not, don't extract
    if (($ITS1Ends > 0) and ($ITS2Ends > 0) and ($ITS1sequence ne "no") and ($ITS2sequence ne "no") and (($ITS1Ends > $ITS2Starts) or ("$ITS1Ends" eq "$ITS2Ends"))) {
	$ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no"; $warning="yes";         # blank all
	$falsepositivescounter++;                                                  # increase the counter for false positives
	if ($warningreason eq "none") { $warningreason="ITS1 and ITS2 overlap. Incorrect extraction."; } } 

    elsif (($ITS1Ends eq "-1") and ($ITS2Ends eq "-1") and ($ITS1sequence ne "no") and ($ITS2sequence ne "no")) {
	$ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no"; $warning="yes";         # blank all
	$falsepositivescounter++;                                                  # increase the counter for false positives
	if ($warningreason eq "none") { $warningreason="ITS1 and ITS2 overlap. Incorrect extraction."; } } 
    
    elsif (($nSSU_End=~/[0-9]+/) and ($n58S_Start=~/[0-9]+/) and ($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < 1000000) and ($n58S_Start < 1000000)) {
                        # checking of the order of nSSU, 5.8S, and nLSU is what we expect it to be; false positives otherwise
	unless ( ($nSSU_End < $n58S_Start) and ($n58S_Start < $n58S_End) and ($n58S_End < $nLSU_Start)) {
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no"; $warning="yes";     # blank all
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    if ($warningreason eq "none") { $warningreason="Order of nSSU, 5.8S, and nLSU incorrect."; } } }         # reason

    elsif (($n58S_Start=~/[0-9]+/) and ($n58S_End eq "1000000") and ($n58S_Start < 1000000) and ($nLSU_Start=~/[0-9]+/)) {
	$falsepositivescounter++;                                              # increase the counter for false positives
	$ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	if ($warningreason eq "none") { $warningreason="Start of 5.8S and LSU detected, but not end of 5.8S."; } }
    
    elsif (($nSSU_End=~/[0-9]+/) and ($n58S_Start=~/[0-9]+/) and ($n58S_End=~/[0-9]+/) and ($n58S_End < 1000000) and ($n58S_Start < 1000000)) {
	unless (($nSSU_End < $n58S_Start) and ($n58S_Start < $n58S_End)) {                   # unless SSU < 58Sstart < 58SEnd
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of nSSU and 5.8S incorrect."; } } }                # reason

    elsif (($n58S_Start=~/[0-9]+/) and ($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < 1000000) and ($n58S_Start < 1000000)) {
	unless (($n58S_Start < $n58S_End) and ($n58S_End < $nLSU_Start)) {                 # unless 5.8Start < 5.8SEnd < nLSU
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of 5.8S and nLSU incorrect."; } } }                # reason

    elsif (($nSSU_End=~/[0-9]+/) and ($n58S_Start=~/[0-9]+/) and ($n58S_Start < 1000000)) {
	unless ($nSSU_End < $n58S_Start) {                                                           # unless SSU < 5.8SStart
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of nSSU and 5.8S incorrect."; } } }                # reason

    elsif (($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < 1000000)) {
	unless ($n58S_End < $nLSU_Start) {                                                            # unless 5.8SEnd < nLSU
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of 5.8S and nLSU incorrect."; } } }                # reason
    
    elsif (($nSSU_End=~/[0-9]+/) and ($n58S_End=~/[0-9]+/) and ($n58S_End < 1000000) ) {
	unless ( $nSSU_End < $n58S_End)  {                                            # unless nSSU < n58SEnd, though strange
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of nSSU and 5.8S incorrect."; } } }                # reason

    elsif (($nSSU_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/)) {
	unless ( ($nSSU_End < $nLSU_Start)) {                                            # unless nSSU < nLSU, though strange
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of nSSU and nLSU incorrect."; } } }                # reason

    elsif (($n58S_Start=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_Start < 1000000)) {
	unless ( $n58S_Start < $nLSU_Start) {                                        # unless 58Sstart < nLSU, though strange
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of 5.8S and nLSU incorrect."; } } }                # reason

    elsif (($n58S_Start=~/[0-9]+/) and ($n58S_End=~/[0-9]+/) and ($n58S_End < 1000000) and ($n58S_Start < 1000000)) {
	unless ($n58S_Start < $n58S_End) {                                            # unless internal order of 5.8S correct
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of 5' and 3' of 5.8S incorrect."; } } }            # reason

    elsif (($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < 1000000)) {
	unless ( ($n58S_End < $nLSU_Start)) {                                                         # unless 5.8Send < nLSU
	    $falsepositivescounter++;                                              # increase the counter for false positives
	    $ITS1sequence="no"; $ITS2sequence="no"; $ITS1extracted="no"; $ITS2extracted="no";  $warning="yes";    # blank all
	    if ($warningreason eq "none") { $warningreason="Order of 5.8S and nLSU incorrect."; } } }

    else { ; }
                                            # end of large block with tests of internal consistency for the extracted regions

    unless ($warning eq "no") {
	print "ERROR: possible false positive. Will not extract. Reason: $warningreason";
	print OUTPUT "ERROR: possible false positive. Will not extract. Reason: $warningreason"; }

    else {

	unless (($ITS1sequence eq "no") or ($ITS2sequence eq "no")) {              # check that both sequences were extracted
                                              # if both end and start of 5.8S were detected, then we can extract 5.8S in full
	    if (($n58S_End=~/[0-9]+/) and ($n58S_End < 1000000) and ($n58S_Start=~/[0-9]+/) and ($n58S_Start < 1000000)) {
		$n58S=substr($Sequence[$i],$n58S_Start-1, ($n58S_End - $n58S_Start +1)); } }          # store 5.8 in variable
	
	if ($ITS1sequence eq "no") {                                                               # if no ITS1 was extracted
	    print OUTPUT "----\t";                                                                     # print this to screen
	    print "----\t"; }                                                                               # and to tsv file
	else {                                                             # else ITS1 was indeed extracted, whole or in part
	    print OUTPUT "ITS1\t";                                                                     # print this to screen
	    print "ITS1\t"; }                                                                               # and to tsv file
	
	if ($ITS2sequence eq "no") {                                                               # if no ITS2 was extracted
	    print OUTPUT "----\t";                                                                     # print this to screen
	    print "----\t"; }                                                                               # and to tsv file
	
	else {                                                             # else ITS2 was indeed extracted, whole or in part
	    print OUTPUT "ITS2\t";	                                                               # print this to screen
	    print "ITS2\t"; }                                                                               # and to tsv file
	
	print "\t\t";   # double tab to make clear distinction between extraction report and absolute positions of subregions
	print OUTPUT "\t\t";                                                                           # and echo to tsv file
	
	if (($ITS1Starts > 0) and ($ITS1Ends > 0) and ($ITS1sequence ne "no")) {               # we have a full ITS1 sequence
	    print OUTPUT "ITS1: $ITS1Starts-$ITS1Ends\t";    # print absolute positions of ITS1 in query sequence to tsv file
	    print "ITS1: $ITS1Starts-$ITS1Ends\t"; }           # print absolute positions of ITS1 in query sequence to screen
	elsif (($ITS1Starts>0) and ($ITS1sequence ne "no") ) {              # well at least we have a start position for ITS1
	    print OUTPUT "ITS1: $ITS1Starts-end\t";	     # print absolute positions of ITS1 in query sequence to tsv file
	    print "ITS1: $ITS1Starts-end\t"; }                 # print absolute positions of ITS1 in query sequence to screen
	elsif (($ITS1Ends>0) and ($ITS1sequence ne "no")) {                       # well at least we have an end for ITS1
	    print OUTPUT "ITS1: 1-$ITS1Ends\t";	             # print absolute positions of ITS1 in query sequence to tsv file
	    print "ITS1: 1-$ITS1Ends\t"; }                     # print absolute positions of ITS1 in query sequence to screen
	else {
	    print OUTPUT "-----      \t";                                          # we have no trace of ITS1; echo to screen
	    print "-----      \t"; }                                             # we have no trace of ITS1; echo to tsv file
	
	if (($ITS2Starts > 0) and ($ITS2Ends > 0) and ($ITS2sequence ne "no")) {                      # as above but for ITS2
	    print OUTPUT "ITS2: $ITS2Starts-$ITS2Ends\t";	
	    print "ITS2: $ITS2Starts-$ITS2Ends\t"; }
	elsif (($ITS2Starts>0) and ($ITS2sequence ne "no")) { 
	    print OUTPUT "ITS2: $ITS2Starts-end\t";
	    print "ITS2: $ITS2Starts-end\t"; }
	elsif (($ITS2Ends>0) and ($ITS2sequence ne "no")) {
	    print OUTPUT "ITS2: 1-$ITS2Ends\t";
	    print "ITS2: 1-$ITS2Ends\t"; }
	else {
	    print OUTPUT "----      \t";
	    print "----      \t"; }
	
	if ($revcomp eq "yes") {  # if "yes" then the query sequence is reverse complementary; indicate this in screen output
	    print OUTPUT "Reverse complementary.";                                                       # ...and to tsv file
	    print "Reverse complementary."; }
	else { ; } }                                                                                        # else do nothing

    print "\n";                        # print a lie break to screen to mark the end of the processing of that query sequence
    print OUTPUT "\n";               # print a lie break to tsv file to mark the end of the processing of that query sequence

    if (($n58S ne "no") and ($ITS1sequence ne "no") and ($ITS2sequence ne "no")) {        # if all three subregions extracted
	$temp="$ITS1sequence" . "$n58S" . "$ITS2sequence";                                       # concatenate the subregions
	print BothFASTA ">$Accno[$number]\n$temp\n"; }                                  # and save sequence to dedicated file

                                                         # below: to print FASTA directly to screen. commented out by default
#    print "GREE>$Accno[$i]full\nGREE$Sequence[$i]\n";
#    unless ($ITS1sequence eq "no") { print "GREE>$Accno[$i]ITS1\nGREE$ITS1sequence\n"; }
#    unless ($n58S eq "no") { print "GREE>$Accno[$i]58S\nGREE$n58S\n"; }
#    unless ($ITS2sequence eq "no") { print "GREE>$Accno[$i]ITS2\nGREE$ITS2sequence\n"; }
#    unless ($temp eq "no") { print "GREE>$Accno[$i]region\nGREE$temp\n"; }
#    print "GREE\nGREE\n";

    if ($ITS1extracted eq "yes") { print ITS1 ">$Accno[$number]\n$ITS1sequence\n"; }#print sccno to dedicated ITS1 FASTA file
    else { print NoITS1 ">$Accno[$number]\n$Sequence[$number]\n"; }                     # else print to file with ITS1 blanks

    if ($ITS2extracted eq "yes") { print ITS2 ">$Accno[$number]\n$ITS2sequence\n"; }#print accno to dedicated ITS2 FASTA file
    else { print NoITS2 ">$Accno[$number]\n$Sequence[$number]\n"; }                     # else print to file with ITS2 blanks

    if (($ITS1extracted eq "yes") and ($ITS2extracted eq "yes")) {    # print accno to file with joint successful extractions
	print Both "$Accno[$number]\n"; }
    elsif (($ITS1extracted eq "no") and ($ITS2extracted eq "no")) { print None "$Accno[$number]\n"; } #file with joint blanks

    else { ; }                                                                                              # else do nothing

    if (($ITS1extracted eq "yes") and ($ITS2extracted eq "yes")) { $ITS1andITS2extractedcounter++; }      # keep track of the
    elsif (($ITS1extracted eq "no") and ($ITS2extracted eq "no")) { $noneextractedcounter++; }       # success of the various
    elsif ($ITS1extracted eq "yes") { $ITS1notITS2extractedcounter++; }                # possible combinations of extractions
    elsif ($ITS2extracted eq "yes") { $notITS1butITS2extractedcounter++; }                               # for summary output
    else { print "WARNING: something is terribly wrong.\n"; }                                           # what just happened?

    $ITS1Starts=0; $ITS1Ends=0; $ITS2Starts=0; $ITS2Ends=0; }                            # reset variables between iterations   

close(ITS1) or &err;                                                                   # below: close the output files or die
close(ITS2) or &err;
close(NoITS1) or &err;
close(NoITS2) or &err;
close(None) or &err;
close(Both) or &err;
close(BothFASTA) or &err;

print "\n------------------------------------------------------------\n";           # now print summary statistics to screen
print "Total number of sequences in input file: $numberofquerysequences\n\n";             # these should be self-explanatory
print "ITS1 and ITS2 were extracted from $ITS1andITS2extractedcounter entries.\n";
print "Only ITS1 was extracted from $ITS1notITS2extractedcounter entries.\n";
print "Only ITS2 was extracted from $notITS1butITS2extractedcounter entries.\n";
print "Neither ITS1 nor ITS2 were extracted from $noneextractedcounter entries. \n";
print " ...including $falsepositivescounter cases of false positives\n";

unless ($ITS1andITS2extractedcounter + $ITS1notITS2extractedcounter + $notITS1butITS2extractedcounter + $noneextractedcounter == $numberofquerysequences) {
        print "WARNING: These do however not sum up to $numberofquerysequences. Something must have gone wrong.\n"; }

print "\nEnd of output.\n";
if (-e "UPsequence.txt") { unlink "UPsequence.txt"; }                       # delete these temporary work files from main dir
if (-e "UPoutput.txt") { unlink "UPoutput.txt"; }

my $endtime=localtime;                                               # get the time again to be able to show the time elapsed

print "Run started: $starttime \n";                                                                   # show the time elapsed
print "Run ended: $endtime \n";


exit;
#############################################################################################################################
################################################### end of main #############################################################
#############################################################################################################################




#############################################################################################################################
sub err {
    print "The script terminated with an error.\n";
    print "Exiting.\n";
    exit; }
#############################################################################################################################




#############################################################################################################################
sub LoadIndata {    # subroutine to read infile to memory. note: this subroutine reads all of the infile directly into memory
                            # with less than 100,000 query sequences, you should be fine even if you run on obsolete hardware
                                                      # 100,000 ITS sequences would take perhaps 80 MB of RAM [plus Perl tax]

    my $accno_set="no";                                   # flag to control whether an accession number has been found or not
    open(INFILE,$infile) or &err;                                                                   # open INFILE for reading
    my @infile=<INFILE>;
    close(INFILE) or &err;

    my $counter=0; my $Accnocounter=0; my ($sequence, $accno);                          # various parameters for internal use
    my $infilelength = scalar @infile -1;                                                     # how many lines in the infile?
    unless ($infile[$infilelength]=~/\n$/) {  $infile[$infilelength]=$infile[$infilelength] . "\n"; }# newline to end of file
    $infile[$infilelength+1]="\n"; $infilelength++;

    for ($counter=0; $counter<=$infilelength; $counter++) {                               # loop over the lines of the infile

	if (!$infile[$counter]) { next; }                                                                  # skip blank lines

	if (($accno_set eq "no") and ($infile[$counter] =~ /^>/)) {   # no accno found, but > makes this line looks promising
	    $accno=$infile[$counter];                                              # let the accession number be this line...
	    $accno=~s/\s//g; $accno=~s/>//g; $accno=~s/-//g;                       # ...after removal of whitespace, >, and -
	    $accno=substr($accno,0,30); # take only the first 20 characters; entry must be unique for the 20 first characters
	    $accno_set="yes"; next;  }                           # mark accession number as set; now go for the sequence data

	elsif (($accno_set eq "yes") and ($counter==$infilelength)) {                       # we have reached the final entry

	    $sequence=~s/\s//g;                                                        # remove whitespace from sequence data
	    $sequence=uc $sequence;                                                                         # make upper-case
	    $sequence=~tr/ATUGCYRSWKMBDHVN//cd;                            # let sequence data contain only C, T, A, G, and N
	    if ((!$sequence) or (length $sequence < 20)) { print "ERROR: Sequence of <20 bp. encountered: $accno ($sequence). Will not proceed.\n"; &err; }
	    $Accno[$Accnocounter]=$accno;                                   # save accession number to accession number array
	    $Sequence[$Accnocounter]=$sequence;                                             # save sequence to sequence array
	    $Accnocounter++;                                              # increase counter of saved query accession numbers
	    if ($counter==$infilelength) { last; } }                                             # we're done if this is true

	elsif (($accno_set eq "yes") and ($infile[$counter+1] =~ /^>/)) {              # if next line is a new entry, then...

	    if (!$sequence) { $sequence=" "; }     # one-line sequence; add temporary space to avoid accessing empty variable
	    $sequence = $sequence . $infile[$counter];        # add present line (which is the last of the entry) to sequence
	    $sequence=~s/\s//g;                                                        # remove whitespace from sequence data
	    $sequence=uc $sequence;                                                                         # make upper-case
	    $sequence=~tr/ATUGCYRSWKMBDHVN//cd;                            # let sequence data contain only C, T, A, G, and N
	    if ((!$sequence) or (length $sequence < 20)) { print "ERROR: Sequence of <20 bp. encountered: $accno ($sequence). Will not proceed.\n"; &err; }
	    $Accno[$Accnocounter]=$accno;                                   # save accession number to accession number array
	    $Sequence[$Accnocounter]=$sequence;                                             # save sequence to sequence array
	    $Accnocounter++;                                              # increase counter of saved query accession numbers
	    $sequence=" "; $accno_set="no"; }                                       # reset sequence parameter, accno tracker

	else { $sequence = $sequence . $infile[$counter]; } }                                   # append to sequence variable
    print "$Accnocounter sequences found.\n\n"; }                                                                 # that's it
#############################################################################################################################




#############################################################################################################################
sub GetITS2 {                                                        # subroutine to extract the ITS2 from the query sequence

    $n58S_End=&Findn58SEnd;                                  # subroutine to find the end of 5.8S uing pre-made HMMER profile
    unless ($n58S_End=~/[0-9]+/) { $n58S_End = &Findn58SEndTwo; }                          # else try again using shorter HMM

    unless ($n58S_End=~/[0-9]+/) { $n58S_End=1000000; }   # if not found, set to very high to be able to tell below (revcomp)

    $nLSU_Start=&FindnLSUStart;                 # subroutine to find the start of nLSU using long HMMER profile for precision
    unless ($nLSU_Start=~/[0-9]+/) { $nLSU_Start=&FindnLSUStartTwo; }            # else, try the shorter HMM (less precision)

    if (($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < $nLSU_Start)) {  # end of 5.8, start of nLSU found
	$ITS2=substr($Sequence[$i],$n58S_End,(($nLSU_Start-$n58S_End)-1));                    # subtract to get the ITS2 data
	if (!$ITS2) {                                   # found anchors at the very extreme of sequence, so no ITS2 extracted
	    return "no"; }                                                               # and don't return any sequence data
	$ITS2extracted="yes";                                       # flag to indicate that we have indeed extracted the ITS2
	return $ITS2; }                                                                           # return ITS2 sequence data

    elsif (($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < 1000000) and ($n58S_End > $nLSU_Start)) {
	$warning="yes";                    # we have a false positive match to either the end of 5.8S or to the start of nLSU
	$warningreason="End of 5.8S detected as being after the start of nLSU.";           # store the reason for the warning
	return "no"; }                        # so we can assume that we will not be able to extract any meaningful ITS2 data

    elsif (($n58S_End=~/[0-9]+/) and ($n58S_End < 1000000 )) {        # if end of 5.8S really set but not start of nLSU found
	$ITS2=substr($Sequence[$i],$n58S_End);                                     # let ITS2 be everything after end of 5.8S
	if (!$ITS2) {                                   # found anchors at the very extreme of sequence, so no ITS2 extracted
	    return "no"; }                                                               # and don't return any sequence data
	$ITS2extracted="yes";                                                             # flag as extracted, though partial
	return $ITS2; }                                                                           # return ITS2 sequence data                      

    elsif ($nLSU_Start=~/[0-9]+/) {                                                       # found nLSU start, not end of 5.8S
	$ITS2=substr($Sequence[$i],0,$nLSU_Start-1);                                             # extract ITS2 as start-nLSU
	if (!$ITS2) {                                   # found anchors at the very extreme of sequence, so no ITS2 extracted
	    return "no"; }                                                               # and don't return any sequence data
	$ITS2extracted="yes";                                                                    # flag successful extraction
	return $ITS2; }

    else {                                                        # neither 5.8S nor nLSU found. trying reverse complementary
	$Sequence[$i]=reverse $Sequence[$i];                                       # reverse sequence (but not complementary)
	$Sequence[$i]=uc($Sequence[$i]);                                                                      # go upper-case
	$Sequence[$i]=~tr/ATUGCYRSWKMBDHVN/TAACGRYSWMKVHDBN/d;                  # this tr handles the complementary operation

	open(FILE,">UPsequence.txt") or &err;          # open sequence file to write the reverse complementary sequence to it
	print FILE ">test\n";
	print FILE "$Sequence[$i]\n";	
	close (FILE) or &err;

	$n58S_End=&Findn58SEnd;                                                 # another attempt at locating the end of 5.8S
	unless ($n58S_End=~/[0-9]+/) { $n58S_End = &Findn58SEndTwo; }         # try again, using shorter HMM (less precision)

	unless ($n58S_End=~/[0-9]+/) { $n58S_End=1000000; }                                               # OK, found nothing

	$nLSU_Start=&FindnLSUStart;                                        # another attempt at locating the start of of nLSU
	unless ($nLSU_Start=~/[0-9]+/) { $nLSU_Start=&FindnLSUStartTwo; }                 # make less precise attempt an nLSU

	if (($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < $nLSU_Start)) { # end 5.8, start of nLSU found
	    $ITS2=substr($Sequence[$i],$n58S_End,(($nLSU_Start-$n58S_End)-1));              # extract accordingly to get ITS2
	    $revcomp="yes";               # since we managed to extract, the sequence is reverse complementary; flag for this
	    if (!$ITS2) {                               # found anchors at the very extreme of sequence, so no ITS2 extracted
		return "no"; }                                                           # and don't return any sequence data
	    $ITS2extracted="yes";                                                            # flag for successful extraction
	    return $ITS2; }                                                                         # and return it to caller

	elsif (($n58S_End=~/[0-9]+/) and ($nLSU_Start=~/[0-9]+/) and ($n58S_End < 1000000) and ($n58S_End > $nLSU_Start)) {
                 # sequence not reverse complementary after all; so reverse it back again an save to file for ITS1 extraction
	    $Sequence[$i]=reverse $Sequence[$i];                                   # reverse sequence (but not complementary)
	    $Sequence[$i]=uc($Sequence[$i]);                                                                     # upper-case
	    $Sequence[$i]=~tr/ATUGCYRSWKMBDHVN/TAACGRYSWMKVHDBN/d;                                       # turn complementary

	    open(FILE,">UPsequence.txt") or &err;      # open sequence file to write the reverse complementary sequence to it
	    print FILE ">test\n";
	    print FILE "$Sequence[$i]\n";	
	    close (FILE) or &err;
	    $warning="yes";                                                  # OK, we have a false positive. set warning flag
	    $warningreason="End of 5.8S detected as being after the start of nLSU.";       # store the reason for the warning
	    return "no"; }                                                                                # no, nothing found

	elsif (($n58S_End=~/[0-9]+/) and ($n58S_End < 1000000 )) {     # if end of 5.8S really found, not start of nLSU found
	    $ITS2=substr($Sequence[$i],$n58S_End);                                 # let ITS2 be everything after end of 5.8S
	    $revcomp="yes";   # flag that the entry really is reverse complementary and leave the revcomp file as-is for ITS1
	    if (!$ITS2) {                               # found anchors at the very extreme of sequence, so no ITS2 extracted
		return "no"; }                                                           # and don't return any sequence data
	    $ITS2extracted="yes";                                                              # flag that we indeed found it
	    return $ITS2;  }                                                                               # return to caller

	elsif ($nLSU_Start=~/[0-9]+/) {                                                   # found nLSU start, not end of 5.8S
		$ITS2=substr($Sequence[$i],0,$nLSU_Start-1);                       # extract sequence data from this position
		$revcomp="yes";                                    # signal that the sequence is indeed reverse complementary
		if (!$ITS2) {                           # found anchors at the very extreme of sequence, so no ITS2 extracted
		    return "no"; }                                                       # and don't return any sequence data
		$ITS2extracted="yes";                                      # flag that we indeed managed to extract something
		return $ITS2; }                                                                            # return to caller
	    
	else {                        # OK, so we found nothing. the sequence is probably not reverse complementary after all
	    $Sequence[$i]=reverse $Sequence[$i];                                   # reverse sequence (but not complementary)
	    $Sequence[$i]=uc($Sequence[$i]);                                                                # turn upper-case
	    $Sequence[$i]=~tr/ATUGCYRSWKMBDHVN/TAACGRYSWMKVHDBN/d;                                         # go complementary
	    
	    open(FILE,">UPsequence.txt") or &err;  # open sequence file to write the non-reverse complementary sequence to it
	    print FILE ">test\n";                                                           # for purposes of ITS1 extraction
	    print FILE "$Sequence[$i]\n";	
	    close (FILE) or &err;
	    
	    return "no"; } } }                                                                              # no trace at all
#############################################################################################################################




#############################################################################################################################
sub Findn58SEnd {                                                                      # subroutine to locate the end of 5.8S

    my @output=`hmmpfam -n -E 0.0015 HMMs/58endjames51.HMM UPsequence.txt`; #target file on HMMER profile; results in @output
    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $containsn58end="no";                                                 # default flag set to "no" = no end of 5.8 found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                         # does not contain end of 5.8S
	    else {$containsn58end="yes"; } } }                                    # else: does indeed match the HMMER profile

    if ($containsn58end eq "no") { return "no"; }                                    # return "no": does not contain 5.8S end

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS2Starts=$2+1;
    return $2; }                       # return the last bp of the HMMER 5.8S profile in the sequence = last position of 5.8S
#############################################################################################################################




#############################################################################################################################
sub Findn58SEndTwo {                                # subroutine to locate the end of 5.8S using shorter HMM (less precision)

    my @output=`hmmpfam -n -E 0.01 HMMs/58endjames23.HMM UPsequence.txt`;  # target file on HMMER profile; results in @output
    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $containsn58end="no";                                                 # default flag set to "no" = no end of 5.8 found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                         # does not contain end of 5.8S
	    else {$containsn58end="yes"; } } }                                    # else: does indeed match the HMMER profile

    if ($containsn58end eq "no") { return "no"; }                                    # return "no": does not contain 5.8S end

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS2Starts=$2+1;
    return $2; }                       # return the last bp of the HMMER 5.8S profile in the sequence = last position of 5.8S
#############################################################################################################################




#############################################################################################################################
sub FindnLSUStartTwo {                                         # subroutine to locate the start of nLSU, less precise version
    my @output=`hmmpfam -n -E 0.009 HMMs/nLSUjamesKH17.HMM UPsequence.txt`; #target file on HMMER profile; results in @output
    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $containsnLSUbeginning="no";                                       # default flag set to "no" = no start of nLSU found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                       # does not contain start of nLSU
	    else {$containsnLSUbeginning="yes"; } } }                             # else: does indeed match the HMMER profile

    if ($containsnLSUbeginning eq "no") { return "no"; }                           # return "no": does not contain nLSU start

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS2Ends=$1-1;
    return $1; }                              # return the first bp of the HMMER nLSU profile in the sequence = start of nLSU
#############################################################################################################################




#############################################################################################################################
sub FindnLSUStart {                                                                  # subroutine to locate the start of nLSU
    my @output=`hmmpfam -n -E 0.005 HMMs/nLSUjamesKH33.HMM UPsequence.txt`; #target file on HMMER profile; results in @output
    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $containsnLSUbeginning="no";                                       # default flag set to "no" = no start of nLSU found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                       # does not contain start of nLSU
	    else {$containsnLSUbeginning="yes"; } } }                             # else: does indeed match the HMMER profile

    if ($containsnLSUbeginning eq "no") { return "no"; }                           # return "no": does not contain nLSU start

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS2Ends=$1-1;
    return $1; }                                 # return the first bp of the HMMER nLSU profile in the sequence = nLSU start
#############################################################################################################################




#############################################################################################################################
sub CheckForCriticalFiles {             # subroutine to test if everything looks OK, which however is a difficult thing to do

    print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";   # guess what for

    print "Nilsson, Abarenkov et al. 2009. An open source software package for automated ITS1 and ITS2 extraction from fungal ITS sequences.\n";

    my @HMMtest=`hmmpfam -h`;                                   # attempt to run HMMER with help switch to see what comes out
    if (!$HMMtest[0]) {                                                                         # HMMER apparently not evoked
	print "\nFATAL ERROR: hmmpfam binary not found in path [with execution rights].\n";                     # not so good
	print "Please install HMMER / copy binary to location covered by the path of your shell.\n";   # so we can't continue
	exit; }
    unless ($HMMtest[0]=~/search one or more sequences against HMM database/) {                 # HMMER apparently not evoked
	print "\nFATAL ERROR: hmmpfam binary not found in path [with execution rights].\n";                  # can't continue
	print "Please install HMMER / copy binary to location covered by the path of your shell.\n";
	exit; }
    else { print " <OK> Found HMMER.\n"; }                                             # HMMER evoked by the command employed

    unless (-e $infile && -s $infile) {                                    # check that infile exists and is of non-zero size
	print "\nFATAL ERROR: indata file not found as 'indata.fasta' in directory 'indata'.\n";     # exit if unsatisfactory
	print "I lack data to work on. Terminating.\n";
	exit; }
    else { print "\nLocated infile as $infile.\n"; }                                              # else looks good; continue

     # the below lines to create a unique working directory for the BLAST run to save the results from ever being overwritten
    @date=gmtime;                                         # get date to create a unique directory in which to save alignments
    $date[99]="1";                                                         # counter to allow for >1 such directories per day
    $date[5]=$date[5]+1900;                                                              # make format of year comprehensible
    $date[4]++;                                                    # increase month by one to get comprehensible month number
    $date[100]=$date[5] . "-" . $date[4] . "-" . $date[3];                                      # construct name of directory
    $date[101]=$date[100] . "-" . $date[99];                                 # append the 1 as the first drectory of that day

    for (my $k=0; $k<=1000000; $k++) {                 # allow for one million directories per day. arbitrary limit of course
	if ((-e "outdata/$date[101]")) {                                             # if directory by this name exists, then
	    $date[99]++;                                                               # increase the trailing counter by one
	    $date[101]=$date[100] . "-" . $date[99]; }                                           # append new counter to name
	else {
	    mkdir "outdata/$date[101]" or &err;                                            # else create the unique directory
	    last; } }

    $outfile="outdata/$date[101]/outdata.csv";                    # this is the name, location of the comma separated outfile

    open(OUTPUT,">","outdata/$date[101]/outdata.csv") or &err;  # open file for writing. there are many of them but they will
    open(ITS1,">","outdata/$date[101]/ITS1.fasta") or &err;       # be needed to cover the needs of the mycological community
    open(NoITS1,">","outdata/$date[101]/NoITS1.fasta") or &err;
    open(ITS2,">","outdata/$date[101]/ITS2.fasta") or &err;
    open(NoITS2,">","outdata/$date[101]/NoITS2.fasta") or &err;
    open(None,">","outdata/$date[101]/None.txt") or &err;
    open(Both,">","outdata/$date[101]/Both.txt") or &err;
    open(BothFASTA,">","outdata/$date[101]/Both.fasta") or &err;

                                                                # below: print information on the file names with the results
    print "A tab-separated echo of the screen output will be saved as outdata/$date[101]/outdata.csv\n";
    print "All extracted ITS1 sequences will be saved as FASTA to outdata/$date[101]/ITS1.fasta\n";
    print "All entries for which ITS1 could not be detected will be saved as FASTA to outdata/$date[101]/NoITS1.fasta\n";
    print "All extracted ITS2 sequences will be saved as FASTA to outdata/$date[101]/ITS2.fasta\n";
    print "All entries for which ITS2 could not be detected will be saved as FASTA to outdata/$date[101]/NoITS2.fasta\n";
    print "Entries for which both of ITS1 and ITS2 were extracted will be summarized as outdata/$date[101]/Both.txt\n";
    print "...and a corresponding FASTA file (ITS1, 5.8S, ITS2) will be saved as outdata/$date[101]/Both.fasta\n";
    print "Entries for which neither ITS1 nor ITS2 were extracted will be summarized as outdata/$date[101]/None.txt\n";


    print "\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
    print "\nProcessing the query sequences... "; }                                                 # where the action starts
#############################################################################################################################




#############################################################################################################################
sub Findn58SStart {                                                   # subroutine to locate the start of 5.8S using long HMM
    my @output=`hmmpfam -n -E 1E-10 HMMs/58james1t54.HMM UPsequence.txt`;                                  # cf. HMMER syntax 
                                                                  # run target file on HMMER profile; save results in @output
    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $contains58beginning="no";                                          # default flag set to "no" = no start of 5.8 found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                                # does not contain 5.8S
	    else {$contains58beginning="yes"; } } }                      # else: OK, sequence was found to contain 5.8S start

    if ($contains58beginning eq "no") { return "no"; }                             # return "no": does not contain 5.8S start

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS1Ends=$1-1;
    return $1; }                              # return the first bp of the HMMER 5.8S profile in the sequence = start of 5.8S
#############################################################################################################################




#############################################################################################################################
sub Findn58SStartTwo {                                              # locate start of 5.8S using shorter HMM (less precision)
    my @output=`hmmpfam -n -E 0.008 HMMs/58james1t25.HMM UPsequence.txt`;                                  # cf. HMMER syntax 
                                                                  # run target file on HMMER profile; save results in @output
    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $contains58beginning="no";                                          # default flag set to "no" = no start of 5.8 found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                                # does not contain 5.8S
	    else {$contains58beginning="yes"; } } }                      # else: OK, sequence was found to contain 5.8S start

    if ($contains58beginning eq "no") { return "no"; }                             # return "no": does not contain 5.8S start

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS1Ends=$1-1;
    return $1; }                              # return the first bp of the HMMER 5.8S profile in the sequence = start of 5.8S
#############################################################################################################################




#############################################################################################################################
sub FindnSSUEndTwo {                               # make another go for the SSU, this time with a more relaxed (shorter) HMM
    my $targetfile;

    my @output=`hmmpfam -n -E 0.011 HMMs/nSSUtehler18bp.HMM UPsequence.txt`; # 0.034                       # cf. HMMER syntax
                                                                  # run target file on HMMER profile; save results in @output
    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $containsn58end="no";                                                # default flag set to "no" = no end of nSSU found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                         # does not contain end of nSSU
	    else {$containsn58end="yes"; } } }                                    # else: does indeed match the HMMER profile

    if ($containsn58end eq "no") { return "no"; }                                    # return "no": does not contain nSSU end

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS1Starts=$2+1;
    return $2; }                             # return the last bp of the HMMER nSSU profile in the sequence = last bp of nSSU
#############################################################################################################################




#############################################################################################################################
sub FindnSSUEnd {                  # find SSU end with (a long) 29 bp HMM (many SSU in GenBank are however shorter than that)
    my $targetfile;
    my @output=`hmmpfam -n -E 0.008 HMMs/nSSUtehler29bp.HMM UPsequence.txt`;                               # cf. HMMER syntax

    my $length = scalar @output-1;                                                     # how many lines of output did we get?
    my $i=0;
    my $containsn58end="no";                                                # default flag set to "no" = no end of nSSU found
    my (@temp, $temp, $temp1, $temp2);                                          # temp variables used in processing of output

    for ($i=0; $i<=$length; $i++) {                                                           # for each line of HMMER output
	unless ($output[$i]=~/^Scores for sequence family classification/) { next; }                     # move through array
	else {                                            # string found. now check if a match was found or not 3 lines ahead
	    if ($output[$i+3] =~ /no hits above thresholds/) { last; }                         # does not contain end of nSSU
	    else {$containsn58end="yes"; } } }                                    # else: does indeed match the HMMER profile

    if ($containsn58end eq "no") { return "no"; }                                    # return "no": does not contain nSSU end

    for ($i=0; $i<=$length; $i++) {                                                             # for all lines of the output

	unless ($output[$i]=~/^Alignments of top-scoring domains/) { next; }                              # move to this line
	else { $temp=$output[$i+1]; $temp2=$temp; last; } }                  # and the take a copy of it for processing below

    $temp2=~/( score [-]*\d+[.]*\d*)/; $temp2=$1; $temp2=~/([-]*\d+[.]*\d*)/; $temp2=$1;
    if ($temp2 < -6) { return "no"; }

    $temp=~ /(from \d+ to \d+)/;                                    # will match, e.g. "from 34 to 85" and store in buffer $1
    $temp=$1;                                                    # copy buffer $1 to $temp, overwriting old contents of $temp
    $temp=~ /(\d+) \D+ (\d+)/;                         # will match, e.g., "34 to 85". 34=start of HMMER 5.8S profile, 85=end

    $ITS1Starts=$2+1;
    return $2; }                             # return the last bp of the HMMER nSSU profile in the sequence = last bp of nSSU
#############################################################################################################################




#############################################################################################################################
sub GetITS1 {                                          # subroutine to oversee the extraction of ITS1 from the query sequence
    $n58S_Start=&Findn58SStart;                            # subroutine to find the start of 5.8S uing pre-made HMMER profile
    unless ($n58S_Start=~/[0-9]+/) { $n58S_Start = &Findn58SStartTwo; }               # else try shorter HMM (less precision)
    
    unless ($n58S_Start=~/[0-9]+/) { $n58S_Start = 1000000; }            # else indicate that the extraction was unsuccessful
    
    $nSSU_End=&FindnSSUEnd;                       # subroutine to find the end of nSSU using long HMMER profile for precision
    unless ($nSSU_End=~/[0-9]+/) { $nSSU_End=&FindnSSUEndTwo; }                # unless found SSU, try the less stringent HMM
    
                                                                                            # start of 5.8, end of nSSU found
    if (($n58S_Start=~/[0-9]+/) and ($nSSU_End=~/[0-9]+/) and ($n58S_Start > $nSSU_End) and ($n58S_Start < 1000000)) {
	$ITS1=substr($Sequence[$i],$nSSU_End,(($n58S_Start-$nSSU_End)-1));                    # subtract to get the ITS1 data
	if (!$ITS1) {                                   # found anchors at the very extreme of sequence, so no ITS1 extracted
	    return "no"; }                                                               # and don't return any sequence data
	$ITS1extracted="yes";                                               # flag that we indeed managed to extract the ITS1
	return $ITS1; }                                                                           # return ITS1 sequence data
    
                                # we have false positives. return "no" for unsuccessful extraction as a precautionary measure
    elsif (($n58S_Start=~/[0-9]+/) and ($nSSU_End=~/[0-9]+/) and ($n58S_Start < 1000000) and ($n58S_Start < $nSSU_End)) {
	$warning="yes";                                                      # OK, we have a false positive. set warning flag
	$warningreason="Start of 5.8S detected as being before the start of nSSU.";        # store the reason for the warning
	return "no"; }                                                     # nothing sensible found, so don't return anything
    
    elsif (($n58S_Start=~/[0-9]+/) and ($n58S_Start < 1000000)) {   # if start of 5.8S really extraced, not end of nSSU found
	$ITS1=substr($Sequence[$i],0,$n58S_Start-1);                           # let ITS1 be everything prior to start of 58S
	if (!$ITS1) {                                   # found anchors at the very extreme of sequence, so no ITS1 extracted
	    return "no"; }                                                               # and don't return any sequence data
	$ITS1extracted="yes";                                                            # flag successful extraction of ITS1
	return $ITS1; }                                                                           # return ITS1 sequence data                      
    
    elsif ($nSSU_End=~/[0-9]+/) {                                # found end of nSSU. assume ITS1 to be everything 3' of that
	$ITS1=substr($Sequence[$i],$nSSU_End);                                                                 # extract ITS1
	if (!$ITS1) {                                   # found anchors at the very extreme of sequence, so no ITS1 extracted
	    return "no"; }                                                               # and don't return any sequence data
	$ITS1extracted="yes";                                                         # flag extraction of ITS1 as successful
	return $ITS1; }                                                                                    # return to caller
    
    elsif ($revcomp eq "no") {                                    # neither 5.8S nor nSSU found. trying reverse complementary
	
	$Sequence[$i]=reverse $Sequence[$i];                                       # reverse sequence (but not complementary)
	$Sequence[$i]=uc($Sequence[$i]);                                                                         # upper-case
	$Sequence[$i]=~tr/ATUGCYRSWKMBDHVN/TAACGRYSWMKVHDBN/d;                  # this tr handles the complementary operation
	
	open(FILE,">UPsequence.txt") or &err;          # open sequence file to write the reverse complementary sequence to it
	print FILE ">test\n";                                                                 # else go reverse complementary
	print FILE "$Sequence[$i]\n";                                                   # print the sequence data to the file
	close (FILE) or &err;                                                              # exit if file could not be closed
	
	$n58S_Start=&Findn58SStart;                                           # another attempt at locating the start of 5.8S
	unless ($n58S_Start=~/[0-9]+/) { $n58S_Start = &Findn58SStartTwo; }                                 # try shorter HMM
	
	unless ($n58S_Start=~/[0-9]+/) { $n58S_Start = 1000000; }                                    # else flag as not found
	
	$nSSU_End=&FindnSSUEnd;                                              # another attempt at locating the end of of nSSU
	unless ($nSSU_End=~/[0-9]+/) { $nSSU_End=&FindnSSUEndTwo; }                       # make less precise attempt an nSSU

	                                                                              # start 5.8 realy found, end nSSU found
	if (($n58S_Start=~/[0-9]+/) and ($nSSU_End=~/[0-9]+/) and ($n58S_Start > $nSSU_End) and ($n58S_Start < 1000000)) {
	    $ITS1=substr($Sequence[$i],$nSSU_End,(($n58S_Start-$nSSU_End)-1));              # extract accordingly to get ITS1
	    $revcomp="yes";                                                             # flag entry as reverse complementary
	    if (!$ITS1) {                               # found anchors at the very extreme of sequence, so no ITS1 extracted
		return "no"; }                                                           # and don't return any sequence data
	    $ITS1extracted="yes";                                                    # flag for successful extraction of ITS1
	    return $ITS1; }                                                                         # and return it to caller
	
                                                                                  	 # 5.S start really found, looks good
	elsif (($n58S_Start=~/[0-9]+/) and ($nSSU_End=~/[0-9]+/) and ($n58S_Start < 1000000) and ($n58S_Start < $nSSU_End)) {
	    $Sequence[$i]=reverse $Sequence[$i];                                   # reverse sequence (but not complementary)
	    $Sequence[$i]=uc($Sequence[$i]);                                                                     # upper-case
	    $Sequence[$i]=~tr/ATUGCYRSWKMBDHVN/TAACGRYSWMKVHDBN/d;                                  # complementary operation
	    $warning="yes";                                                  # OK, we have a false positive. set warning flag
	    $warningreason="Start of 5.8S detected as being before the start of nSSU.";    # store the reason for the warning
	    return "no"; }                                                          # return no; nothing meaningful extracted
	
	elsif (($n58S_Start=~/[0-9]+/) and ($n58S_Start < 1000000)) {  # if start of 5.8S really found, not end of nSSU found
	    $ITS1=substr($Sequence[$i],0,$n58S_Start-1);                               # let ITS1 be everything prior to 5.8S
	    $revcomp="yes";                                                    # entry is reverse complementary; flag as such
	    if (!$ITS1) {                               # found anchors at the very extreme of sequence, so no ITS1 extracted
		return "no"; }                                                           # and don't return any sequence data
	    $ITS1extracted="yes";                                                 # flag that ITS1 was extracted successfully
	    return $ITS1; }                                                                                # return to caller
	
	elsif ($nSSU_End=~/[0-9]+/) {                                                                # only end of nSSU found
	    $ITS1=substr($Sequence[$i],$nSSU_End);          # which however is good enough for purposes of extraction of ITS1
	    $revcomp="yes";                                                         # flag the entry as reverse complementary
	    if (!$ITS1) {                               # found anchors at the very extreme of sequence, so no ITS1 extracted
		return "no"; }                                                           # and don't return any sequence data
	    $ITS1extracted="yes";                                                # flag for successful extraction of the ITS1
	    return $ITS1; }                                                    	                       # and return to caller
	else { 
	    $Sequence[$i]=reverse $Sequence[$i];                                   # reverse sequence (but not complementary)
	    $Sequence[$i]=uc($Sequence[$i]);                                                                  # go upper-case
	    $Sequence[$i]=~tr/ATUGCYRSWKMBDHVN/TAACGRYSWMKVHDBN/d;              # this tr handles the complementary operation
	    return "no"; } }                                                                          # no, nothing extracted
    
    else { return "no"; } }                                                                           # no, nothing extracted
#############################################################################################################################

