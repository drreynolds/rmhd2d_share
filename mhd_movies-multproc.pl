#!/usr/bin/perl -w

###########################################################################
#  This script works with single or multi-processor output from the 2D 
#  resistive MHD simulations to do the following tasks:
#      1. If necessary, it first contatenates the output files from all 
#         processors into a single file for each time step.  It then 
#         deletes the individual processor output files.
#      2. It creates movies of the output as 3D surfaces over 2D contours
###########################################################################

# input the example numbers
my @exnumvec;
my $excount = 0;
if (scalar(@ARGV) > 0 ) {
    print ("\nCommand-line interface chosen.\n");
    print ("The following example numbers will be used to create movies:\n");
    my $i;
    for ($i=0;$i<scalar(@ARGV);$i++) {
	$exnumvec[$i] = $ARGV[$i];
	print ("   example $exnumvec[$i] \n");
	$excount++;
    }
}
else {
    print ("\nInput the example number to create movies from.\n");
    print ("Current results directory listing:\n\n");
    system '\ls --color results';
    print ("\nEnter each example number followed by the [enter] key.\n");
    print ("When finished, enter 0.\n");
    my $ext = 0;
    while ( $ext == 0 ) {
	my $choice;
	print "   choice: ";
	chomp($choice = <STDIN>);
	if ( int($choice) == $choice ) {
	    if ( $choice == 0 ) {
		$ext = 1;
		print "thank you\n"; }
	    else {
		$exnumvec[$excount] = $choice;
		$excount++;	}
	}
	else {
	    print "invalid choice $choice, choose again\n"; 
	}
    }
}

my $i;
for ($i=0;$i<$excount;$i++) {
  
    my $extest;
    if ($exnumvec[$i] < 10) { 
	$extest = "results/ex0$exnumvec[$i]/energyHistory.txt"; }
    else {
	$extest = "results/ex$exnumvec[$i]/energyHistory.txt";
    }
    
    if (-e $extest) {
	print "\nCreating movies for example $exnumvec[$i]";	
	exmovies($exnumvec[$i]);
    }
    else {
	print "\nERROR: example $exnumvec[$i] does not exist\n";
    }
}




########################################################
####    SUBROUTINES
########################################################


# exmovies(exnum)
#   subroutine to create movies for a particular example
sub exmovies {
    my $exnum = shift;

    # get the output directory
    my $exdir;
    if ($exnum < 10) {  # example number one digit long
	$exdir = "results/ex0$exnum/"; }
    else {              # example number two digits long
	$exdir = "results/ex$exnum/"; 
    }

    # force flush after every print statement
    $| = 1;
    
    # load in timelevel vector
    my @times = loadtimes($exdir.'energyHistory.txt');
    print "\n  There are $#times timelevels, ";

    # find dump interval
    my $ndump = get_ndump($exdir);
    print "with dump interval $ndump ";

    # find number of subdomains
    my $nsubd = get_nsubd($exdir);
    print "and $nsubd subdomains.\n";

    # if nsubd > 1, merge timelevel output files
    if ( $nsubd > 1 ) {
	merge_subd($ndump,$#times,$nsubd,$exdir);
    }

    # get min, max values for each variable
    my @mmvals = extreme_vals($exdir,$exnum,$ndump,@times);

    ###########  MULTI-PLOT MOVIE  ###########
    my $xlabel = "x";
    my $ylabel = "y";
    my $movname = $exdir."ex".$exnum."_allmovies.mpg";
    mult_movie( $xlabel, $ylabel, $mmvals[0], $mmvals[1], $mmvals[2], 
		$mmvals[3], $mmvals[4], $mmvals[5], $mmvals[6], 
		$mmvals[7], $mmvals[14], $mmvals[15], $mmvals[8], 
		$mmvals[9], $mmvals[10], $mmvals[11], $mmvals[12], 
		$mmvals[13], $mmvals[16], $mmvals[17], $mmvals[18], 
		$mmvals[19], $movname, $exdir, $ndump, @times );

    ###########  DENSITY MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   my $zlabel = "rho";
    my $fnumber = 3;
    my $title = "Density";
    $movname = $exdir."ex".$exnum."_density.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[0], $mmvals[1], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

    ###########  X-VELOCITY MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   $zlabel = "vx";
    $fnumber = 4;
    $title = "X-Velocity";
    $movname = $exdir."ex".$exnum."_x-velocity.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[2], $mmvals[3], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

    ###########  Y-VELOCITY MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   $zlabel = "vy";
    $fnumber = 5;
    $title = "Y-Velocity";
    $movname = $exdir."ex".$exnum."_y-velocity.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[4], $mmvals[5], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

#     ###########  Z-VELOCITY MOVIE  ###########
#     $xlabel = "x";   $ylabel = "y";   $zlabel = "vz";
#     $fnumber = 6;
#     $title = "Z-Velocity";
#     $movname = $exdir."ex".$exnum."_z-velocity.mpg";
#     single_movie( $xlabel, $ylabel, $zlabel, $mmvals[6], $mmvals[7], $title, 
#                   $fnumber, $movname, $exdir, $ndump, @times);

    ###########  X-MAGNETIC FIELD MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   $zlabel = "Bx";
    $fnumber = 7;
    $title = "X-Magnetic Field";
    $movname = $exdir."ex".$exnum."_x-mag_field.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[8], $mmvals[9], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

    ###########  Y-MAGNETIC FIELD MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   $zlabel = "By";
    $fnumber = 8;
    $title = "Y-Magnetic Field";
    $movname = $exdir."ex".$exnum."_y-mag_field.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[10], $mmvals[11], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

#     ###########  Z-MAGNETIC FIELD MOVIE  ###########
#     $xlabel = "x";   $ylabel = "y";   $zlabel = "Bz";
#     $fnumber = 9;
#     $title = "Z-Magnetic Field";
#     $movname = $exdir."ex".$exnum."_z-mag_field.mpg";
#     single_movie( $xlabel, $ylabel, $zlabel, $mmvals[12], $mmvals[13], $title, 
#                   $fnumber, $movname, $exdir, $ndump, @times);

    ###########  PRESSURE MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   $zlabel = "P";
    $fnumber = 10;
    $title = "Pressure";
    $movname = $exdir."ex".$exnum."_pressure.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[14], $mmvals[15], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

    ###########  DIVB MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   $zlabel = "divB";
    $fnumber = 11;
    $title = "divB";
    $movname = $exdir."ex".$exnum."_divB.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[16], $mmvals[17], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

    ###########  CURRENT MOVIE  ###########
    $xlabel = "x";   $ylabel = "y";   $zlabel = "J";
    $fnumber = 12;
    $title = "Current";
    $movname = $exdir."ex".$exnum."_current.mpg";
    single_movie( $xlabel, $ylabel, $zlabel, $mmvals[18], $mmvals[19], $title, 
		  $fnumber, $movname, $exdir, $ndump, @times);

    # finished
    print ("\n  Movies finished, written to $exdir.\n\n");
    return 0;
}


# @mmvals = extreme_vals(exdir,exnum,ndump,@times)
#   subroutine to create multi-plot movie of all fields for example exnum
sub extreme_vals {
    my $exdir = shift;
    my $exnum = shift;
    my $ndump = shift;
    my @times = @_;

    # set minimum and maximum deformations
    print ("\n  setting axis limits\n");
    my @minmaxvals = var_minmax($#times,$exdir,$ndump);
    my $rhomin = $minmaxvals[0];    my $rhomax = $minmaxvals[1];
    my $umin   = $minmaxvals[2];    my $umax   = $minmaxvals[3];
    my $vmin   = $minmaxvals[4];    my $vmax   = $minmaxvals[5];
    my $wmin   = $minmaxvals[6];    my $wmax   = $minmaxvals[7];
    my $Bxmin  = $minmaxvals[8];    my $Bxmax  = $minmaxvals[9];
    my $Bymin  = $minmaxvals[10];   my $Bymax  = $minmaxvals[11];
    my $Bzmin  = $minmaxvals[12];   my $Bzmax  = $minmaxvals[13];
    my $Pmin   = $minmaxvals[14];   my $Pmax   = $minmaxvals[15];
    my $dBmin  = $minmaxvals[16];   my $dBmax  = $minmaxvals[17];
    my $Jmin   = $minmaxvals[18];   my $Jmax   = $minmaxvals[19];
#     print "\n Var_minmax reports the following axis limits:\n";
#     print "   rho (min,max) = ($rhomin,$rhomax)\n";
#     print "     u (min,max) = ($umin,$umax)\n";
#     print "     v (min,max) = ($vmin,$vmax)\n";
#     print "    Bx (min,max) = ($Bxmin,$Bxmax)\n";
#     print "    By (min,max) = ($Bymin,$Bymax)\n";
#     print "     P (min,max) = ($Pmin,$Pmax)\n";
#     print "  divB (min,max) = ($dBmin,$dBmax)\n";
#     print "     J (min,max) = ($Jmin,$Jmax)\n";


    # ensure that min/max values are not the same
    if ( $rhomin == $rhomax ) {
	# density is always positive
	$rhomin = 0.999*$rhomin; 
	$rhomax = 1.001*$rhomax;
    }
    if ( $umin == $umax ) {
	# velocities span zero
	$umin = $umin - 0.00000000000001;
	$umax = $umax + 0.00000000000001;
    }
    if ( $vmin == $vmax ) {
	# velocities span zero
	$vmin = $vmin - 0.00000000000001;
	$vmax = $vmax + 0.00000000000001;
    }
    if ( $wmin == $wmax ) {
	# velocities apan zero
	$wmin = $wmin - 0.00000000000001;
	$wmax = $wmax + 0.00000000000001;
    }
    if ( $Pmin == $Pmax ) {
	# pressure is always positive
	$Pmin = 0.999*$Pmin; 
	$Pmax = 1.001*$Pmax;
    }
    if ( $Bxmin == $Bxmax ) {
	# magnetic fields span zero
	$Bxmin = $Bxmin - 0.00000000000001;
	$Bxmax = $Bxmax + 0.00000000000001;
    }
    if ( $Bymin == $Bymax ) {
	# magnetic fields span zero
	$Bymin = $Bymin - 0.00000000000001;
	$Bymax = $Bymax + 0.00000000000001;
    }
    if ( $Bzmin == $Bzmax ) {
	# magnetic fields span zero
	$Bzmin = $Bzmin - 0.00000000000001;
	$Bzmax = $Bzmax + 0.00000000000001;
    }
    if ( $dBmin == $dBmax ) {
	# divB always spans zero
	$dBmin = $dBmin - 0.00000000000001;
	$dBmax = $dBmax + 0.00000000000001;
    }
    if ( $Jmin == $Jmax ) {
	# current always spans zero
	$Jmin = $Jmin - 0.00000000000001;
	$Jmax = $Jmax + 0.00000000000001;
    }

    # expand min/max values
    $rhomin = $rhomin - 0.01*($rhomax-$rhomin);
    $rhomax = $rhomax + 0.01*($rhomax-$rhomin);
    $umin   = $umin   - 0.01*($umax-$umin);
    $umax   = $umax   + 0.01*($umax-$umin);
    $vmin   = $vmin   - 0.01*($vmax-$vmin);
    $vmax   = $vmax   + 0.01*($vmax-$vmin);
    $wmin   = $wmin   - 0.01*($wmax-$wmin);
    $wmax   = $wmax   + 0.01*($wmax-$wmin);
    $Pmin   = $Pmin   - 0.01*($Pmax-$Pmin);
    $Pmax   = $Pmax   + 0.01*($Pmax-$Pmin);
    $Bxmin  = $Bxmin  - 0.01*($Bxmax-$Bxmin);
    $Bxmax  = $Bxmax  + 0.01*($Bxmax-$Bxmin);
    $Bymin  = $Bymin  - 0.01*($Bymax-$Bymin);
    $Bymax  = $Bymax  + 0.01*($Bymax-$Bymin);
    $Bzmin  = $Bzmin  - 0.01*($Bzmax-$Bzmin);
    $Bzmax  = $Bzmax  + 0.01*($Bzmax-$Bzmin);
    $dBmin  = $dBmin  - 0.01*($dBmax-$dBmin);
    $dBmax  = $dBmax  + 0.01*($dBmax-$dBmin);
    $Jmin   = $Jmin   - 0.01*($Jmax-$Jmin);
    $Jmax   = $Jmax   + 0.01*($Jmax-$Jmin);

    # put those values back into the minmaxvals array
    $minmaxvals[0]  = $rhomin;   $minmaxvals[1]  = $rhomax;
    $minmaxvals[2]  = $umin;     $minmaxvals[3]  = $umax;
    $minmaxvals[4]  = $vmin;     $minmaxvals[5]  = $vmax;
    $minmaxvals[6]  = $wmin;     $minmaxvals[7]  = $wmax;
    $minmaxvals[8]  = $Bxmin;    $minmaxvals[9]  = $Bxmax;
    $minmaxvals[10] = $Bymin;    $minmaxvals[11] = $Bymax;
    $minmaxvals[12] = $Bzmin;    $minmaxvals[13] = $Bzmax;
    $minmaxvals[14] = $Pmin;     $minmaxvals[15] = $Pmax;
    $minmaxvals[16] = $dBmin;    $minmaxvals[17] = $dBmax;
    $minmaxvals[18] = $Jmin;     $minmaxvals[19] = $Jmax;
    print "\n Using the following axis limits:\n";
    print "   rho (min,max) = ($rhomin,$rhomax)\n";
    print "     u (min,max) = ($umin,$umax)\n";
    print "     v (min,max) = ($vmin,$vmax)\n";
#    print "     w (min,max) = ($wmin,$wmax)\n";
    print "    Bx (min,max) = ($Bxmin,$Bxmax)\n";
    print "    By (min,max) = ($Bymin,$Bymax)\n";
#    print "    Bz (min,max) = ($Bzmin,$Bzmax)\n";
    print "     P (min,max) = ($Pmin,$Pmax)\n";
    print "  divB (min,max) = ($dBmin,$dBmax)\n";
    print "     J (min,max) = ($Jmin,$Jmax)\n";
    
    return @minmaxvals;
}




# $min = min(x1,x2)
#   subroutine to find the minimum element of two elements
sub min {
    return ($_[0] > $_[1]) ? $_[1] : $_[0];
}


# $max = max(x1,x2)
#   subroutine to find the maximum element of two elements
sub max {
    return ($_[0] < $_[1]) ? $_[1] : $_[0];
}


# @outvec = minmax( \@vec ), where @outvec = (min,max) 
#   subroutine to find the min and max elements of an array 
#   (array passed by reference)
sub minmax {
    my $vec = shift;
    my @out;
    $out[0] = @$vec[0];  # $out[0] = min
    $out[1] = @$vec[0];  # $out[1] = max
    foreach (@$vec) {
	my $val = $_;
	$out[0] = ($out[0] < $val) ? $out[0] : $val;
	$out[1] = ($out[1] > $val) ? $out[1] : $val;
    }
    return @out;
}


# @res = vecmult(scalar, \@vec)
#   subroutine to multiply a vector by a scalar 
#   (vector passed by reference)
sub vecmult {
    my $scal = shift;
    my $vec  = shift;
    my @res;
    my $i = 0;
    foreach (@$vec) {
	$res[$i] = $_ * $scal;
	$i++;
    }
    return @res;
}


# @res = vecsub(vec1,vec2)
#   subtracts vec2 from vec1 and returns result
sub vecsub {
    my @vec = @_;
    my @res;
    my $half = int ($#vec + 1) / 2;
    my $i;
    for ($i=0; $i<$half; $i++) { 
	$res[$i] = $vec[$i] - $vec[$i+$half];
    }
    return @res;
}


# @res = vecadd(vec1,vec2)
#   subroutine to add two vectors
sub vecadd {
    my @vec = @_;
    my @res;
    my $half = int ($#vec + 1) / 2;
    my $i;
    for ($i=0; $i<$half; $i++) { 
	$res[$i] = $vec[$i] + $vec[$i+$half];
    }
    return @res;
}


# @vec = loadtimes(filename)
#   subroutine to read in timelevels from an MHD output file
sub loadtimes {
    my $filename = $_[0];
    my $j = 0;
    my @tmpvec;
    my @timevec;
    open(INFILE, $filename) or die 
	"\nloadtimes error: couldn't open $filename!\n";
    while (<INFILE>) {
	my $inrow = $_;
	chomp($inrow);
	@tmpvec = split(" ", $inrow);
	$timevec[$j] = $tmpvec[0];
	$j++;
    }
    close(INFILE);
    return @timevec;
}


# ndump = get_ndump(exdir)
#   subroutine to output dump frequency of example in exdir
sub get_ndump {
    my $exdir = $_[0];
    my $j = 0;
    my $ndump = 0;
    my $testfile;
    until ( $ndump > 0 ) {
	$j++;
	if ( $j < 10 ) {    $testfile  = $exdir."output.001.00000".$j; }
	elsif ( $j < 100 ) { $testfile  = $exdir."output.001.0000".$j; }
	elsif ( $j < 1000 ) { $testfile  = $exdir."output.001.000".$j; }
	elsif ( $j < 10000 ) { $testfile  = $exdir."output.001.00".$j; }
	elsif ( $j < 100000 ) { $testfile  = $exdir."output.001.0".$j; }
	else  {                  $testfile  = $exdir."output.001.".$j; }
	if ( -e $testfile ) {
	    $ndump = $j;
	    return $ndump;
	}
    }
}


# nsubd = get_nsubd(exdir)
#   subroutine to output number of subdomains of example in exdir
sub get_nsubd {
    my $exdir = $_[0];
    my $j = 0;
    my $nsubd = 0;
    my $testfile;
    until ( $nsubd > 0 ) {
	$j++;
	if ( $j < 10 ) {    $testfile  = $exdir."output.00".$j.".000000"; }
	elsif ( $j < 100 ) { $testfile  = $exdir."output.0".$j.".000000"; }
	else {                $testfile  = $exdir."output.".$j.".000000"; }
	if ( -e $testfile ) {
	}
	else {
	    $nsubd = $j-1;
	    return $nsubd;
	}
    }
}


# merge_subd(ndump,Nt,nsubd,exdir)
#   subroutine to merge multiple subdomain output into single files
sub merge_subd {
    my $ndump = $_[0];
    my $Nt    = $_[1];
    my $nsubd = $_[2];
    my $exdir = $_[3];
    my $tj = 0;

    print "  merging output files in $exdir. Timelevel: ";
    # iterate over time levels
    for ($tj=0; $tj*$ndump<=$Nt; $tj++) {
	my $outnum = $tj*$ndump;
	print "$outnum, ";

	# set timelevel number for output files
	my $fnum;
	if    ($outnum < 10) {$fnum = ".00000".$outnum;}
	elsif ($outnum < 100) {$fnum = ".0000".$outnum;}
	elsif ($outnum < 1000) {$fnum = ".000".$outnum;}
	elsif ($outnum < 10000) {$fnum = ".00".$outnum;}
	elsif ($outnum < 100000) {$fnum = ".0".$outnum;}
	else                      {$fnum = ".".$outnum;}

	# iterate over subdomain output files for this time step
	my $tmploc = 0;
	my @tmpvec; $tmpvec[0] = "";
	my $dj = 0;
	for ($dj=1; $dj<=$nsubd; $dj++) {
	    
	    # set subdomain number and input filename for this subdomain
	    my $infile;
	    if ( $dj < 10 ) {    $infile  = $exdir."output.00".$dj.$fnum; }
	    elsif ( $dj < 100 ) { $infile  = $exdir."output.0".$dj.$fnum; }
	    else {                 $infile  = $exdir."output.".$dj.$fnum; }

	    # open subdomain timelevel file
	    open(INFILE, $infile) or die
		"\n  merge_subd error: couldn't open $infile, quitting.\n";

	    # go through subdomain timelevel file
	    while (<INFILE>) {
		# skip comment lines
		next if /^#/;

		# split the row into its usable parts
		my $inrow = $_;  # inrow holds all fields at one point
		chomp($inrow);   # remove newline characters
		my @tmpvec2 = (split " ", $inrow, 5);  # extract data at point

		# on non-empty rows, append the data to tmpvec[tmploc]
		if ( $#tmpvec2 > 3 ) {
		    $tmpvec[$tmploc] = $tmpvec[$tmploc].$inrow."\n";
		}

		# on empty rows, update tmploc and append new element to tmpvec
		else {
		    if ( ($tmpvec[$tmploc] ne "") ) {
			$tmploc=$tmploc+1;
			push @tmpvec, "";
		    }
		}
	    }
	    close(INFILE);

	    # delete individual subdomain file now that we have its data
	    unlink $infile;
	}
	# get rid of the remaining empty entry at the end of @tmpvec
	pop(@tmpvec);


	# We now have all timelevel data stored in @tmpvec, with each array 
	# element consisting of all field data at a line of spatial points.
	# We must now sort these points to iterate smoothly over the entire 
	# domain

	#    Create indexing array, so that sort doesn't need to move
	#    around all of the output data through memory
	@tmpindx = linspace(0,$#tmpvec,$#tmpvec+1);

	#    Cache-based sort, where we do a first sweep over the data to
	#    cache the indx location with subkeys based on corresponding
	#    location in tmpvec
	keys my %cache1 = @tmpindx;
	keys my %cache2 = @tmpindx;
	($cache1{$_}, $cache2{$_}) = 
	    map { (split " ", $tmpvec[$_], 3)[1], 
		  (split " ", $tmpvec[$_], 2)[0] } $_ 
	    for @tmpindx;

	# Multi-key sort through cached data
	@tmpindx = sort {
	    # primary comparison sorts yvals into ascending numerical order
	    $cache1{$a} <=> $cache1{$b}
		||
            # if yvals are equal, secondary comparison sorts xvals
	    $cache2{$a} <=> $cache2{$b}
	} @tmpindx;


	# Now that tmpvec is sorted, output to the merged timelevel file, 
	my $outfile = ">".$exdir."output.001".$fnum;
	open(OUTFILE, $outfile) or die
	    "\n  merge_subd error: couldn't open file $outfile, quitting.\n";
	
	# output header information and first row to outfile
	print OUTFILE "# Merged Gnuplot output file for RMHD simulation:\n";
	print OUTFILE "#       x             y            rho            u";
	print OUTFILE "             v             w             Bx        ";
	print OUTFILE "    By            Bz         pressure     ";
	print OUTFILE "divergence     jcurrent      tenergy\n";
	print OUTFILE $tmpvec[0];
	
	# print rows to outfile, skipping a line when switching yvals
	my $jpt;
	for ($jpt=1; $jpt<=$#tmpvec; $jpt++) {
	    my @tmpvec2 = split(" ", $tmpvec[$tmpindx[$jpt-1]], 3);
	    my @tmpvec3 = split(" ", $tmpvec[$tmpindx[$jpt]],   3);
	    if ($tmpvec2[1] == $tmpvec3[1]) {
		print OUTFILE $tmpvec[$tmpindx[$jpt]];
	    }
	    else {
		print OUTFILE " \n";
		print OUTFILE $tmpvec[$tmpindx[$jpt]];
	    }
	}
	close(OUTFILE);
    }
    return 0;
}


# @minmaxvec = tlevel_minmax(fname)
#   subroutine to find min/max values for each variable in a timelevel state
sub tlevel_minmax {
    my $fname = $_[0];
    my $j = 0;
    my @mmvec;
    $mmvec[0]  = 1e200;   $mmvec[1]  = -1e200;  # 0,1 corresp. to rho
    $mmvec[2]  = 1e200;   $mmvec[3]  = -1e200;  # 2,3 corresp. to vx
    $mmvec[4]  = 1e200;   $mmvec[5]  = -1e200;  # 4,5 corresp. to vy
    $mmvec[6]  = 1e200;   $mmvec[7]  = -1e200;  # 6,7 corresp. to vz
    $mmvec[8]  = 1e200;   $mmvec[9]  = -1e200;  # 8,9 corresp. to Bx
    $mmvec[10] = 1e200;   $mmvec[11] = -1e200;  # 10,11 corresp. to By
    $mmvec[12] = 1e200;   $mmvec[13] = -1e200;  # 12,13 corresp. to Bz
    $mmvec[14] = 1e200;   $mmvec[15] = -1e200;  # 14,15 corresp. to p
    $mmvec[16] = 1e200;   $mmvec[17] = -1e200;  # 16,17 corresp. to divB
    $mmvec[18] = 1e200;   $mmvec[19] = -1e200;  # 18,19 corresp. to J
    my @tmpvec;
    open(INFILE, $fname) or die 
	"\n  tlevel_minmax error: couldn't open $fname, quitting.\n";
    my $vartmp;
    while (<INFILE>) {
	# skip lines that start with a '#'
	next if /^#/;

	# split the row into its usable parts
	my $inrow = $_;
	chomp($inrow);
	@tmpvec = split(" ", $inrow);

	# skip rows with less than 4 entries
	next if ($#tmpvec < 4);

	# record the min, max values for each variable
	$vartmp = 1.0*$tmpvec[2];
	$mmvec[0] = ($mmvec[0] > $vartmp) ? $vartmp : $mmvec[0];
	$mmvec[1] = ($mmvec[1] < $vartmp) ? $vartmp : $mmvec[1];
	$vartmp = 1.0*$tmpvec[3];
	$mmvec[2] = ($mmvec[2] > $vartmp) ? $vartmp : $mmvec[2];
	$mmvec[3] = ($mmvec[3] < $vartmp) ? $vartmp : $mmvec[3];
	$vartmp = 1.0*$tmpvec[4];
	$mmvec[4] = ($mmvec[4] > $vartmp) ? $vartmp : $mmvec[4];
	$mmvec[5] = ($mmvec[5] < $vartmp) ? $vartmp : $mmvec[5];
	$vartmp = 1.0*$tmpvec[5];
	$mmvec[6] = ($mmvec[6] > $vartmp) ? $vartmp : $mmvec[6];
	$mmvec[7] = ($mmvec[7] < $vartmp) ? $vartmp : $mmvec[7];
	$vartmp = 1.0*$tmpvec[6];
	$mmvec[8] = ($mmvec[8] > $vartmp) ? $vartmp : $mmvec[8];
	$mmvec[9] = ($mmvec[9] < $vartmp) ? $vartmp : $mmvec[9];
	$vartmp = 1.0*$tmpvec[7];
	$mmvec[10] = ($mmvec[10] > $vartmp) ? $vartmp : $mmvec[10];
	$mmvec[11] = ($mmvec[11] < $vartmp) ? $vartmp : $mmvec[11];
	$vartmp = 1.0*$tmpvec[8];
	$mmvec[12] = ($mmvec[12] > $vartmp) ? $vartmp : $mmvec[12];
	$mmvec[13] = ($mmvec[13] < $vartmp) ? $vartmp : $mmvec[13];
	$vartmp = 1.0*$tmpvec[9];
	$mmvec[14] = ($mmvec[14] > $vartmp) ? $vartmp : $mmvec[14];
	$mmvec[15] = ($mmvec[15] < $vartmp) ? $vartmp : $mmvec[15];
	$vartmp = 1.0*$tmpvec[10];
	$mmvec[16] = ($mmvec[16] > $vartmp) ? $vartmp : $mmvec[16];
	$mmvec[17] = ($mmvec[17] < $vartmp) ? $vartmp : $mmvec[17];
	$vartmp = 1.0*$tmpvec[11];
	$mmvec[18] = ($mmvec[18] > $vartmp) ? $vartmp : $mmvec[18];
	$mmvec[19] = ($mmvec[19] < $vartmp) ? $vartmp : $mmvec[19];
    }
    close(INFILE);
    return @mmvec;
}


# @minmaxvals = var_mimnax(Nt,exdir,ndump)
#   subroutine to read in matrix rows as strings in a vector
sub var_minmax {

    # calling and local variables
    my $Nt    = $_[0];
    my $exdir = $_[1];
    my $ndump = $_[2];
    my $j;

    # force flush after every print statement
    $| = 1;
    
    # start with the initial state
    my $fname = $exdir."output.001.000000";
    print "  examining output files in $exdir :  0";
    my @tmpminmax = tlevel_minmax($fname);
    my @minmaxvals;
    $minmaxvals[0]  = $tmpminmax[0];  
    $minmaxvals[1]  = $tmpminmax[1];
    $minmaxvals[2]  = $tmpminmax[2];  
    $minmaxvals[3]  = $tmpminmax[3];
    $minmaxvals[4]  = $tmpminmax[4];  
    $minmaxvals[5]  = $tmpminmax[5];
    $minmaxvals[6]  = $tmpminmax[6];  
    $minmaxvals[7]  = $tmpminmax[7];
    $minmaxvals[8]  = $tmpminmax[8];  
    $minmaxvals[9]  = $tmpminmax[9];
    $minmaxvals[10] = $tmpminmax[10];  
    $minmaxvals[11] = $tmpminmax[11];
    $minmaxvals[12] = $tmpminmax[12];  
    $minmaxvals[13] = $tmpminmax[13];
    $minmaxvals[14] = $tmpminmax[14];  
    $minmaxvals[15] = $tmpminmax[15];
    $minmaxvals[16] = $tmpminmax[16];  
    $minmaxvals[17] = $tmpminmax[17];
    $minmaxvals[18] = $tmpminmax[18];  
    $minmaxvals[19] = $tmpminmax[19];

    # go through the multiple files and pick out minmax values
    for ($j=1; $j*$ndump<=$Nt; $j++) {
	my $fnum = $j*$ndump;
	print ", $fnum";
	if    ($fnum < 10) {$fname = $exdir."output.001.00000".$fnum;}
	elsif ($fnum < 100) {$fname = $exdir."output.001.0000".$fnum;}
	elsif ($fnum < 1000) {$fname = $exdir."output.001.000".$fnum;}
	elsif ($fnum < 10000) {$fname = $exdir."output.001.00".$fnum;}
	elsif ($fnum < 100000) {$fname = $exdir."output.001.0".$fnum;}
	else                    {$fname = $exdir."output.001.".$fnum;}
	@tmpminmax = tlevel_minmax($fname);
	$minmaxvals[0]  = min($minmaxvals[0],  $tmpminmax[0]);
	$minmaxvals[1]  = max($minmaxvals[1],  $tmpminmax[1]);
	$minmaxvals[2]  = min($minmaxvals[2],  $tmpminmax[2]);
	$minmaxvals[3]  = max($minmaxvals[3],  $tmpminmax[3]);
	$minmaxvals[4]  = min($minmaxvals[4],  $tmpminmax[4]);
	$minmaxvals[5]  = max($minmaxvals[5],  $tmpminmax[5]);
	$minmaxvals[6]  = min($minmaxvals[6],  $tmpminmax[6]);
	$minmaxvals[7]  = max($minmaxvals[7],  $tmpminmax[7]);
	$minmaxvals[8]  = min($minmaxvals[8],  $tmpminmax[8]);
	$minmaxvals[9]  = max($minmaxvals[9],  $tmpminmax[9]);
	$minmaxvals[10] = min($minmaxvals[10], $tmpminmax[10]);
	$minmaxvals[11] = max($minmaxvals[11], $tmpminmax[11]);
	$minmaxvals[12] = min($minmaxvals[12], $tmpminmax[12]);
	$minmaxvals[13] = max($minmaxvals[13], $tmpminmax[13]);
	$minmaxvals[14] = min($minmaxvals[14], $tmpminmax[14]);
	$minmaxvals[15] = max($minmaxvals[15], $tmpminmax[15]);
	$minmaxvals[16] = min($minmaxvals[16], $tmpminmax[16]);
	$minmaxvals[17] = max($minmaxvals[17], $tmpminmax[17]);
	$minmaxvals[18] = min($minmaxvals[18], $tmpminmax[18]);
	$minmaxvals[19] = max($minmaxvals[19], $tmpminmax[19]);
    }
    return @minmaxvals;
}


# @vec = linspace(lowlimit,toplimit,numsteps)
#   subroutine to create linear span into array
sub linspace {
    my $lowlimit = $_[0];
    my $toplimit = $_[1];
    my $numsteps = $_[2];
    my $j;
    my $dx = ($toplimit - $lowlimit) / ($numsteps - 1);
    my @outvec;
    for ($j=0; $j<$numsteps; $j++) {
	$outvec[$j] = $lowlimit + $dx * $j;
    }
    return @outvec;
}


# multiplot(infile,outfile,time,xlabel,ylabel,rhomin,rhomax,umin,umax,
#      vmin,vmax,wmin,wmax,Pmin,Pmax,Bxmin,Bxmax,Bymin,Bymax,Bzmin,
#      Bzmax,dBmin,dBmax,Jmin,Jmax)
#   subroutine to plot files using gnuplot.  Note most arguments passed as
#   the string "_" will be left out or taken care of automatically.
sub multiplot {
    # get plot arguments
    my $ifname = $_[0];   my $ofname = $_[1];
    my $time   = $_[2];
    my $xlabel = $_[3];   my $ylabel = $_[4];
    my $rhomin = $_[5];   my $rhomax = $_[6];
    my $umin   = $_[7];   my $umax   = $_[8];
    my $vmin   = $_[9];   my $vmax   = $_[10];
    my $wmin   = $_[11];  my $wmax   = $_[12];
    my $Pmin   = $_[13];  my $Pmax   = $_[14];
    my $Bxmin  = $_[15];  my $Bxmax  = $_[16];
    my $Bymin  = $_[17];  my $Bymax  = $_[18];
    my $Bzmin  = $_[19];  my $Bzmax  = $_[20];
    my $dBmin  = $_[21];  my $dBmax  = $_[22];
    my $Jmin   = $_[23];  my $Jmax   = $_[24];

    # check to see if bounds taken care of
    my $xminmax = "[ ]";    
    my $yminmax = "[ ]";
    my $rhominmax = ($rhomin eq "_" || $rhomax eq "_") ? "[ ]" : "[$rhomin:$rhomax]";
    my $uminmax   = ($umin   eq "_" || $umax   eq "_") ? "[ ]" : "[$umin:$umax]";
    my $vminmax   = ($vmin   eq "_" || $vmax   eq "_") ? "[ ]" : "[$vmin:$vmax]";
    my $Pminmax   = ($Pmin   eq "_" || $Pmax   eq "_") ? "[ ]" : "[$Pmin:$Pmax]";
    my $Bxminmax  = ($Bxmin  eq "_" || $Bxmax  eq "_") ? "[ ]" : "[$Bxmin:$Bxmax]";
    my $Byminmax  = ($Bymin  eq "_" || $Bymax  eq "_") ? "[ ]" : "[$Bymin:$Bymax]";
    my $dBminmax  = ($dBmin  eq "_" || $dBmax  eq "_") ? "[ ]" : "[$dBmin:$dBmax]";
    my $Jminmax   = ($Jmin   eq "_" || $Jmax   eq "_") ? "[ ]" : "[$Jmin:$Jmax]";

    # set x,y tickmarks
    my $xtics = "set xtics 8; \n";
    my $ytics = "set ytics 4; \n";

    # check if labels taken care of;
#    $xlabel = ($xlabel eq "_") ? "" : "    set xlabel \"$xlabel\" \n";
#    $ylabel = ($ylabel eq "_") ? "" : "    set ylabel \"$ylabel\" \n";
    $xlabel = "";
    $ylabel = "";

    # set variable labels
    my $rholabel = "    set zlabel \"rho\" \n";
    my $ulabel   = "    set zlabel \"u\" \n";
    my $vlabel   = "    set zlabel \"v\" \n";
    my $Plabel   = "    set zlabel \"p\" \n";
    my $Bxlabel  = "    set zlabel \"Bx\" \n";
    my $Bylabel  = "    set zlabel \"By\" \n";
    my $dBlabel  = "    set zlabel \"divB\" \n";
    my $Jlabel   = "    set zlabel \"J\" \n";

    # set titles
    my $title1 = "    set title \"\\n\\n\\n Density\" \n";
    my $title2 = "    set title \"t = $time";
    $title2 = $title2."\\n\\n\\nX-Velocity\" \n";
    my $title3 = "    set title \"\\n\\n\\n Y-Velocity\" \n";
    my $title5 = "    set title \"\\n\\n Pressure\" \n";
    my $title6 = "    set title \"\\n\\n X-Magnetic Field\" \n";
    my $title7 = "    set title \"\\n\\n Y-Magnetic Field\" \n";
    my $title8 = "    set title \"\\n\\n Div(Magnetic Field)\" \n";
    my $title9 = "    set title \"\\n\\n Current\" \n";

    # set text size
    my $textsize = "small";  # choices are small, medium and large

    # set contour location
    my $contourloc = "base";
#    my $contourloc = "surface";
#    my $contourloc = "both";

    # set contour parameters
    my $cntrparam = "set cntrparam bspline; \n";
    $cntrparam = $cntrparam."set cntrparam levels auto 7; \n";
#     my $cntrparam = "set cntrparam bspline; \n";
#     $cntrparam = $cntrparam."set cntrparam points 7; \n";
#     $cntrparam = $cntrparam."set cntrparam order 10; \n";
#     $cntrparam = $cntrparam."set cntrparam levels auto 7; \n";

    # set 3D view
    my $surface = "";
#    my $surface = "set hidden3d; \n";

    # set data sampling for coarser mesh
    my $every   = "every 3:3";
    my $dBevery = "every 6:6";

    # set colors for dark background
    my $bgcolor     = "x161616";
    my $bordercolor = "xededed";
    my $axiscolor   = "xffffff";
    my $color1      = "x456dff"; 
    my $color2      = "xfaf320"; 
    my $color3      = "x00ff00"; 
    my $color4      = "xff7b00";
    my $color5      = "x1af0ff"; 
    my $color6      = "xff00bb"; 
    my $color7      = "x87bdff"; 
    my $color8      = "xff0000"; 

#     # set colors for light background
#     my $bgcolor     = "xffffff";
#     my $bordercolor = "x2b2b2b";
#     my $axiscolor   = "x000000";
#     my $color1      = "x3452c0"; 
#     my $color2      = "xdfd91d"; 
#     my $color3      = "x00c700"; 
#     my $color4      = "xe66f00";
#     my $color8      = "xff0000"; 
#     my $color5      = "x18dae8"; 
#     my $color6      = "xf200b1"; 
#     my $color7      = "x73a2da"; 

    # set desired output size, gnuplot png default is 640x480 (xsize x ysize)
    my $xsize = 1280;  # edit this
    my $ysize = 640;   # edit this
#    my $xsize = 1000;  # edit this
#    my $ysize = 650;   # edit this
    my $xfac  = $xsize/640;  # leave this
    my $yfac  = $ysize/480;  # leave this
    my $xfac14 = $xfac/4;    # leave this
    my $xfac12 = $xfac/2;    # leave this
    my $xfac34 = $xfac*3/4;  # leave this
    my $yfac12 = $yfac/2;    # leave this

    # pipe commands to gnuplot
    open( OUTFILE, '|gnuplot' );
    print OUTFILE "set terminal png $textsize $bgcolor $bordercolor $axiscolor $color1 $color2 $color3 $color4 $color5 $color6 $color7 $color8; \n";
    print OUTFILE "set size $xfac,$yfac; \n";
    print OUTFILE "set output \"$ofname\"; \n";
    print OUTFILE "set data style lines; \n";
    print OUTFILE "set contour $contourloc \n";
    print OUTFILE $cntrparam;
    print OUTFILE $surface;
    print OUTFILE $xtics;
    print OUTFILE $ytics;
    print OUTFILE "set nokey; \n";
    print OUTFILE "set multiplot; \n";
    print OUTFILE "set size $xfac14,$yfac12; \n\n";
    print OUTFILE "set origin 0.0,$yfac12; \n";
    print OUTFILE $title1;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $rholabel;
    print OUTFILE "splot $xminmax $yminmax $rhominmax \"$ifname\" $every using 1:2:3 \n\n";
    print OUTFILE "set origin $xfac14,$yfac12; \n";
    print OUTFILE $title2;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $ulabel;
    print OUTFILE "splot $xminmax $yminmax $uminmax \"$ifname\" $every using 1:2:4 \n\n";
    print OUTFILE "set origin $xfac12,$yfac12; \n";
    print OUTFILE $title3;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $vlabel;
    print OUTFILE "splot $xminmax $yminmax $vminmax \"$ifname\" $every using 1:2:5 \n\n";
    print OUTFILE "set origin $xfac34,$yfac12; \n";
    print OUTFILE $title8;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $dBlabel;
    print OUTFILE "splot $xminmax $yminmax $dBminmax \"$ifname\" $dBevery using 1:2:11 \n\n";
    print OUTFILE "set origin 0.0,0.0; \n";
    print OUTFILE $title5;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $Plabel;
    print OUTFILE "splot $xminmax $yminmax $Pminmax \"$ifname\" $every using 1:2:10 \n\n";
    print OUTFILE "set origin $xfac14,0.0; \n";
    print OUTFILE $title6;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $Bxlabel;
    print OUTFILE "splot $xminmax $yminmax $Bxminmax \"$ifname\" $every using 1:2:7 \n\n";
    print OUTFILE "set origin $xfac12,0.0; \n";
    print OUTFILE $title7;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $Bylabel;
    print OUTFILE "splot $xminmax $yminmax $Byminmax \"$ifname\" $every using 1:2:8 \n\n";
    print OUTFILE "set origin $xfac34,0.0; \n";
    print OUTFILE $title9;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $Jlabel;
    print OUTFILE "splot $xminmax $yminmax $Jminmax \"$ifname\" $every using 1:2:12 \n\n";
    print OUTFILE "set nomultiplot; \n";
    close(OUTFILE);
    
    return 0;
}


# mult_movie( xlabel, ylabel, rhomin, rhomax, umin, umax, vmin, vmax, 
#        wmin, wmax, Pmin, Pmax, Bxmin, Bxmax, Bymin, Bymax, Bzmin, 
#        Bzmax, dBmin, dBmax, Jmin, Jmax, movname, exdir, ndump, @times )
#   subroutine to make movies using gnuplot for the plots, and 
#   mpeg_encode for the movie encoding.
sub mult_movie {
    # get plot arguments
    my $xlabel  = shift;  my $ylabel  = shift;
    my $rhomin  = shift;  my $rhomax  = shift;
    my $umin    = shift;  my $umax    = shift;
    my $vmin    = shift;  my $vmax    = shift;
    my $wmin    = shift;  my $wmax    = shift;
    my $Pmin    = shift;  my $Pmax    = shift;
    my $Bxmin   = shift;  my $Bxmax   = shift;
    my $Bymin   = shift;  my $Bymax   = shift;
    my $Bzmin   = shift;  my $Bzmax   = shift;
    my $dBmin   = shift;  my $dBmax   = shift;
    my $Jmin    = shift;  my $Jmax    = shift;
    my $movname = shift;  my $exdir   = shift;
    my $ndump   = shift;  my @times   = @_;
    my $Nt = $#times;
    
    # movie parameters
    my $PATTERN = 'IBBP';
    my $GOP_SIZE = '10';
    my $SLICES_PER_FRAME = '1';
    my $PIXEL = 'HALF';
    my $RANGE = '10';
    my $PSEARCH_ALG = 'TWOLEVEL';
    my $BSEARCH_ALG = 'CROSS2';
    my $IQSCALE = '5';
    my $PQSCALE = '15';
    my $BQSCALE = '25';
    my $REFERENCE_FRAME = 'ORIGINAL';

    # get rid of old temporary files
    unlink glob "$exdir/tmp*.pnm";

    # start making the movie frames
    print ("  making $movname movie frames\n");
    for (my $j=0; $j*$ndump<=$Nt; $j++) {
	my $fnum = $j*$ndump;
	
	# get time of snapshot
	my $time = $times[$j];

	# plot file
	my $infile;
	if ( $fnum < 10 ) {    $infile  = $exdir."output.001.00000".$fnum; }
	elsif ( $fnum < 100 ) { $infile  = $exdir."output.001.0000".$fnum; }
	elsif ( $fnum < 1000 ) { $infile  = $exdir."output.001.000".$fnum; }
	elsif ( $fnum < 10000 ) { $infile  = $exdir."output.001.00".$fnum; }
	elsif ( $fnum < 100000 ) { $infile  = $exdir."output.001.0".$fnum; }
	else {                      $infile  = $exdir."output.001.".$fnum; }
	my $picfile = $exdir."tmp".$j.".png";
	multiplot( $infile, $picfile, $time, $xlabel, $ylabel, $rhomin, 
		   $rhomax, $umin, $umax, $vmin, $vmax, $wmin, $wmax, $Pmin, 
		   $Pmax, $Bxmin, $Bxmax, $Bymin, $Bymax, $Bzmin, $Bzmax, 
		   $dBmin, $dBmax, $Jmin, $Jmax );
    }
    
    # create movie parameter file
    open( MOVFILE, '>defmovie' );
    print MOVFILE "# automated parameter file for mpeg_encode, \n";
    print MOVFILE "# created by perl script movies.pl \n\n";
    print MOVFILE "PATTERN           $PATTERN \n";
    print MOVFILE "OUTPUT            $movname \n\n";
    print MOVFILE "BASE_FILE_FORMAT  PNM \n";
    print MOVFILE "INPUT_CONVERT     pngtopnm * \n";
    print MOVFILE "GOP_SIZE          $GOP_SIZE \n";
    print MOVFILE "SLICES_PER_FRAME  $SLICES_PER_FRAME \n";
    print MOVFILE "FRAME_RATE        23.976 \n";
    print MOVFILE "ASPECT_RATIO      1.2015 \n";
    print MOVFILE "FORCE_ENCODE_LAST_FRAME \n\n";
    print MOVFILE "INPUT_DIR         ./$exdir \n\n";
    print MOVFILE "INPUT \n";
    print MOVFILE "tmp*.png [0-$Nt]\n";
    print MOVFILE "END_INPUT \n\n";
    print MOVFILE "PIXEL             $PIXEL \n";
    print MOVFILE "RANGE             $RANGE \n\n";
    print MOVFILE "PSEARCH_ALG       $PSEARCH_ALG \n";
    print MOVFILE "BSEARCH_ALG       $BSEARCH_ALG \n\n";
    print MOVFILE "IQSCALE           $IQSCALE \n";
    print MOVFILE "PQSCALE           $PQSCALE \n";
    print MOVFILE "BQSCALE           $BQSCALE \n\n";
    print MOVFILE "REFERENCE_FRAME   $REFERENCE_FRAME \n";
    print MOVFILE "USER_DATA         /dev/null \n";
    close(MOVFILE);
    
    # create movie
    print ("  encoding movie $movname\n");
    system "mpeg_encode -realquiet defmovie";
 
    # get rid of temporary files
    unlink glob "defmovie $exdir/tmp*.png";

    return 0;
}


# single_plot(infile,outfile,time,xlabel,ylabel,zlabel,zmin,zmax,title,fnumber)
#   subroutine to plot fluid number 'fnumber' using gnuplot.  Note most 
#   arguments passed as the string "_" will be left out or taken care of 
#   automatically.
sub single_plot {
    # get plot arguments
    my $ifname = $_[0];   my $ofname  = $_[1];
    my $time   = $_[2];   my $xlabel  = $_[3];
    my $ylabel = $_[4];   my $zlabel  = $_[5];
    my $zmin   = $_[6];   my $zmax    = $_[7];
    my $title  = $_[8];   my $fnumber = $_[9];

    # check to see if bounds taken care of
    my $xminmax = "[ ]";    
    my $yminmax = "[ ]";
    my $zminmax = ($zmin eq "_" || $zmax eq "_") ? "[ ]" : "[$zmin:$zmax]";

    # set x,y tickmarks
    my $xtics = "set xtics 4; \n";
    my $ytics = "set ytics 2; \n";

    # check if labels taken care of;
    $xlabel = ($xlabel eq "_") ? "" : "    set xlabel \"$xlabel\" \n";
    $ylabel = ($ylabel eq "_") ? "" : "    set ylabel \"$ylabel\" \n";
    $zlabel = ($zlabel eq "_") ? "" : "    set zlabel \"$zlabel\" \n";
    $title  = ($title eq "_")  ? "    set title \"t = $time\" \n" :
	"    set title \"$title,   t = $time\" \n";

    # set text size
    my $textsize = "small";  # choices are small, medium and large

    # set contour location
    my $contourloc = "base";
#    my $contourloc = "surface";
#    my $contourloc = "both";

    # set contour parameters
    my $cntrparam = "set cntrparam bspline; \n";
    $cntrparam = $cntrparam."set cntrparam points 7; \n";
    $cntrparam = $cntrparam."set cntrparam order 10; \n";
    $cntrparam = $cntrparam."set cntrparam levels auto 7; \n";

    # set 3D view
    my $surface = "";
#    my $surface = "set hidden3d; \n";

    # set data sampling for coarser mesh
    my $every = "every 2:2";

    # set colors for dark background
    my $bgcolor     = "x161616";
    my $bordercolor = "xededed";
    my $axiscolor   = "xffffff";
    my $color1      = "x456dff"; 
    my $color2      = "xfaf320"; 
    my $color3      = "x00ff00"; 
    my $color4      = "xff7b00";
    my $color5      = "x1af0ff"; 
    my $color6      = "xff00bb"; 
    my $color7      = "x87bdff"; 
    my $color8      = "xff0000"; 

#     # set colors for light background
#     my $bgcolor     = "xffffff";
#     my $bordercolor = "x2b2b2b";
#     my $axiscolor   = "x000000";
#     my $color1      = "x3452c0"; 
#     my $color2      = "xdfd91d"; 
#     my $color3      = "x00c700"; 
#     my $color4      = "xe66f00";
#     my $color8      = "xff0000"; 
#     my $color5      = "x18dae8"; 
#     my $color6      = "xf200b1"; 
#     my $color7      = "x73a2da"; 

    # set desired output size, gnuplot png default is 640x480 (xsize x ysize)
    my $xsize = 640;   # edit this
    my $ysize = 480;   # edit this
    my $xfac  = $xsize/640;  # leave this
    my $yfac  = $ysize/480;  # leave this

    # create temporary output files
    open( OUTFILE, '|gnuplot' );
    print OUTFILE "set terminal png $textsize $bgcolor $bordercolor $axiscolor $color1 $color2 $color3 $color4 $color5 $color6 $color7 $color8; \n";
    print OUTFILE "set size $xfac,$yfac; \n";
    print OUTFILE "set output \"$ofname\"; \n";
    print OUTFILE "set data style lines; \n";
    print OUTFILE "set contour $contourloc \n";
    print OUTFILE $cntrparam;
    print OUTFILE $surface;
    print OUTFILE $xtics;
    print OUTFILE $ytics;
    print OUTFILE "set nokey; \n";
    print OUTFILE $title;
    print OUTFILE $xlabel;
    print OUTFILE $ylabel;
    print OUTFILE $zlabel;
    print OUTFILE "splot $xminmax $yminmax $zminmax \"$ifname\" $every using 1:2:$fnumber \n\n";
    close(OUTFILE);
    
    return 0;
}



# single_movie( xlabel, ylabel, zlabel, zmin, zmax, title, fnumber, 
#               movname, exdir, ndump, @times )
#   subroutine to make single-fluid movie using gnuplot for the plots, 
#   and mpeg_encode for the movie encoding.
sub single_movie {
    # get plot arguments
    my $xlabel  = shift;  my $ylabel  = shift;
    my $zlabel  = shift;  my $zmin    = shift;  
    my $zmax    = shift;  my $title   = shift;
    my $fnumber = shift;  my $movname = shift;  
    my $exdir   = shift;  my $ndump   = shift;  
    my @times   = @_;
    my $Nt = $#times;
    
    # movie parameters
    my $PATTERN = 'IBBP';
    my $GOP_SIZE = '10';
    my $SLICES_PER_FRAME = '1';
    my $PIXEL = 'HALF';
    my $RANGE = '10';
    my $PSEARCH_ALG = 'TWOLEVEL';
    my $BSEARCH_ALG = 'CROSS2';
    my $IQSCALE = '5';
    my $PQSCALE = '15';
    my $BQSCALE = '25';
    my $REFERENCE_FRAME = 'ORIGINAL';

    # get rid of old temporary files
    unlink glob "$exdir/tmp*.pnm";

    # start making the movie frames
    print ("  making $movname movie frames\n");
    for (my $j=0; $j*$ndump<=$Nt; $j++) {
	
	# get time of snapshot
	my $time = $times[$j];

	# plot file
	my $infile;
	my $jnd = $j*$ndump;
	if ( $jnd < 10 ) {    $infile  = $exdir."output.001.00000".$jnd; }
	elsif ( $jnd < 100 ) { $infile  = $exdir."output.001.0000".$jnd; }
	elsif ( $jnd < 1000 ) { $infile  = $exdir."output.001.000".$jnd; }
	elsif ( $jnd < 10000 ) { $infile  = $exdir."output.001.00".$jnd; }
	elsif ( $jnd < 100000 ) { $infile  = $exdir."output.001.0".$jnd; }
	else {                     $infile  = $exdir."output.001.".$jnd; }
	my $picfile = $exdir."tmp".$j.".png";
	single_plot( $infile, $picfile, $time, $xlabel, $ylabel, $zlabel, 
		     $zmin, $zmax, $title, $fnumber );
    }
    
    # create movie parameter file
    open( MOVFILE, '>defmovie' );
    print MOVFILE "# automated parameter file for mpeg_encode, \n";
    print MOVFILE "# created by perl script movies.pl \n\n";
    print MOVFILE "PATTERN           $PATTERN \n";
    print MOVFILE "OUTPUT            $movname \n\n";
    print MOVFILE "BASE_FILE_FORMAT  PNM \n";
    print MOVFILE "INPUT_CONVERT     pngtopnm * \n";
    print MOVFILE "GOP_SIZE          $GOP_SIZE \n";
    print MOVFILE "SLICES_PER_FRAME  $SLICES_PER_FRAME \n";
    print MOVFILE "FRAME_RATE        23.976 \n";
    print MOVFILE "ASPECT_RATIO      1.2015 \n";
    print MOVFILE "FORCE_ENCODE_LAST_FRAME \n\n";
    print MOVFILE "INPUT_DIR         ./$exdir \n\n";
    print MOVFILE "INPUT \n";
    print MOVFILE "tmp*.png [0-$Nt]\n";
    print MOVFILE "END_INPUT \n\n";
    print MOVFILE "PIXEL             $PIXEL \n";
    print MOVFILE "RANGE             $RANGE \n\n";
    print MOVFILE "PSEARCH_ALG       $PSEARCH_ALG \n";
    print MOVFILE "BSEARCH_ALG       $BSEARCH_ALG \n\n";
    print MOVFILE "IQSCALE           $IQSCALE \n";
    print MOVFILE "PQSCALE           $PQSCALE \n";
    print MOVFILE "BQSCALE           $BQSCALE \n\n";
    print MOVFILE "REFERENCE_FRAME   $REFERENCE_FRAME \n";
    print MOVFILE "USER_DATA         /dev/null \n";
    close(MOVFILE);
    
    # create movie
    print ("  encoding movie $movname\n");
    system "mpeg_encode -realquiet defmovie";
 
    # get rid of temporary files
    unlink glob "defmovie $exdir/tmp*.png";

    return 0;
}



######################################################
#  End of script
######################################################
