#!/usr/bin/perl
############################################
# BATCH_DYN -- batch-processes a group of  #
# CM1 files using dyn.                     #
# -----------------------------------------#
# Ryan Hastings, 13 January 2013           #
#------------------------------------------#
# Cleaned up and commented, 25 Dec 2016,   #
# Ryan Hastings                            #
############################################

#-------- configuration variables ---------#
# files to process.  this assumes the output
# files are in the form dyn_000001.nc (i.e.,
# "dyn_" followed by six digits padded by
# zeroes, followed by ".nc"
$f1 =  642; # number of first file
$f2 = 1036; # number of last file
$fstep = 1; # step between files

# dyn configuration
$dz = 500; # yes, it's called dz, but actually put in
           # dx or dy
$output_dl = 1; # output linear dynamic pressure fields?
$output_dn = 1; # output nonlinear dynamic pressure fields?
$output_b  = 1; # output buoyancy pressure fields?

$dyn="/convect/s1/rmh265/cm1_pp_v0.9/bin/dyn"; # location of dyn on your system
$ptype = 5; # microphysics type

#------------- loop through files ---------#
for ($f=$f1;$f<=$f2;$f+=$fstep) {

        $tag = sprintf('%06d',$f);
	$cm1name = "cm1out_$tag.nc";
	$inname  = "dyn_$tag.input";
	$dynname = "dyn_$tag.nc";

	#---- write out input file --------#
	open(INPUT,">$inname");

	print INPUT "&inputs\n  infile = \'$cm1name\',\n";
	print INPUT "  outfile = \'$dynname\',\n  dz = $dz,\n";
        print INPUT "  ptype = \'$ptype\',\n/\n";

	print INPUT "\n&outputs\n  output_dl = $output_dl,\n";
  	print INPUT "  output_dn = $output_dn,\n";
        print INPUT "  output_b  = $output_b,\n/\n";

	close(INPUT);

	#------- execute dyn --------------#
	system("mv $inname dyn.input");
	system("$dyn");
	system("mv dyn.input $inname");

}


exit;
