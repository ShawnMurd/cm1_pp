#!/opt/local/bin/perl
##############################################################
#           TOFFSET_INDEXER.PL                               #
#                                                            #
# For cm1r17 and earlier releases, if a model run is started #
# from a restart file, but run with a different time step,   #
# the model time in the output files will be incorrect. This #
# script produces an index of which file is at which time.   #
# The first column is the incorrect time, the second column  #
# is the correct time, and the third column is the filename. #
##############################################################
$f1 = 27; # number of first file
$f2 = 2626; # number of last file

$t1=78; # time of first file
$tsb=3003; # time it should be
$dt = 3; # time increment

############################

$toffset = $tsb-$t1;

open(FILE,">time_index.txt");

$t=$t1;
for($f=$f1;$f<=$f2;$f++) {
  $ftag = sprintf('%06d',$f);
  $t2 = $t+$toffset;
  print FILE "cm1out_$ftag.nc\t$t\t$t2\n";
  $t+=$dt;
}
close(FILE);
exit;  
