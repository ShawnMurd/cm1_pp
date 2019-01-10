#!/usr/bin/perl
######################################################
##             PARCEL_LIST_MAKER.PL                 ##
##                                                  ##
## Create a list of parcels as input for traj.  The ##
## parcels comprise a box from x1-x2, y1-y2 and     ##
## z1-z2.                                           ##
##==================================================##
## v1.0, 25 Dec 2016, Ryan Hastings                 ##
######################################################

#------------------ configuration variables ---------#

# NOTE:  The x-, y-, and z-spacing CANNOT equal zero,
# even if you're doing it on a 2-D plane.  If, for
# example, you're doing a horizontal plane, set $z1 and
# $z2 to the same value, and set $dz to any nonzero number.

# x-dimensions
$x1 = -15000; # western edge of box
$x2 =      0; # eastern edge of box
$dx =    500; # x-spacing

$y1 =      0; # southern edge of box
$y2 =  15000; # northern edge of box
$dy =    500; # y-spacing

$z1 =   25; # bottom edge of box
$z2 =   25; # top edge of box
$dz =    150; # z-spacing

$outfile = "parcels.input"; # output filename

#------------------ execution part -------------------#

$np = (($x2-$x1)/$dx+1)*(($y2-$y1)/$dy+1)*(($z2-$z1)/$dz+1);

open(PARCELS,">$outfile");
print PARCELS "$np\n";

for ($x=$x1;$x<=$x2;$x+=$dx) {
for ($y=$y1;$y<=$y2;$y+=$dy) {
for ($z=$z1;$z<=$z2;$z+=$dz) {
	print PARCELS "$x\t$y\t$z\n";
}
}
}

close(PARCELS);

exit;
