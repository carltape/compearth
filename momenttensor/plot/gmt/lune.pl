#!/usr/bin/perl -w

#==========================================================
#
#  Perl script to write a shell script for GMT plotting.
#  Carl Tape, carltape@gi.alaska.edu, May 2012
#
#  Example plot of representing moment tensors on the funamental lune.
#  Reference: W. Tape and C. Tape, GJI, 2012
#
#  Last tested 5-25-2012 with GMT 4.5.3
#  
#==========================================================

use Math::Trig;

$cshfile = "lune.csh";

$fontno = "1"; 

open(CSH,">$cshfile");
print CSH "gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.3c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p\n";

$R = "-R-30/30/-90/90";
#$R = "-R-30/30/-90/90r";
$origin = "-X2 -Y1";
$xtick1 = 10; $ytick1 = 10;
$xtick2 = 5; $ytick2 = 5;
#$B = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":WesN";
$B = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":wesn";

$wid = "3i";
#$J = "-JA0/0/$wid"; $title = "Lambert equal-area ($J)"; $ftag = "lambert";
#$J = "-JB0/0/0/90/$wid"; $title = "Albers equal-area ($J)"; $ftag = "albers";

#$J = "-JY0/$wid";  $title1 = "Cylindrical equal-area ($J)"; $ftag = "cylindrical";
#$J = "-JI0/$wid";  $title1 = "Sinusoidal equal-area ($J)"; $ftag = "sinusoidal";
#$J = "-JKf0/$wid"; $title1 = "Eckert IV equal-area ($J)"; $ftag = "eckert4";
#$J = "-JK0/$wid";  $title1 = "Eckert VI equal-area ($J)"; $ftag = "eckert6";
#$J = "-JW0/$wid";  $title1 = "Mollweide equal-area ($J)"; $ftag = "mollewide";
$J = "-JH0/$wid";   $title1 = "Hammer equal-area ($J)"; $ftag = "hammer";

$title1 = "Representation of source types on the fundamental lune";
$title2 = "(W. Tape and C. Tape, 2012, GJI, \"A geometric setting for moment tensors\")";

# colors
#$magenta = "148/0/211";
$magenta = "160/32/240";
$orange = "255/165/0";
$red = "255/0/0";
$blue = "30/144/255";
$cyan = "0/255/255";
$green = "50/205/50";

$lgray = 200;
$dgray = 120;

# KEY COMMAND
$iplot = 1;  # =0 (reference lune), =1 (dots from punlished studies), =2 (reference beachballs)
$kplot = 1;  # orientation of MT at center of lune (iplot=2 only)
$psfile = "lune_${ftag}_iplot${iplot}_kplot${kplot}.ps";

print CSH "psbasemap $J $R $B -G$lgray -K -V -P $origin > $psfile\n"; # START

$pdir = "./dfiles/";

# plot patches
$fname = "$pdir/beach_patch_01.lonlat";
print CSH "psxy $fname -G$dgray -J -R -K -O -V >>$psfile\n";
$fname = "$pdir/beach_patch_02.lonlat";
print CSH "psxy $fname -G255 -J -R -K -O -V >>$psfile\n";

print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";

# plot arcs
# dev, iso+DC, iso, iso, CDC
$lwid = 3;
$fname1 = "$pdir/beach_arc_01.lonlat";
$fname2 = "$pdir/beach_arc_02.lonlat";
$fname3 = "$pdir/beach_arc_03.lonlat";
$fname4 = "$pdir/beach_arc_04.lonlat";
$fname5 = "$pdir/beach_arc_05.lonlat";
if ($iplot==1) {
  $W = "-W2p,0,--";
  print CSH "psxy $fname1 $W -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname2 $W -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname3 $W -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname4 $W -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname5 $W -J -R -K -O -V >>$psfile\n";
} else {
  print CSH "psxy $fname1 -W${lwid}p,$magenta -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname2 -W${lwid}p,$red -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname3 -W${lwid}p,$blue -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname4 -W${lwid}p,$blue -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname5 -W${lwid}p,0/0/0 -J -R -K -O -V >>$psfile\n";
}

if ($iplot != 2) {
  # plot points
  $csize = 12;
  $fname = "$pdir/beach_points.lonlat";
  print CSH "psxy $fname -N -Sc${csize}p -W1p,0/0/0 -G255 -J -R -K -O -V >>$psfile\n";

  # plot labels
  $fsize = 14;
  $fontno = 1;
  print "$fname\n";
  open(IN,$fname); @plines = <IN>; close(IN);
  for ($i = 1; $i <= @plines; $i++) {
    ($plon,$plat,$plab,$Dx,$Dy) = split(" ",$plines[$i-1]);
    #print "\n--$plon -- $plat-- $plab --";
    $D = "-D${Dx}p/${Dy}p";
    print CSH "pstext -N -J -R -K -O -V $D >>$psfile<<EOF\n$plon $plat $fsize 0 $fontno CM $plab\nEOF\n";
  }
  #print CSH "awk '{print \$1,\$2,$fsize,0,0,\"CM\",\$3}' $fname | pstext -N -J -R -K -O -V >> $psfile\n";
}

if ($iplot==1) {
  # moment tensors from various studies
  $csize = 8;
  $fname = "$pdir/beachpts_Ford2009_points.dat";
  print CSH "psxy $fname -N -Sc${csize}p -W0.5p,0/0/0 -G$red -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/beachpts_Foulger2004_points.dat";
  print CSH "psxy $fname -N -Sc${csize}p -W0.5p,0/0/0 -G$orange -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/beachpts_Minson2007_points.dat";
  print CSH "psxy $fname -N -Sc${csize}p -W0.5p,0/0/0 -G$green -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/beachpts_Walter2009_points.dat";
  print CSH "psxy $fname -N -Sc${csize}p -W0.5p,0/0/0 -G$cyan -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/beachpts_Walter2010_points.dat";
  print CSH "psxy $fname -N -Sc${csize}p -W0.5p,0/0/0 -G$magenta -J -R -K -O -V >>$psfile\n";

} elsif ($iplot==2) {
  # reference beachballs on the lune
  $cmtinfo = "-Sm0.5 -L0.5p/0/0/0 -G255/0/0 -N";
  $cmtfile = sprintf("$pdir/beachballs_%i_psmeca",$kplot);
  print CSH "psmeca $cmtfile $J $R $cmtinfo -K -O -V >> $psfile\n";
} 

#-----------------------------

$J_title = "-JX1i";  # -JM7i
$R_title = "-R0/1/0/1";
$otitle1 = "-Xa3.5 -Ya6";

# legend for plotting published studies
if($iplot==1) {
print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G$magenta $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n0 1.2\nEOF\n";
print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0.2 1.2 12 0 $fontno LM Walter 2010 (n=14)\nEOF\n";

print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G$cyan $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n0 0.9\nEOF\n";
print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0.2 0.9 12 0 $fontno LM Walter 2009 (n=13)\nEOF\n";

print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G$green $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n0 0.6\nEOF\n";
print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0.2 0.6 12 0 $fontno LM Minson 2007 (n=18)\nEOF\n";

print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G$orange $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n0 0.3\nEOF\n";
print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0.2 0.3 12 0 $fontno LM Foulger 2004 (n=26)\nEOF\n";

print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G$red $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n0 0\nEOF\n";
print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0.2 0 12 0 $fontno LM Ford 2009 (n=32)\nEOF\n";

#$x = -30; $y = rad2deg(asin(1/sqrt(3)));
#print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G255/165/0 -J -R -K -O -V >>$psfile<<EOF\n$x $y\nEOF\n";
}

#-----------------------------

$otitle1 = "-Xa-1 -Ya9.0";
$otitle2 = "-Xa-1 -Ya8.7";
print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0 0 14 0 $fontno LM $title1\nEOF\n";
print CSH "pstext -N $R_title $J_title $otitle2 -O -V >>$psfile<<EOF\n 0 0 11 0 $fontno LM $title2\nEOF\n";

close (CSH);
system("csh -f $cshfile");

system("ps2pdf $psfile");

# you may need to install gv to view (or use something else)
system("gv $psfile &");

#==================================================
