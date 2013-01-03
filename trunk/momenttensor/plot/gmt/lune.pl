#!/usr/bin/perl -w

#==========================================================
#
#  Perl script to write a shell script for GMT plotting.
#  Carl Tape, carltape@gi.alaska.edu, May 2012
#
#  Example plot of representing moment tensors on the funamental lune.
#  The basic concepts behind this representation of moment tensors can be found in
#  W. Tape and C. Tape, "A geometric setting for moment tensors," Geophysical J. International, 2012
#
#  Last tested 8-22-2012 with GMT 4.5.8, but with a custom psmeca by Doug Dreger
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

# colors
#$magenta = "148/0/211";
$magenta = "160/32/240";
$orange = "255/165/0";
$red = "255/0/0";
$blue = "30/144/255";
$cyan = "0/255/255";
$green = "50/205/50";
$sienna = "160/82/45";
$brown = "139/69/16";
$black = "0/0/0";
$white = "255/255/255";

$lgray = 200;
$dgray = 120;

# KEY COMMAND
$iplot = 1;  # =0 (reference lune), =1 (dots from published studies), =2 (reference beachballs)
$lplot = 1;  # =1-2: reference MTs on the lune (iplot=2 only)
$kplot = 3;  # =1-4: orientation of MT at center of lune (iplot=2 only)
if($iplot==2) {
  $psfile = "lune_${ftag}_iplot${iplot}_lplot${lplot}_kplot${kplot}.ps";
} else {
  $psfile = "lune_${ftag}_iplot${iplot}.ps";
}
$ipatch = 1;  # three shaded patches on the lune
$icrack = 1;  # nu=0.25 arc between crack points
if ($iplot==0) {$plot_ref_points = 1; $plot_ref_labels = 1;}
if ($iplot==1) {$plot_ref_points = 0; $plot_ref_labels = 1;}
if ($iplot==2) {$plot_ref_points = 0; $plot_ref_labels = 0;}

$clune = $lgray;
if($ipatch==0) {$clune = $sienna;}

print CSH "psbasemap $J $R $B -G$clune -K -V -P $origin > $psfile\n"; # START

$pdir = "./dfiles/";

# plot patches
if ($ipatch==1) {
  $fname = "$pdir/beach_patch_01.lonlat";
  print CSH "psxy $fname -G$dgray -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/beach_patch_02.lonlat";
  print CSH "psxy $fname -G255 -J -R -K -O -V >>$psfile\n";
  print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";
}

# plot arcs
# dev, iso+DC, iso, iso, CDC nu=0.25, CDC nu=0
$lwid = 3;
$fname1 = "$pdir/beach_arc_01.lonlat";
$fname2 = "$pdir/beach_arc_02.lonlat";
$fname3 = "$pdir/beach_arc_03.lonlat";
$fname4 = "$pdir/beach_arc_04.lonlat";
$fname5 = "$pdir/beach_arc_05.lonlat";
$fname6 = "$pdir/beach_arc_06.lonlat";
if ($iplot==1) {
  $W = "-W2p,0,--";
  print CSH "psxy $fname1 $W -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname2 $W -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname3 $W -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname4 $W -J -R -K -O -V >>$psfile\n";
  if($icrack==1) {print CSH "psxy $fname5 $W -J -R -K -O -V >>$psfile\n";}
  print CSH "psxy $fname6 $W -J -R -K -O -V >>$psfile\n";
} else {
  if($ipatch==1) {@cols = ($magenta,$red,$blue,$blue,$black,$blue);}
  else           {@cols = ($magenta,$orange,$red,$white,$black,$blue);}
  if($ipatch==1) {print CSH "psxy $fname1 -W${lwid}p,$cols[0] -J -R -K -O -V >>$psfile\n";}
  if($ipatch==1) {print CSH "psxy $fname2 -W${lwid}p,$cols[1] -J -R -K -O -V >>$psfile\n";}
  print CSH "psxy $fname3 -W${lwid}p,$cols[2] -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname4 -W${lwid}p,$cols[3] -J -R -K -O -V >>$psfile\n";
  if($icrack==1) {print CSH "psxy $fname5 -W${lwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  print CSH "psxy $fname6 -W${lwid}p,$cols[5] -J -R -K -O -V >>$psfile\n";
}

# plot lune reference points and labels
if ($plot_ref_points || $plot_ref_labels) {
  $fname = "$pdir/beach_points.lonlat";
  if ($plot_ref_points) {
    $csize = 12;
    print CSH "psxy $fname -N -Sc${csize}p -W1p,0/0/0 -G255 -J -R -K -O -V >>$psfile\n";
  }
  if ($plot_ref_labels) {
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
}

if ($iplot==1) {
  # moment tensors from various studies
  $csize = 8;  # size of dots
  @csizes = ($csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize/2);
  @cols = ($red,$orange,$green,$white,$cyan,$magenta,$black,$green);
  @ftags = ("Ford2009","Foulger2004","Minson2007","Minson2008","Walter2009","Walter2010","Pesicek2012","Baig2010");
  @ftits = ("Ford 2009 (n=32)","Foulger 2004 (n=26)","Minson 2007 (n=18)","Minson 2008 (n=7)","Walter 2009 (n=13)","Walter 2010 (n=14)","Pesicek 2012 (n=7)","Baig 2010 (n=577)");
  @inds = (1..7);
  #@inds = (8,1..7);
  #@inds = 8;

  for ($i = 1; $i <= @inds; $i++) {
    $j = $inds[$i-1];
    $cz = $csizes[$j-1];
    $fname = sprintf("$pdir/beachpts_%s_points.dat",$ftags[$j-1]);
    print CSH "psxy $fname -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] -J -R -K -O -V >>$psfile\n";
  }

} elsif ($iplot==2) {
  # reference beachballs on the lune
  $cmtinfo = "-Sm0.5 -L0.5p/0/0/0 -G255/0/0 -N";
  $cmtfile = sprintf("$pdir/beachballs_ilune%i_iref%i_psmeca",$lplot,$kplot);
  print CSH "psmeca $cmtfile $J $R $cmtinfo -K -O -V >> $psfile\n";
} 

#-----------------------------

$J_title = "-JX1i";  # -JM7i
$R_title = "-R0/1/0/1";
$otitle1 = "-Xa3.5 -Ya6";

# legend for plotting published studies
if($iplot==1) {
  $x0 = 0; $y0 = 1.2; $dy = 0.3;
  for ($i = 1; $i <= @inds; $i++) {
    $j = $inds[$i-1];
    $cz = $csizes[$j-1];
    $x = $x0 + 0.2;
    $y = $y0 - ($i-1)*$dy;
    print CSH "psxy -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n$x0 $y\nEOF\n";
    print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n $x $y 12 0 $fontno LM $ftits[$j-1]\nEOF\n";
  }

#$x = -30; $y = rad2deg(asin(1/sqrt(3)));
#print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G255/165/0 -J -R -K -O -V >>$psfile<<EOF\n$x $y\nEOF\n";
}

#-----------------------------

# optional: plot a title
$otitle1 = "-Xa-1 -Ya9.0";
$otitle2 = "-Xa-1 -Ya8.7";
if (0==1) {
  $title1 = "Representation of source types on the fundamental lune";
  $title2 = "(W. Tape and C. Tape, 2012, GJI, \"A geometric setting for moment tensors\")";
  print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0 0 14 0 $fontno LM $title1\nEOF\n";
  print CSH "pstext -N $R_title $J_title $otitle2 -O -V >>$psfile<<EOF\n 0 0 11 0 $fontno LM $title2\nEOF\n";
} else {
  # any command with no -K should work (this one will not plot anything)
  print CSH "pstext $R_title $J_title $otitle2 -O -V >>$psfile<<EOF\n -1 -1 11 0 $fontno LM TEST\nEOF\n"; 
}

close (CSH);
system("csh -f $cshfile");

system("ps2pdf $psfile");

# you may need to install gv to view (or use something else)
system("gv $psfile &");

# to make composite pdf file:
# for file in `ls lune_hammer_iplot2_*.ps` ; do ps2pdf $file ; done ; pdcat -r lune_hammer_iplot2*pdf all_lune_hammer_iplot2.pdf
# pdcat -r lune_hammer_iplot0_kplot1.pdf lune_hammer_iplot1_kplot1.pdf all_lune_hammer_iplot2.pdf all_lune.pdf

#==================================================
