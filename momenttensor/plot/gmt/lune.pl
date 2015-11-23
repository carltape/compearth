#!/usr/bin/perl -w

#===================================================================================
#
#  lune.pl
#  Perl script to write a shell script for GMT plotting.
#  Carl Tape, ctape@alaska.edu, May 2012
#
#  Plot moment tensors (or source types) on the fundamental lune.
#  The basic concepts behind this representation of moment tensors can be found in TT2012
#  See also examples in TT2013
#
#  Set kmax=2 to also plot source types on the vw rectangle plot of TT2015
#  
#  Abbreviations
#    gCDC = generalized crack-plus-double-couple model (TT2013)
#
#  References:
#  TT2012: W. Tape and C. Tape, "A geometric setting for moment tensors," Geophysical J. International, 2012
#  TT2013: W. Tape and C. Tape, "The classical model for moment tensors," Geophysical J. International, 2013
#  TT2015: W. Tape and C. Tape, "A uniform parameterization for seismic moment tensors," Geophysical J. International, 2015
#
#  Last tested 11-08-2015 with GMT 4.5.3
#  Our psmeca (for plotting beachballs) uses a modified version of utilmeca.c by Doug Dreger.
#  
#===================================================================================

# all we use this for is pi
use Math::Trig;
use Math::Trig ':pi';

$cshfile = "lune.csh";

$fontno = "1"; 
$ticklen = "0.0c";    # =0 for no ticks

open(CSH,">$cshfile");
# set GMT default parameters
print CSH "gmtset PAPER_MEDIA letter PAGE_ORIENTATION landscape MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $ticklen LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10 HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p\n";

#----------
# lune source-type plot

$Rlune = "-R-30/30/-90/90";
#$Rlune = "-R-30/30/-90/90r";
$origin = "-X2.5 -Y0.5";
$xtick1 = 10; $ytick1 = 10;
$xtick2 = 5; $ytick2 = 5;
#$Blune = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":WesN";
$Blune = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":wesn";
#$Blune = "-Ba${xtick1}f${xtick2}:\" \":/a${ytick1}f${xtick2}:\" \":wesn";   # no gridlines

$wid = "2.8i";
#$Jlune = "-JA0/0/$wid"; $title = "Lambert equal-area ($Jlune)"; $ftag = "lambert";
#$Jlune = "-JB0/0/0/90/$wid"; $title = "Albers equal-area ($Jlune)"; $ftag = "albers";

#$Jlune = "-JY0/$wid";  $title1 = "Cylindrical equal-area ($Jlune)"; $ftag = "cylindrical";
#$Jlune = "-JI0/$wid";  $title1 = "Sinusoidal equal-area ($Jlune)"; $ftag = "sinusoidal";
#$Jlune = "-JKf0/$wid"; $title1 = "Eckert IV equal-area ($Jlune)"; $ftag = "eckert4";
#$Jlune = "-JK0/$wid";  $title1 = "Eckert VI equal-area ($Jlune)"; $ftag = "eckert6";
#$Jlune = "-JW0/$wid";  $title1 = "Mollweide equal-area ($Jlune)"; $ftag = "mollewide";
$Jlune = "-JH0/$wid";   $title1 = "Hammer equal-area ($Jlune)"; $ftag = "hammer";

#----------
# vw rectangle source-type plot

$vmax = 1/3;
$vmin = -$vmax;
$wmax = 3*pi/8;
$wmin = -$wmax;
$vran = $vmax - $vmin;
$wran = $wmax - $wmin;

$Rrect = "-R$vmin/$vmax/$wmin/$wmax";
$xtick1 = $vmax; $ytick1 = $wmax;
$xtick2 = 0.1; $ytick2 = 0.1;
$Brect = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":wesn";

$wid = 2.1;
$ywid = $wid*$wran/$vran;
$Jrect = "-JX${wid}i/${ywid}i";

#----------

# colors
#$magenta = "148/0/211";
$magenta = "160/32/240";
$orange = "255/165/0";
$red = "255/0/0";
$pink = "255/150/150";
$blue = "30/144/255";
$cyan = "0/255/255";
$green = "50/205/50";
$sienna = "160/82/45";
$brown = "139/69/16";
$yellow = "255/255/0";
$lred = "255/150/150";
$lpurple = "132/112/255";

$black = "0/0/0";
$white = "255/255/255";
$lgray = 200;
$dgray = 120;

# Poisson values considered for gCDC model
@nus = (0.25,0.36);

# PLOTTING OPTIONS
$iplot = 1;  # =0 (reference lune), =1 (dots from published studies), =2 (reference beachballs), =3 (catalog of full moment tensors)
$lplot = 1;  # =1-2: reference MTs on the lune (iplot=2 only)
$kplot = 1;  # =1-4: orientation of MT at center of lune (iplot=2 only)
if($iplot==2) {
  $psfile = "lune_${ftag}_iplot${iplot}_lplot${lplot}_kplot${kplot}.ps";
} else {
  $psfile = "lune_${ftag}_iplot${iplot}.ps";
}
# points
$plot_ref_points = 1;  # plot reference points on lune (ISO, DC, etc)
$plot_ref_labels = 1;  # reference labels: =0 (none), =1 (eigs), =2 (ISO, DC, etc)
$iiplotDC = 0;
# patches for lam_i = 0 regions
$ipatch = 1;
# arcs
$idev = 0;             # deviatoric arc
$inup25 = 1;           # nu=0.25 arc between crack points
$inup36 = 0;           # nu=0.36 arc between crack points
$ilam3 = 1;            # lam3 = 0 arc
$ilam2 = 1;            # lam2 = 0 arc between dipoles
$ilam1 = 1;            # lam1 = 0 arc
# gCDC options
$ipatchgcdc = 0;       # patches for gCDC model (inugcdc > 0)
$iarcgcdc = 0;         # boundary arcs for gCDC model (inugcdc > 0)
@nus = (0.25,0.36);    # Poisson values considered for gCDC model
$inugcdc = 1;          #   =1 for nu=0.25, =2 fo0.36 (ice)
# other options
$ilegend = 0;          # legend for data points
$ititle = 0;           # title (see below)

# iplot options will override some of the above specifications
#if ($iplot==0) {$plot_ref_points = 1; $plot_ref_labels = 1;}
#if ($iplot==1) {$plot_ref_points = 0; $plot_ref_labels = 1;}
#if ($iplot==2) {$plot_ref_points = 0; $plot_ref_labels = 0;}
if ($iplot==1) {$ilegend = 1;}

$clune = $white;
#$clune = $lgray;
#$clune = $sienna;
if($ipatch==1) {$clune = $lgray;}

# size of markers
$msize_ref = 8;
$msize_data = 6;

# directories with data files
# pdir2 is for a user's directory of files that are outside of compearth
$pdir = "./dfiles/";
$pdir2 = "/home/carltape/PROJECTS/cmt/figs/bfiles/";

# loop over lune, then vw-rectangle
$kmax = 1;
for ($k = 1; $k <= $kmax; $k++) {

if($k==1) {
  $R = $Rlune; $B = $Blune; $J = $Jlune;
  print CSH "psbasemap $J $R $B -G$clune -K -V $origin > $psfile\n"; # START
} else {
  $R = $Rrect; $B = $Brect; $J = $Jrect;
  print CSH "psbasemap $J $R $B -G$clune -K -O -V -X5 >> $psfile\n";
  $ilegend = 0;
}

# columns 1-2 for lune, columns 3-4 for rectangle
$col1 = 2*$k - 1;
$col2 = 2*$k;

# plot patches
if ($ipatch==1) {
  $fname = "$pdir/sourcetype_patch_01.dat";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname | psxy -G$dgray -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/sourcetype_patch_02.dat";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname | psxy -G255 -J -R -K -O -V >>$psfile\n";
  print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";
}
if ( ($ipatchgcdc > 0 || $iarcgcdc > 0) && $inugcdc > 0 ) {
  $fname1 = sprintf("$pdir/sourcetype_patch_nu0p%2.2i_01.dat",100*$nus[$inugcdc-1]);
  $fname2 = sprintf("$pdir/sourcetype_patch_nu0p%2.2i_02.dat",100*$nus[$inugcdc-1]);
  if($ipatchgcdc > 0) {
  print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | psxy -G$red -J -R -K -O -V >>$psfile\n";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname2 | psxy -G$red -J -R -K -O -V >>$psfile\n";
  print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";
}
  if($iarcgcdc > 0) {
  print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | psxy -W2p,$red -J -R -K -O -V >>$psfile\n";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname2 | psxy -W2p,$red -J -R -K -O -V >>$psfile\n";
}
}

# plot arcs
# 1 deviatoric (equator)
# 2 iso+DC (center longitude)
# 3 lam3=0 bottom of +ISO patch
# 4 lam1=0 top of -ISO patch
# 5 CDC nu=0.25
# 6 lam2=0 (CDC nu=0) between dipoles
# 7 CDC nu=0.36
$lwid = 3;
$fname1 = "$pdir/sourcetype_arc_01.dat";
$fname2 = "$pdir/sourcetype_arc_02.dat";
$fname3 = "$pdir/sourcetype_arc_03.dat";
$fname4 = "$pdir/sourcetype_arc_04.dat";
$fname5 = "$pdir/sourcetype_arc_05.dat";
$fname6 = "$pdir/sourcetype_arc_06.dat";
$fname7 = "$pdir/sourcetype_arc_07.dat";
if ($iplot==1) {
  $W = "-W2p,0,--";
  $W = "-W2p,0";
  if($idev==1)   {print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | psxy $W -W2p,$blue -J -R -K -O -V >>$psfile\n";}
  #print CSH "psxy $fname2 $W -J -R -K -O -V >>$psfile\n";
  if($ilam3==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname3 | psxy $W -J -R -K -O -V >>$psfile\n";}
  if($ilam1==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname4 | psxy $W -J -R -K -O -V >>$psfile\n";}
  if($inup25==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname5 | psxy $W -W2p,$green -J -R -K -O -V >>$psfile\n";}
  if($inup36==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname7 | psxy $W -W2p,$green -J -R -K -O -V >>$psfile\n";}
  if($ilam2==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname6 | psxy $W -J -R -K -O -V >>$psfile\n";}
} else {
  if($ipatch==1) {@cols = ($magenta,$red,$blue,$blue,$black,$blue);}
  else           {@cols = ($magenta,$orange,$black,$black,$blue,$black);}
  if($idev==1)   {print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | psxy -W${lwid}p,$cols[0] -J -R -K -O -V >>$psfile\n";}
  #if($ipatch==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname2 | psxy -W${lwid}p,$cols[1] -J -R -K -O -V >>$psfile\n";}
  if($ilam3==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname3 | psxy -W${lwid}p,$cols[2] -J -R -K -O -V >>$psfile\n";}
  if($ilam1==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname4 | psxy -W${lwid}p,$cols[3] -J -R -K -O -V >>$psfile\n";}
  if($inup25==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname5 | psxy -W${lwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  if($inup36==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname7 | psxy -W${lwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  if($ilam2==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname6 | psxy -W${lwid}p,$cols[5] -J -R -K -O -V >>$psfile\n";}
}
if($iarcgcdc==1 && $iplot==1 && $inup25==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname5 | psxy -W2p,$blue -J -R -K -O -V >>$psfile\n";}

# plot lune reference points and labels
if ($plot_ref_points || $plot_ref_labels) {
  if($k==1) {
     $fname = "$pdir/sourcetype_points_lune.dat";
  } else {
     $fname = "$pdir/sourcetype_points_rect.dat";
  }
  #if ($plot_ref_points) {
  #  $csize = 12;
  #  print CSH "psxy $fname -N -Sc${csize}p -W1p,0/0/0 -G255 -J -R -K -O -V >>$psfile\n";
  #}
  $pinfo = "-N -Sc${msize_ref}p -W1p,0/0/0 -G0";
  $pinfo = "-N -Sp${msize_ref}p -W1p,0";
  $fsize = 12;
  $fontno = 1;
  # full set of reference points
  open(IN,$fname); @plines = <IN>; close(IN);
  $nplot0 = @plines;
  # subset of reference points
  @iiplot = (1..8);  # default points (note: might want 10 and 11 as default)
  if($inup25==1) {@iiplot = (@iiplot,10,11)}
  if($inup36==1) {@iiplot = (@iiplot,16,17)}
  if($inugcdc==1 && ($ipatchgcdc==1 || $iarcgcdc==1)) {@iiplot = (@iiplot,10..15)}
  if($inugcdc==2 && ($ipatchgcdc==1 || $iarcgcdc==1)) {@iiplot = (@iiplot,16..21)}
  if($iiplotDC==1) {@iiplot = (@iiplot,9)}
  #@iiplot = (1..$nplot0);  # all points

  $nplot = @iiplot;
  print "plotting $nplot (out of $nplot0) reference points/labels\n";
  #print "\n @plines \n";
  for ($h = 1; $h <= $nplot; $h++) {
    $i = $iiplot[$h-1];
    ($plon,$plat,$plab,$plab2,$align,$Dx,$Dy) = split(" ",$plines[$i-1]);
    #print "\n--$plon -- $plat-- $plab -- $plab2";
    $D = "-D${Dx}p/${Dy}p";
    if (${plot_ref_points}) {
      print CSH "psxy $pinfo -J -R -K -O -V >>$psfile<<EOF\n$plon $plat\nEOF\n";
    }
    if (${plot_ref_labels} > 0) {
      if(${plot_ref_labels}==1) {$ptext=$plab;} else {$ptext=$plab2;}
      print CSH "pstext -N -J -R -K -O -V $D >>$psfile<<EOF\n$plon $plat $fsize 0 $fontno $align $ptext\nEOF\n";
    }
  }
  #print CSH "awk '{print \$1,\$2,$fsize,0,0,\"CM\",\$3}' $fname | pstext -N -J -R -K -O -V >> $psfile\n";

# # plot some test points 
# if(0==1) {
# #beta = [1.571 1.437 1.306  1.178 1.056 0.944 0.845 0.767 0.716];
# #gamma =  [0 0.113 0.23   0.355 0.495 0.654 0.839 1.056 1.303];
# print CSH "psxy -J -R -N $pinfo -K -O -V >>$psfile<<EOF
#  0 0
# -2.747 9.619
# -5.657 19.215
# -8.93 28.759
# -12.864 38.208
# -17.971 47.487
# -25.239 56.443
# -36.788 64.718
# -57.062 71.375 
# -90. 74.207
# }
# EOF\n";
# }

}

if ($iplot==1) {
  # moment tensors from various studies

  # NOTE: NOT ALL OF THESE DATA SETS ARE AVAILABLE HERE;
  #       ONLY THE ONES THAT ARE LISTED IN PUBLICATIONS ARE AVAILABLE (see ./dfiles/README/).
  @ftags = ("Ford2009","Ford2009nuclear","Ford2009earthquake","Ford2009mine","Foulger2004","Minson2007","Minson2008","Walter2009","Walter2010","Pesicek2012",
      "Miller1996phd","Baig2010","Sileny2006","Sileny2008","Sileny2009","Dreger2012","Julian2010","Pesicek2012_238Fig14","Ross1996","Ross1996phd",
      "Vavrycuk2001","Vavrycuk2011","Ortega2014","Alvizuri2015");

  $csize = $msize_data;  # size of dots
  @csizes = ($csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize,
      $csize,$csize/2,$csize,$csize,$csize,$csize/2,$csize,$csize,$csize,$csize,
      $csize,$csize,$csize,$csize);
  @cols = ($lgray,$lred,$lpurple,$green,$orange,$cyan,$blue,$green,$magenta,$lgray,
       $green,$green,$red,$orange,$magenta,$cyan,$cyan,$cyan,$cyan,$dgray,
       $magenta,$brown,$cyan,$red);

  @inds = (9,8,1,6,5);            # TapeTape2012 figure 25
  #@inds = (11,5,22,21,7,6,10,24); # TapeTape2013 figure 14b: mostly volcanic and geothermal
  #@inds = (2..4);                 # TapeTape2013 figure 14c: Ford2009 (nu = 0.25)
  #@inds = (8,9);                  # TapeTape2013 figure 14d: Walter2009,2010 (nu = 0.36)
  #@inds = (12,20,17,13,14,15,4);  # TapeTape2013 figure S14: induced events (lam2 = 0)

  #@inds = (17,13,14,15,4);        # induced -- no Baig
  #@inds = (11,20);                # Foulger2004, Figure 8ab (not c)
  #@inds = (16);                   # Dreger2012 (excluding Long Valley and Geysers regions)
  #@inds = (13..15);               # Sileny
  #@inds = 14;

        for ($i = 1; $i <= @inds; $i++) {
          $j = $inds[$i-1];
          $cz = $csizes[$j-1];
          $fname = sprintf("$pdir/sourcetype_gdvw_%s.dat",$ftags[$j-1]);
          if (not -f $fname) {
              print "input file $fname does not exist\n";
              # check files in second directory
              $fname = sprintf("$pdir2/sourcetype_gdvw_%s.dat",$ftags[$j-1]);
              if (not -f $fname) {die("\n check if input file $fname exists\n")}
          }
          #print CSH "psxy $fname -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] -J -R -K -O -V >>$psfile\n";
          print CSH "awk '{print \$${col1},\$${col2}}' $fname | psxy -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] -J -R -K -O -V >> $psfile\n";
        }

#       # Minson vs GCMT (comment out default block above)
#       @ftags = ("Minson2007 (n=14)","GCMT (n=14)");
#       @cols = ($red,$cyan);
#       @inds = (1,2);
#       $fname1 = "/home/carltape/papers/SOURCE/DATA/MinsonGCMT_Minson.dat";
#       $fname2 = "/home/carltape/papers/SOURCE/DATA/MinsonGCMT_GCMT.dat";
#       print CSH "awk '{print \$8,\$9}' $fname1 | psxy -N -Sc${csize}p -W0.5p,0/0/0 -G$cols[0] -J -R -K -O -V >> $psfile\n";
#       print CSH "awk '{print \$8,\$9}' $fname2 | psxy -N -Sc${csize}p -W0.5p,0/0/0 -G$cols[1] -J -R -K -O -V >> $psfile\n";

      # # Dreger et al. 2012, colored by F-test significance (comment out default block above)
      # $cptfile = "color.cpt";
      # #print CSH "makecpt -Crainbow -T50/100/5 -D > $cptfile\n";
      # print CSH "makecpt -Cseis -T50/100/5 -D -I > $cptfile\n";
      # $ilegend = 0;
      # $fname = "/home/carltape/papers/SOURCE/DATA/bsldreger_fmt_lune.dat";
      # $Fmin = 40;  # try 40,70,90
      # `awk '\$3 > $Fmin' $fname > dtemp`;
      # $nplot = `wc dtemp | awk '{print \$1}'`; chomp($nplot);
      # print "\n$nplot FMTs with F > $Fmin\n";
      # print CSH "awk '{print \$1,\$2,\$3}' dtemp | psxy -N -Sc${csize}p -W0.5p,0/0/0 -C$cptfile -J -R -K -O -V >> $psfile\n";
      # $Dscale = "-D0/1/2/0.2";
      # $Bscale = "-B10f5:\"F-test significance\": -Eb10p";
      # print CSH "psscale -C$cptfile $Dscale $Bscale -Xa3.5 -Ya6 -V -K -O >> $psfile\n";

} elsif ($iplot==2) {
  if($k==2) {die("beachball plotting not yet implemented for rectangle plot -- set kmax = 1\n")}
  # reference beachballs on the lune
  $cmtinfo = "-Sm0.45 -L0.5p/0/0/0 -G255/0/0 -N";
  $cmtfile = sprintf("$pdir/beachballs_ipts%i_iref%i_lune_psmeca",$lplot,$kplot);
  print CSH "psmeca $cmtfile $J $R $cmtinfo -K -O -V >> $psfile\n";

} elsif ($iplot==3) {
  if($k==2) {die("beachball plotting not yet implemented for rectangle plot -- set kmax = 1\n")}
  # catalog of moment tensors on the lune
  $cmtfile = "/home/alvizuri/REPOSITORIES/manuscripts/alvizuri/papers/2014fmt/data/utuhalf_P01_V10_R01_S10_lune_fixedmag_psmeca";
  if (not -f $cmtfile) {die("\n check if cmt file $cmtfile exists\n")}
  $cptfile = "color.cpt";
  print CSH "makecpt -Chaxby -T0/3/0.5 -D > $cptfile\n";
  $cmtinfo = "-Sm0.2 -L0.5p/0/0/0 -N";
  print CSH "psmeca $cmtfile $J $R $cmtinfo -Z$cptfile -K -O -V >> $psfile\n";
  $Dscale = "-D3/7/1.5/0.25";
  $Bscale = "-B1f0.5g0.5:\"Mw\": -Al";
  print CSH "gmtset TICK_LENGTH 8p\n";
  print CSH "psscale -C$cptfile $Dscale $Bscale -V -K -O >> $psfile\n";
} 

#-----------------------------

$J_title = "-JX1i";  # -JM7i
$R_title = "-R0/1/0/1";
$olegend = "-Xa-2.0 -Ya6.5";

# legend for plotting published studies
if($iplot==1 && $ilegend==1) {
  $x0 = 0; $y0 = 1.2; $dy = 0.3;
  for ($i = 1; $i <= @inds; $i++) {
    $j = $inds[$i-1];
    $cz = $csizes[$j-1];
    $x = $x0 + 0.2;
    $y = $y0 - ($i-1)*$dy;
    # number of points (for legend)
    $fname = sprintf("$pdir/sourcetype_gdvw_%s.dat",$ftags[$j-1]);
    if (not -f $fname) {
       $fname = sprintf("$pdir2/sourcetype_gdvw_%s.dat",$ftags[$j-1]);
       if (not -f $fname) {die("\n check if input file $fname exists\n")}
    }
    $nplot = `wc $fname | awk '{print \$1}'`; chomp($nplot);
    $lab = sprintf("%s (n=%i)",$ftags[$j-1],$nplot);
    #$lab = $ftags[$j-1];      # Minson vs GCMT
    print CSH "psxy -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] $R_title $J_title $olegend -K -O -V >>$psfile<<EOF\n$x0 $y\nEOF\n";
    print CSH "pstext -N $R_title $J_title $olegend -K -O -V >>$psfile<<EOF\n $x $y 12 0 $fontno LM $lab\nEOF\n";
  }

#$x = -30; $y = rad2deg(asin(1/sqrt(3)));
#print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G255/165/0 -J -R -K -O -V >>$psfile<<EOF\n$x $y\nEOF\n";
}

#-----------------------------

# final command must have no -K
if($k==$kmax) {$xtag = ""} else {$xtag = "-K"}

# optional: plot a title
$otitle1 = "-Xa2 -Ya9.0"; $fsize1 = 16;
$otitle2 = "-Xa2 -Ya8.7"; $fsize2 = 12;
if ($ititle==1) {
  $title1 = "Representation of source types on the fundamental lune";
  $title2 = "(W. Tape and C. Tape, 2012, GJI, \"A geometric setting for moment tensors\")";
  #$title1 = "Berkeley Seismological Laboratory full moment tensor catalog (n = $nplot, Fsig > $Fmin)";
  #$title2 = "Dreger, Chiang, Ford, Walter, 2012, Monitoring Research Review";
  #$title1 = "Moment tensors for non-induced events"; $title2 = "";
  #$title1 = "Moment tensors for induced events"; $title2 = "";

  print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0 0 $fsize1 0 $fontno CM $title1\nEOF\n";
  print CSH "pstext -N $R_title $J_title $otitle2 $xtag -O -V >>$psfile<<EOF\n 0 0 $fsize2 0 $fontno CM $title2\nEOF\n";

} else {
  # any command with no -K should work
  print CSH "psxy -R -J $xtag -O -T -V >> $psfile\n";
}

}  # for loop over $k

close (CSH);
system("csh -f $cshfile");

# create a pdf file
#system("ps2pdf $psfile");

# you may need to install gv to view (or use something else)
system("gv $psfile &");

# to make composite pdf file:
# for file in `ls lune_hammer_iplot2_*.ps` ; do ps2pdf $file ; done ; pdcat -r lune_hammer_iplot2*pdf all_lune_hammer_iplot2.pdf
# pdcat -r lune_hammer_iplot0_kplot1.pdf lune_hammer_iplot1_kplot1.pdf all_lune_hammer_iplot2.pdf all_lune.pdf

#==================================================
