#!/usr/bin/perl -w

#==========================================================
#
#  lune.pl
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

#use Math::Trig;

$cshfile = "lune.csh";

$fontno = "1"; 
$ticklen = "0.0c";    # =0 for no ticks

open(CSH,">$cshfile");
print CSH "gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $ticklen LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p\n";

$R = "-R-30/30/-90/90";
#$R = "-R-30/30/-90/90r";
$origin = "-X2 -Y1";
$xtick1 = 10; $ytick1 = 10;
$xtick2 = 5; $ytick2 = 5;
#$B = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":WesN";
$B = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":wesn";
#$B = "-Ba${xtick1}f${xtick2}:\" \":/a${ytick1}f${xtick2}:\" \":wesn";   # no gridlines

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

# PLOTTING OPTIONS
$iplot = 1;  # =0 (reference lune), =1 (dots from published studies), =2 (reference beachballs)
$lplot = 1;  # =1-2: reference MTs on the lune (iplot=2 only)
$kplot = 3;  # =1-4: orientation of MT at center of lune (iplot=2 only)
if($iplot==2) {
  $psfile = "lune_${ftag}_iplot${iplot}_lplot${lplot}_kplot${kplot}.ps";
} else {
  $psfile = "lune_${ftag}_iplot${iplot}.ps";
}
$ipatch = 0;           # two shaded patches on the lune near the ISO regions
$ipatchgcdc = 0;       # patches for gCDC model for nu=0.25
$iarcgcdc = 0;         # boundary arcs for gCDC model for nu=0.25

$idev = 0;             # deviatoric arc
$inup25 = 1;           # nu=0.25 arc between crack points
$inup36 = 0;           # nu=0.36 arc between crack points
$ilam3 = 1;            # lam3 = 0 arc
$ilam2 = 1;            # lam2 = 0 arc between dipoles
$ilam1 = 1;            # lam1 = 0 arc

$ilegend = 0;          # legend for data points
$ititle = 0;           # title (see below)
$plot_ref_points = 1;  # plot reference points on lune (ISO, DC, etc)
$plot_ref_labels = 1;  # reference labels: =0 (none), =1 (eigs), =2 (ISO, DC, etc)

# iplot options will override some of the above specifications
#if ($iplot==0) {$plot_ref_points = 1; $plot_ref_labels = 1;}
#if ($iplot==1) {$plot_ref_points = 0; $plot_ref_labels = 1;}
#if ($iplot==2) {$plot_ref_points = 0; $plot_ref_labels = 0;}
if ($iplot==1) {$ilegend = 1;}

$clune = $white;
#$clune = $lgray;
#$clune = $sienna;
if($ipatch==1) {$clune = $lgray;}

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
if ($ipatchgcdc==1 || $iarcgcdc==1) {
  #$fname1 = "$pdir/beach_patch_nu0p25_01.lonlat";
  #$fname2 = "$pdir/beach_patch_nu0p25_02.lonlat";
  $fname1 = "$pdir/beach_patch_nu0p36_01.lonlat";
  $fname2 = "$pdir/beach_patch_nu0p36_02.lonlat";
  if($ipatchgcdc==1) {
  print CSH "psxy $fname1 -G$red -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname2 -G$red -J -R -K -O -V >>$psfile\n";
  print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";
}
  if($iarcgcdc==1) {
  print CSH "psxy $fname1 -W2p,$red -J -R -K -O -V >>$psfile\n";
  print CSH "psxy $fname2 -W2p,$red -J -R -K -O -V >>$psfile\n";
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
$fname1 = "$pdir/beach_arc_01.lonlat";
$fname2 = "$pdir/beach_arc_02.lonlat";
$fname3 = "$pdir/beach_arc_03.lonlat";
$fname4 = "$pdir/beach_arc_04.lonlat";
$fname5 = "$pdir/beach_arc_05.lonlat";
$fname6 = "$pdir/beach_arc_06.lonlat";
$fname7 = "$pdir/beach_arc_07.lonlat";
if ($iplot==1) {
  $W = "-W2p,0,--";
  $W = "-W2p,0";
  if($idev==1) {print CSH "psxy $fname1 $W -W2p,$blue -J -R -K -O -V >>$psfile\n";}
  #print CSH "psxy $fname2 $W -J -R -K -O -V >>$psfile\n";
  if($ilam3==1) {print CSH "psxy $fname3 $W -J -R -K -O -V >>$psfile\n";}
  if($ilam1==1) {print CSH "psxy $fname4 $W -J -R -K -O -V >>$psfile\n";}
  if($inup25==1) {print CSH "psxy $fname5 $W -W2p,$green -J -R -K -O -V >>$psfile\n";}
  if($inup36==1) {print CSH "psxy $fname7 $W -W2p,$green -J -R -K -O -V >>$psfile\n";}
  if($ilam2==1) {print CSH "psxy $fname6 $W -J -R -K -O -V >>$psfile\n";}
} else {
  if($ipatch==1) {@cols = ($magenta,$red,$blue,$blue,$black,$blue);}
  else           {@cols = ($magenta,$orange,$black,$black,$blue,$black);}
  if($idev==1) {print CSH "psxy $fname1 -W${lwid}p,$cols[0] -J -R -K -O -V >>$psfile\n";}
  if($ipatch==1) {print CSH "psxy $fname2 -W${lwid}p,$cols[1] -J -R -K -O -V >>$psfile\n";}
  if($ilam3==1) {print CSH "psxy $fname3 -W${lwid}p,$cols[2] -J -R -K -O -V >>$psfile\n";}
  if($ilam1==1) {print CSH "psxy $fname4 -W${lwid}p,$cols[3] -J -R -K -O -V >>$psfile\n";}
  if($inup25==1) {print CSH "psxy $fname5 -W${lwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  if($inup36==1) {print CSH "psxy $fname7 -W${lwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  if($ilam2==1) {print CSH "psxy $fname6 -W${lwid}p,$cols[5] -J -R -K -O -V >>$psfile\n";}
}
if($iarcgcdc==1 && $iplot==1 && $inup25==1) {print CSH "psxy $fname5 -W2p,$blue -J -R -K -O -V >>$psfile\n";}

# plot lune reference points and labels
if ($plot_ref_points || $plot_ref_labels) {
  $fname = "$pdir/beach_points.lonlat";
  #if ($plot_ref_points) {
  #  $csize = 12;
  #  print CSH "psxy $fname -N -Sc${csize}p -W1p,0/0/0 -G255 -J -R -K -O -V >>$psfile\n";
  #}
  $csize = 8;
  $pinfo = "-N -Sc${csize}p -W1p,0/0/0 -G0";
  $pinfo = "-N -Sp${csize}p -W1p,0";
  $fsize = 14;
  $fontno = 1;
  open(IN,$fname); @plines = <IN>; close(IN);
  if($ipatchgcdc==1 || $iarcgcdc==1 || $inup36==1) {
      $nplot = @plines;
  } else {
     $nplot=9;
     $nplot=8;  # no DC point
     if($inup25==1) {$nplot=11}
  };
  print "plotting $nplot reference points/labels\n";
  for ($i = 1; $i <= $nplot; $i++) {
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

  # NOTE: NOTE ALL OF THESE DATA SETS ARE AVAILABLE HERE
  @ftags = ("Ford2009","Ford2009nuclear","Ford2009earthquake","Ford2009mine","Foulger2004","Minson2007","Minson2008","Walter2009","Walter2010","Pesicek2012","Miller1996phd",
      "Baig2010","Sileny2006","Sileny2008","Sileny2009","Dreger2012","Julian2010","Pesicek2012_238Fig14","Ross1996","Ross1996phd","Vavrycuk2001","Vavrycuk2011");

  $csize = 8;  # size of dots
  @csizes = ($csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize,
      $csize/2,$csize,$csize,$csize,$csize/2,$csize,$csize,$csize,$csize,$csize,$csize);
  @cols = ($lgray,$lred,$lpurple,$green,$orange,$red,$blue,$green,$magenta,$lgray,$green,
       $green,$red,$orange,$magenta,$cyan,$cyan,$cyan,$cyan,$dgray,$magenta,$brown);

  @inds = (9,8,1,6,5);            # TapeTape2012 figure 25
  #@inds = (11,5,22,21,7,6,10);    # non-induced
  #@inds = (12,20,17,13,14,15,4);  # induced
  #@inds = (17,13,14,15,4);        # induced -- no Baig
  #@inds = (11,20);                # Foulger2004, Figure 8
  #@inds = (16);                   # Dreger2012 (excluding Long Valley and Geysers regions)
  #@inds = (13..15);               # Sileny     
  #@inds = (2..4);                 # Ford2009
  #@inds = (8,9);                   # Walter2009,2010
  #@inds = 14;

        for ($i = 1; $i <= @inds; $i++) {
          $j = $inds[$i-1];
          $cz = $csizes[$j-1];
          $fname = sprintf("$pdir/beachpts_%s_points.dat",$ftags[$j-1]);
          if (not -f $fname) {die("\n check if input file $fname exists\n")}
          print CSH "psxy $fname -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] -J -R -K -O -V >>$psfile\n";
        }

#       # Minson vs GCMT (comment out default block above)
#       @ftags = ("Minson2007 (n=14)","GCMT (n=14)");
#       @cols = ($red,$cyan);
#       @inds = (1,2);
#       $fname1 = "/home/carltape/papers/SOURCE_INVERSION/DATA/MinsonGCMT_Minson.dat";
#       $fname2 = "/home/carltape/papers/SOURCE_INVERSION/DATA/MinsonGCMT_GCMT.dat";
#       print CSH "awk '{print \$8,\$9}' $fname1 | psxy -N -Sc${csize}p -W0.5p,0/0/0 -G$cols[0] -J -R -K -O -V >> $psfile\n";
#       print CSH "awk '{print \$8,\$9}' $fname2 | psxy -N -Sc${csize}p -W0.5p,0/0/0 -G$cols[1] -J -R -K -O -V >> $psfile\n";

#      # Dreger et al. 2012, colored by F-test significance (comment out default block above)
#      $cptfile = "color.cpt";
#      #print CSH "makecpt -Crainbow -T50/100/5 -D > $cptfile\n";
#      print CSH "makecpt -Cseis -T50/100/5 -D -I > $cptfile\n";
#      $ilegend = 0;
#      $fname = "/home/carltape/papers/SOURCE_INVERSION/DATA/bsldreger_fmt_lune.dat";
#      $Fmin = 40;  # try 40,70,90
#      `awk '\$3 > $Fmin' $fname > dtemp`;
#      $nplot = `wc dtemp | awk '{print \$1}'`; chomp($nplot);
#      print "\n$nplot FMTs with F > $Fmin\n";
#      print CSH "awk '{print \$1,\$2,\$3}' dtemp | psxy -N -Sc${csize}p -W0.5p,0/0/0 -C$cptfile -J -R -K -O -V >> $psfile\n";
#      $Dscale = "-D0/1/2/0.2";
#      $Bscale = "-B10f5:\"F-test significance\": -Eb10p";
#      print CSH "psscale -C$cptfile $Dscale $Bscale -Xa3.5 -Ya6 -V -K -O >> $psfile\n";

} elsif ($iplot==2) {
  # reference beachballs on the lune
  $cmtinfo = "-Sm0.5 -L0.5p/0/0/0 -G255/0/0 -N";
  $cmtfile = sprintf("$pdir/beachballs_ilune%i_iref%i_psmeca",$lplot,$kplot);
  print CSH "psmeca $cmtfile $J $R $cmtinfo -K -O -V >> $psfile\n";
} 

#-----------------------------

$J_title = "-JX1i";  # -JM7i
$R_title = "-R0/1/0/1";
$olegend = "-Xa3.0 -Ya7.2";

# legend for plotting published studies
if($iplot==1 && $ilegend==1) {
  $x0 = 0; $y0 = 1.2; $dy = 0.3;
  for ($i = 1; $i <= @inds; $i++) {
    $j = $inds[$i-1];
    $cz = $csizes[$j-1];
    $x = $x0 + 0.2;
    $y = $y0 - ($i-1)*$dy;
    # number of points (for legend)
    $fname = sprintf("$pdir/beachpts_%s_points.dat",$ftags[$j-1]);
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
  print CSH "pstext -N $R_title $J_title $otitle2 -O -V >>$psfile<<EOF\n 0 0 $fsize2 0 $fontno CM $title2\nEOF\n";

} else {
  # any command with no -K should work (this one will not plot anything)
  print CSH "pstext $R_title $J_title $otitle2 -O -V >>$psfile<<EOF\n -1 -1 $fsize2 0 $fontno LM TEST\nEOF\n"; 
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
