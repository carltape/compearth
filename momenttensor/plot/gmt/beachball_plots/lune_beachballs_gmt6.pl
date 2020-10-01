#!/usr/bin/perl -w

# This is a custom script adapted from lune.pl on 2019-04-19
# It is intended for plotting a reference set of moment tensors
# and to make sure that the plotting capabilities are correct.
# An abbreviated version of the shell script is lune_beachballs_min.csh
# Please use lune.pl for research purposes.
#
# EXAMPLE A: lune_beachballs_gmt6.pl gmt611
# EXAMPLE B: uncomment the lines with modx below, then:
#            lune_beachballs_gmt6.pl gmt611_psmeca_flags
#

($gmttag) = @ARGV;
if (@ARGV ne 1) {error("lune_beachballs_gmt6.pl gmttag")}

$Rlune = "-R-30/30/-90/90";
#$Rlune = "-R-30/30/-90/90r";
$xtick1 = 10; $ytick1 = 10;
$xtick2 = 5; $ytick2 = 5;
#$Blune = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":WesN";
$Blune = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":wesn";
#$Blune = "-Ba${xtick1}f${xtick2}:\" \":/a${ytick1}f${xtick2}:\" \":wesn";   # no gridlines

$Blune = "-Bxa${xtick1}f${xtick2}g${xtick1} -Bya${ytick1}f${ytick2}g${ytick1} -Bwesn";  # 6.0.1

$wid = 2.8;  # KEY: width of subplot, in inches
$Jlune = "-JH0/${wid}i"; $title1 = "Hammer equal-area ($Jlune)"; $ftag = "hammer";

# UNCOMMENT THESE TWO LINES TO GET A RECTANGLE INSTEAD OF A LUNE
#$rwid = $wid*0.8; $rhgt = $wid*2.8;
#$Jlune = "-JX${rwid}i/${rhgt}i"; $title1 = "Non-geographical ($Jlune)"; $ftag = "ngeo";

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
$iplot = 2;  # =0 (reference lune)
             # =1 (dots from published studies)
             # =2 (reference beachballs)
             # =3 (catalog of full moment tensors)
$lplot = 1;  # =1-4: reference MTs on the lune (iplot=2 only)
$kplot = 1;  # =1-5: orientation of MT at center of lune (iplot=2 only)
$splot = 0;  # =0-5: text labels above reference beachballs (iplot=2 only, lplot=3 only)
@splotlabs = ("","_lam","_gammadelta","_alphanu","_zetaphi","_vw");
$slabel = $splotlabs[$splot];
if($iplot==2) {
  $psfile = "lune_${ftag}_iplot${iplot}_lplot${lplot}_kplot${kplot}_${gmttag}.ps";
} else {
  $psfile = "lune_${ftag}_iplot${iplot}_${gmttag}.ps";
}
# points
$plot_ref_points = 0;  # plot reference points on lune (ISO, DC, etc)
$plot_ref_labels = 0;  # reference labels: =0 (none), =1 (eigs), =2 (ISO, DC, etc)
$iiplotDC = 1;         # plot reference point (and label) at DC
# patches for lam_i = 0 regions
$ipatch = 1;
# arcs
$idev = 1;             # deviatoric arc (delta = 0)
$iiso = 1;             # DC+ISO arc (gamma = 0)
$inup25 = 0;           # nu=0.25 arc between crack points
$inup36 = 0;           # nu=0.36 arc between crack points
$ilam3 = 1;            # lam3 = 0 arc
$ilam2 = 1;            # lam2 = 0 arc between dipoles
$ilam1 = 1;            # lam1 = 0 arc
$arcwid = 3;           # plotting width of arcs (in points)
# gCDC options
$ipatchgcdc = 0;       # patches for gCDC model (inugcdc > 0)
$iarcgcdc = 0;         # boundary arcs for gCDC model (inugcdc > 0)
@nus = (0.25,0.36);    # Poisson values considered for gCDC model
$inugcdc = 0;          #   =1 for nu=0.25, =2 for nu=0.36 (ice)
# other options
$ilegend = 0;          # legend for data points
$ititle = 1;           # title (see below)

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

$fontno = "1";
#$ticklen = "0.2c";
$ticklen = "0.0c";    # no ticks
$tpen = "2p";
$fpen = "2p";

# modifications for plotting beachball text labels (iplot=2)
if($iplot==2 && $splot > 0) {
    $tpen = "0.5p"; $fpen = "0.5p"; $arcwid = 0.5;
    $idev = 0; $iiso = 0; $inup25 = 0; $inup36 = 0;
}
if($iplot==1) {$arcwid = 2;}

# KEY: number of lune plots; paper size
$X0 = 0.5; $Y0 = 0.5; $origin = "-X$X0 -Y$Y0";
$xmax = 6;
#$xmax = 3;    # modx
$pwidth = 9.5;
$dX = $wid + 0.7;
$pheight = $X0 + $dX*$xmax;

$cshfile = "lune_beachballs_${gmttag}.csh";
open(CSH,">$cshfile");
# set GMT default parameters
print CSH "gmt gmtset PS_MEDIA Custom_${pwidth}ix${pheight}i PS_PAGE_ORIENTATION landscape PROJ_LENGTH_UNIT inch MAP_FRAME_TYPE plain MAP_TICK_LENGTH $ticklen MAP_FRAME_PEN $fpen MAP_TICK_PEN $tpen FONT_ANNOT_PRIMARY 12p,Helevetica,black FONT_HEADING 18p,Helvetiva,black FONT_LABEL 10p,Helvetiva,black\n";

# directories with data files
$pdir = "../dfiles/";

$R = $Rlune; $B = $Blune; $J = $Jlune;

# loop over different lune subplots
for ($x = 1; $x <= $xmax; $x++) {

if($x==$xmax) {$kplot=3; $lplot=2;} else {$kplot=$x; $lplot=1;}
#$kplot = 1; $lplot = 1;  # modx

if($x==1) {
  print CSH "gmt psbasemap $J $R ${B}+g$clune -K -V $origin > $psfile\n"; # START
} else {
  print CSH "gmt psbasemap $J $R ${B}+g$clune -K -O -V -X$dX >> $psfile\n";
}

# columns 1-2 for lune, columns 3-4 for rectangle
$col1 = 1;
$col2 = 2;

# plot patches
if ($ipatch==1) {
  $fname = "$pdir/sourcetype_patch_01.dat";
  #print CSH "awk '{print \$${col1},\$${col2}}' $fname > temp.txt\n";
  #print CSH "gmt psxy -R -J -O -K temp.txt -L+x0 -G$dgray >> $psfile\n";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname | gmt psxy -G$dgray -L+x0 -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/sourcetype_patch_02.dat";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname | gmt psxy -G255 -L+x0 -J -R -K -O -V >>$psfile\n";
  print CSH "gmt psbasemap $J $R -Bxa10f5g10 -Bya10f5g10 -Bwesn -K -O -V >> $psfile\n";
}
if ( ($ipatchgcdc > 0 || $iarcgcdc > 0) && $inugcdc > 0 ) {
  $fname1 = sprintf("$pdir/sourcetype_patch_nu0p%2.2i_01.dat",100*$nus[$inugcdc-1]);
  $fname2 = sprintf("$pdir/sourcetype_patch_nu0p%2.2i_02.dat",100*$nus[$inugcdc-1]);
  if($ipatchgcdc > 0) {
  print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | gmt psxy -G$red -J -R -K -O -V >>$psfile\n";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname2 | gmt psxy -G$red -J -R -K -O -V >>$psfile\n";
  print CSH "gmt psbasemap $J $R $B -K -O -V >> $psfile\n";
}
  if($iarcgcdc > 0) {
  print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | gmt psxy -W2p,$red -J -R -K -O -V >>$psfile\n";
  print CSH "awk '{print \$${col1},\$${col2}}' $fname2 | gmt psxy -W2p,$red -J -R -K -O -V >>$psfile\n";
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
$fname1 = "$pdir/sourcetype_arc_01.dat";
$fname2 = "$pdir/sourcetype_arc_02.dat";
$fname3 = "$pdir/sourcetype_arc_03.dat";
$fname4 = "$pdir/sourcetype_arc_04.dat";
$fname5 = "$pdir/sourcetype_arc_05.dat";
$fname6 = "$pdir/sourcetype_arc_06.dat";
$fname7 = "$pdir/sourcetype_arc_07.dat";
if ($iplot==1) {
  # default formatting for arcs
  $W = "-W${arcwid}p,0,--";
  #$W = "-W${arcwid}p,0";
  if($idev==1)   {print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | gmt psxy $W -J -R -K -O -V >>$psfile\n";}  # -W${arcwid}p,$blue
  if($iiso==1)   {print CSH "awk '{print \$${col1},\$${col2}}' $fname2 | gmt psxy $W -J -R -K -O -V >>$psfile\n";}  # -W${arcwid}p,$orange
  if($ilam3==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname3 | gmt psxy $W -J -R -K -O -V >>$psfile\n";}
  if($ilam1==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname4 | gmt psxy $W -J -R -K -O -V >>$psfile\n";}
  if($inup25==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname5 | gmt psxy $W -J -R -K -O -V >>$psfile\n";}  # -W${arcwid}p,$green
  if($inup36==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname7 | gmt psxy $W -J -R -K -O -V >>$psfile\n";}  # -W${arcwid}p,$green
  if($ilam2==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname6 | gmt psxy $W -J -R -K -O -V >>$psfile\n";}

} else {
  if($ipatch==1) {@cols = ($magenta,$red,$blue,$blue,$black,$blue);}
  else           {@cols = ($magenta,$orange,$black,$black,$blue,$black);}
  if($idev==1)   {print CSH "awk '{print \$${col1},\$${col2}}' $fname1 | gmt psxy -W${arcwid}p,$cols[0] -J -R -K -O -V >>$psfile\n";}
  if($iiso==1)   {print CSH "awk '{print \$${col1},\$${col2}}' $fname2 | gmt psxy -W${arcwid}p,$cols[1] -J -R -K -O -V >>$psfile\n";}
  if($ilam3==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname3 | gmt psxy -W${arcwid}p,$cols[2] -J -R -K -O -V >>$psfile\n";}
  if($ilam1==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname4 | gmt psxy -W${arcwid}p,$cols[3] -J -R -K -O -V >>$psfile\n";}
  if($inup25==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname5 | gmt psxy -W${arcwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  if($inup36==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname7 | gmt psxy -W${arcwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  if($ilam2==1)  {print CSH "awk '{print \$${col1},\$${col2}}' $fname6 | gmt psxy -W${arcwid}p,$cols[5] -J -R -K -O -V >>$psfile\n";}
}
if($iarcgcdc==1 && $iplot==1 && $inup25==1) {print CSH "awk '{print \$${col1},\$${col2}}' $fname5 | gmt psxy -W2p,$blue -J -R -K -O -V >>$psfile\n";}

# plot lune reference points and labels
if ($plot_ref_points || $plot_ref_labels) {
     $fname = "$pdir/sourcetype_points_lune.dat";
  $pinfo = "-N -Sc${msize_ref}p -W1p,0/0/0 -G0";
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
      print CSH "gmt psxy $pinfo -J -R -K -O -V >>$psfile<<EOF\n$plon $plat\nEOF\n";
    }
    if (${plot_ref_labels} > 0) {
      if(${plot_ref_labels}==1) {$ptext=$plab;} else {$ptext=$plab2;}
      print CSH "gmt pstext -N -J -R -K -O -V $D >>$psfile<<EOF\n$plon $plat $fsize 0 $fontno $align $ptext\nEOF\n";
    }
  }
  #print CSH "awk '{print \$1,\$2,$fsize,0,0,\"CM\",\$3}' $fname | gmt pstext -N -J -R -K -O -V >> $psfile\n";
}
if ($iplot==2) {
  $xtag = "lune";
  # reference beachballs on the lune
  $beachballfontsize = "8p"; #if($splot==1) {$beachballfontsize = "6p";}
  $letter = "m";
  #if($x==2) {$letter = "z"}  # modx
  #if($x==3) {$letter = "d"}  # modx
  $cmtinfo = "-S${letter}0.45/$beachballfontsize -L0.5p,0/0/0 -G255/0/0 -N";
  $cmtfile = sprintf("$pdir/beachballs_ipts%i_iref%i_%s_psmeca%s",$lplot,$kplot,$xtag,$slabel);
  if (not -f $cmtfile) {die("\n check if cmt file $cmtfile exists\n");}
  print CSH "gmt psmeca $cmtfile $J $R $cmtinfo -K -O -V >> $psfile\n";
  if($splot > 0) {
    $textfile = sprintf("$pdir/pstext%s",$slabel);
    if (not -f $textfile) {die("\n check if text file $textfile exists\n");}
    print CSH "gmt pstext $textfile -JX11i/8.5i -R0/11/0/8.5 -K -O -V -Xa3.5i -Ya0i >> $psfile\n";
  }
}

#-----------------------------

$J_title = "-JX1i";  # -JM7i
$R_title = "-R0/1/0/1";

# optional: plot a title
$tya = $pwidth - 0.9; $tyb = $tya - 0.3;
$otitle1 = "-Xa0.5 -Ya$tya"; $fsize1 = 16;
$otitle2 = "-Xa0.5 -Ya$tyb"; $fsize2 = 12;
if ($ititle==1 && $x==1) {
  $title1 = "Reference sets of moment tensors (input files available in carltape compearth github repository)";
  #$title1 = "Examples using the psmeca flags -Sm (full), -Sz (deviatoric), and -Sd (double couple)";  # modx
  $title2 = "Plotted using $gmttag";
  print CSH "gmt pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0 0 $fsize1 0 $fontno LM $title1\nEOF\n";
  print CSH "gmt pstext -N $R_title $J_title $otitle2 -K -O -V >>$psfile<<EOF\n 0 0 $fsize2 0 $fontno LM $title2\nEOF\n";
}

}  # for loop over $x

# any command with no -K should work
print CSH "psxy -R -J -O -T -V >> $psfile\n";

close (CSH);
system("csh -f $cshfile");

# create a pdf file
#system("ps2pdf $psfile");

# you may need to install gv to view (or use something else)
system("gv $psfile &");

#==================================================
