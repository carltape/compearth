#!/usr/bin/perl -w

#==========================================================
#
#  socal_rot.pl
#  Carl Tape
#  23-Aug-2007
#  
#  This script inputs a surface velocity field aninfod strain scalar fields
#  computed in test_platemodel2strain.m and outputs a figure.
#  
#==========================================================

$cshfile = "socal_rot.csh";

$icolor = 1;    # ccc

$dir0 = "/home/carltape/gmt";

# plates and faults
$plate_dir = "$dir0/plates";
$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
if (not -f ${plate_file}) { die("Check if ${plate_file} exist or not\n") }
$kcf_file     = "$dir0/faults/kcf.xy";
$fault_file   = "$dir0/faults/jennings.xy";
if (not -f $fault_file) { die("Check if $fault_file exist or not\n") }

$fault_info_k = "-M -W0.5p,0/0/0";
$fault_info_r = "-M -W0.5p,255/0/0";

$sinfo = "-Sc4p -G255 -W0.5p,0/0/0";
$poleinfo1 = "-Sa20p -G0/255/255 -W0.75p,0/0/0";
$circleinfo = "-Sc25p -W1.0p,0/0/0,--";

# velocity, strain rate, etc
$vel_dir0   = "${plate_dir}/surface_velocities";

#----------------------------------

# file name, ticks, range for velocities

$iregion = 3;      # region (1=west_us, 2=cal, 3=socal, 6=cascadia)
#$idata = 83;        # choose GPS dataset (1 = NASA REASON; 2 = CCMM; 1X = strike-slip; 2X = rotational)
#$ndim = 3;         # ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 1;         # min grid order used
$qmax = 7;         # max grid order used
$basis = 2;        # 1 for splines; 2 for wavelets
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$imask  = 1;       # plot mask
$ieuler = 1;       # plot euler poles

@labs = ("(a)  Rotation rate","(b)  Euler vectors","(c)  Rotation rate","(d)  Euler vectors","(e)  Rotation rate","(f)  Euler vectors");
@iBs = (2,2,14,14,2,2);
@iplates = (0,0,0,1,0,1);
@ifaults = (0,0,1,0,1,0);
$origin = "-X0.5 -Y7"; $wid = 2.0; 
$dX = $wid + 0.5; $dY = 2.6;

$t = 0;
#======================================================
for ($k = 1; $k <= 3; $k = $k+1) {

if ($k == 1 ) {
   $idata = 13; $ndim = 2; $imask = 1;
   $shift1 = " "; $shift2 = "-Y-$dY";

} elsif ($k == 2) {
   $idata = 1; $ndim = 2; $imask = 1;
   $shift1 = "-X$dX -Y$dY"; $shift2 = "-Y-$dY";

} elsif ($k == 3) {
   $idata = 1; $ndim = 3; $imask = 1;
   $shift1 = "-X$dX -Y$dY"; $shift2 = "-Y-$dY";
}

$t = $t+1;
$ib = $iBs[$t-1];
$ifault = @ifaults[$t-1];
$iplate = @iplates[$t-1];

$itag = "socal";
$z_title = 1.15;
$xtick1 = 2; $xtick2 = 1;
$vscale = 0.005;
$xlab = 0.95*$wid; $ylab = 1.0*$wid;
$xlegend = 4; $ylegend = -1.0;
$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";

# VELOCITY FIELD DATASET
$gps_dir  = "$dir0/gps_data";
$dlab = sprintf("d%2.2i",$idata);
if ($idata==1) {
  $tag = "NASA REASON data set (cGPS)";
  $gps_pts = "${gps_dir}/US/reason_fixed_NAM_subset_points.dat";
  $gps_vec = "${gps_dir}/US/reason_fixed_NAM_subset_psvelo.dat";
  $igc = 0; $vec_ref = 50; $vec_conf = 0.9; $vscale = 0.01;
  $cmin = 0; $cmax = 40; $ctick = 10;
  $vumin = -2; $vumax = 2; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;

} elsif ($idata==2) {
  $tag = "California Crustal Motion Map, v1.0  (Shen et al., SCEC-2006)";
  $gps_pts = "${gps_dir}/US/california/socal_vfield_4p0_points.dat";
  $gps_vec = "${gps_dir}/US/california/socal_vfield_4p0_psvelo.dat";
  $igc = 0; $vec_ref = 50; $vec_conf = 0.9; $vscale = 0.01;
  $cmin = 0; $cmax = 40; $ctick = 10;
  $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;
}

  # strike-slip field
  if ($idata >= 10 and $idata < 20) {
    $gps_pts = "${gps_dir}/synthetic/syn_vfield_${dlab}_points.dat";

    $igc = 1; $ifault = 0; $iplate = 0;
    $vec_ref = 50; $vec_conf = 0.9; $vscale = 0.01;
    $cmin = 0; $cmax = 40; $ctick = 10;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 2.0; $escale = 5;
    $fac_min = 0.50;   # 1.0 for uniform color scaling
    $fac_max = 1.0;    # 1.0 default
    if ($idata == 10) {
      $tag = "synthetic infinite strike-slip fault (uniform points, no errors)"; $imask = 0;
    } elsif ($idata == 11) {
      $tag = "synthetic infinite strike-slip fault (uniform points, with errors)"; $imask = 0;
    } elsif ($idata == 12) {
      $tag = "synthetic infinite strike-slip fault (REASON points, no errors)";
    } elsif ($idata == 13) {
      $tag = "synthetic infinite strike-slip fault (REASON points, with errors)";
    }
  }
if (not -f $gps_pts) { die("Check if $gps_pts exist or not\n") }

$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);

#--------------------------------------------

# plotting specifications
$fsize0 = "20";
$fsize1 = "18";
$fsize2 = "14";
$fsize3 = "10";
$fontno = "1";    # 1 or 4
$tick   = "0.2c";
$fpen   = "1p";
$tpen   = "1p";

$phi = "\@~\146\@~";
$theta = "\@~\161\@~";

# plotting specificiations

# A : smallest feature plotted, in km^2; D : resolution
$coast_res     = "-A1000 -Dl";
$coast_infoK   = "$coast_res -W1.5p,0/0/0 -Na/1.0p";
$coast_infoW   = "$coast_res -W1.5p,255/255/255 -Na/1p,255/255/255,t";
$coast_infoR   = "$coast_res -W1.5p,255/0/0 -Na/1p,255/0/0,t";

$plate_infoG     = "-M -W1.5p,0/255/0";
$plate_infoR     = "-M -W1.5p,255/0/0";
$plate_infoK     = "-M -W1.5p,0/0/0";
$plate_infoW     = "-M -W1.5p,255/255/255";
$plate_infoGr    = "-M -W3p,200";

$textinfo = "-G255 -S1p";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

#-----------------------------------------------
# BOOLEAN: WHICH FIGURES TO PLOT
$ixv    = 1;
$ipdf   = 0;

#==================================================

$glab = sprintf("%2.2i",$igc);
$qlab = sprintf("q%2.2i_q%2.2i",$qmin,$qmax);
$blab = sprintf("b%1.1i",$basis);
$nlab = sprintf("%1iD",$ndim);
$slab = sprintf("s%1i",$iLmat);

$name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}";

print "\n $name \n";

if($igc==1){
  $gc_boundary  = "$dir0/gps_data/synthetic/gps_gc_${dlab}.dat";
  if (not -f $gc_boundary) { die("Check if $gc_boundary exist or not\n") }
}

# components and magnitudes of v-field
$strain_mag        = "${vel_dir0}/misc/${name}_spline_strain.dat";
$bounds            = "${vel_dir0}/misc/${name}_spline_bounds.dat";
$colors            = "${vel_dir0}/misc/${name}_spline_colors.dat";
$spline_centers    = "${vel_dir0}/misc/${name}_spline_gridpoints.dat";
$euler_poles       = "${vel_dir0}/misc/${name}_euler_vector_psxy.dat";
$euler_anti_poles  = "${vel_dir0}/misc/${name}_euler_anti_vector_psxy.dat";
$euler_poles_scale = "${vel_dir0}/misc/${name}_euler_vector_scale_psxy.dat";

if (not -f $strain_mag)     { die("Check if $strain_mag exist or not\n") }
if (not -f $bounds)         { die("Check if $bounds exist or not\n") }
if (not -f $colors)         { die("Check if $colors exist or not\n") }
if (not -f ${spline_centers}) { die("Check if ${spline_centers} exist or not\n") }
if (not -f ${euler_poles}) { die("Check if ${euler_poles} exist or not\n") }
if (not -f ${euler_anti_poles}) { die("Check if ${euler_anti_poles} exist or not\n") }
#if (not -f ${euler_poles_scale}) { die("Check if ${euler_poles_scale} exist or not\n") }

# get bounds of the region
open(IN,"$bounds"); @lines = <IN>;
($lonmin,$lonmax,$latmin,$latmax) = split(" ",$lines[0]);
$latmin = 32; $latmax = 38; $lonmin = -122; $lonmax = -114;

$R1 = "-R$lonmin/$lonmax/$latmin/$latmax";
print "\n bounds are $R1 \n\n";
#die("testing");

# get color limits for the plots: dilatation, strain, rotation, maxlam
open(IN,"$colors"); @lines = <IN>;
($cminrot,$cmaxrot,$cpwr4) = split(" ",$lines[2]);  # rotation
$cmaxrot = 2.4;
$norm2 = 10**$cpwr2; $norm3 = 10**$cpwr3; $norm4 = 10**$cpwr4; $norm5 = 10**$cpwr5;
print "\n $cpwr2 $cpwr3 $cpwr4 $cpwr5 $norm2 $norm3 $norm4 $norm5 \n\n";

$title_rot   = sprintf("Estimated rotation, 10\@+%1i\@+ yr\@+-1\@+",$cpwr4);
$Bscale_rot  = sprintf("-B%2.2f:\"Rotation rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxrot*0.25,$cpwr4);
$unit_rot    = sprintf("10\@+%i\@+ yr\@+-1\@+",$cpwr4);

$J1 = "-JM${wid}i";

# plot title
$J_title = "-JX${wid}";
$R_title = "-R0/1/0/1";
$x_title = 0.5;

$fsize_title = 16;

# open CSH file
if($k == 1) {
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
}

#-----------------------------------------
# make color files

# rotation
$cran = $cmaxrot - $cminrot;
$dc = $cran/100;
$Trot = "-T$cminrot/$cmaxrot/$dc";

$colorbar = "rainbow";
$c_romania = "$dir0/color_maps/Romanian_flag_smooth.cpt";

# 170     0      0  to 255   255    255

# -I to invert the color scale

$colorbar = "seis";      # -I to flip
$colorbar = $c_romania;
$cptrot = "color4.cpt";
print CSH "makecpt -C$colorbar $Trot -D > $cptrot\n";

#-----------------------------------------

$B0 = "-Ba${xtick1}f${xtick2}d::";

# KEY: commands for interpolation and masking
$interp = "-I0.1";
$interp_surf = "-I0.04 -T0 -S0";
#if($iregion==1) {$interp_surf = "-I0.2 -T0 -S0";}
$mask_info = "-I0.1 -G200 -S0.2";
$grdfile = "temp.grd";
$mask_file = "${vel_dir0}/misc/${name}_masked_pts.dat";
if (not -f ${mask_file}) { die("Check if ${mask_file} exist or not\n") }

$pname = $name;

#-----------------------------------------

#==============================
# from here on out, use portrait plotting with 2-column figures

$fault_info_k    = "-M -W0.5p,0/0/0";
$fault_info_r    = "-M -W0.5p,255/0/0";

$hwid = $wid/2;
$fsize_title = 12;
$Dlen = 1.7;
$Dscale = "-D2/-3/$Dlen/0.15h";
$J1 = "-JM${wid}i";
$J_title = "-JX${wid}";

print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH 0.15c LABEL_FONT_SIZE $fsize3 ANOT_FONT_SIZE $fsize3  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize3 FRAME_PEN 1p TICK_PEN 1p\n";

#==============================

  $fname = "socal_rot";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
  $B = $B0.$Bopts[$ib];

#================================

#$title = sprintf("$labs[$t-1]  Rotation rate, 10\@+%i\@+ yr\@+-1",$cpwr4);
#$Bscale  = sprintf("-B%2.2f:\"  \": -Ef10p",$cmaxrot*0.25);
$title = $labs[$t-1];
$Bscale = $Bscale_rot;

if($k == 1) {
  print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
} else {
  print CSH "psbasemap $J1 $R1 $B -K -O -V $shift1 >> $psfile\n";
}

if($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$6/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1 \n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptrot -T -K -O -V >> $psfile\n";
}
if ($imask==1) {
  print CSH "echo mask file is ${mask_file}\n";
  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile \n";
  print CSH "psmask -C -K -O -V >> $psfile \n";
}
if($k==1) {print CSH "psscale -C$cptrot $Dscale $Bscale -K -O -V >> $psfile \n";}
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}

if ($idata >= 60 and $idata < 80) {
  $epoles = "${vel_dir0}/misc/${name}_epole_mean_points.dat";
  if (not -f $epoles) { die("Check if $epoles exist or not\n") }
  print CSH "psxy $epoles $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile\n";
}

# plot spline centers
#print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

#================================

$t = $t+1;
$ib = $iBs[$t-1];
$ifault = @ifaults[$t-1];
$iplate = @iplates[$t-1];
$title = "$labs[$t-1]";

$R1 = "-R-126/-110/30/42";
#$coast_infoK = "-Dl -W0.5p,0/0/0 -A1000 -S200";
$B0 = "-Ba5f1d::";
$B = $B0.$Bopts[$ib];

print CSH "psbasemap $J1 $R1 $B -K -O -V $shift2 >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";

if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoR -K -O -V >> $psfile\n";}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoR -K -O -V >> $psfile\n";}

# plot euler poles and anti-poles
if ($ieuler == 1) {
  $pinfo      = "-Sc4p -G255/255/255 -W0.25p,0/0/0";
  $pinfo_anti = "-Sc4p -G255/0/0 -W0.25p,0/0/0";
  $seis_info3 = "-Scp -W0.5p/255/255/255";

  # make colorpoint file
  $cptfile = "color.cpt"; $cmin = -0.01; $cmax = 0.01; $dc = ($cmax-$cmin)/2;
  $T = "-T$cmin/$cmax/$dc -D -Z "; $colorbar = "-Cpolar";
  print CSH "makecpt $colorbar $T > $cptfile \n";

  #print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_poles} > poles_test\n"; die("TESTING");

  # KEY: plot euler poles
  #print CSH "psxy ${euler_poles} $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
  #print CSH "psxy ${euler_anti_poles} $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_anti_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";

  #-----------------------------

  if ($k == 1) {

    $xdots = $xlegend;
    $ydots = $ylegend - 0.3;
    $origin_dots = "-Xa${xdots} -Ya${ydots}";

    # plot euler vector scale
    @vec_ref = ($evec_ref*0.75, $evec_ref, 1.25*$evec_ref);
    $sv1 = sprintf("%.2f",$vec_ref[0]);
    $sv2 = sprintf("%.2f",$vec_ref[1]);
    $sv3 = sprintf("%.2f",$vec_ref[2]);
    #open(IN,${euler_poles_scale}); $escale = <IN>; chomp($escale);
    #if (not -f ${euler_poles_scale}) {die("Check if ${euler_poles_scale} exist or not\n");}

    # MAKE SURE THAT THIS FORMULA MATCHES WHAT IS USED IN PLOTTING THE POLES ABOVE
    $dsize[0] = $escale*abs($vec_ref[0]);
    $dsize[1] = $escale*abs($vec_ref[1]);
    $dsize[2] = $escale*abs($vec_ref[2]);
    $dy = 0.25;
    for ($ik = 0; $ik < 5; $ik = $ik+1) {
      $ytxt[$ik] = $ik*$dy;
    }
    $x1 = 0;
    $x2 = $x1 + 1.00;
    print CSH "psxy -JX1 -R0/1/0/1 -C$cptfile -N ${seis_info3} ${origin_dots} -K -O -V >>$psfile<<EOF
$x1 $ytxt[3] $vec_ref[2] $dsize[2]
$x1 $ytxt[2] $vec_ref[1] $dsize[1]
$x1 $ytxt[1] $vec_ref[0] $dsize[0]
$x2 $ytxt[3] -$vec_ref[2] $dsize[2]
$x2 $ytxt[2] -$vec_ref[1] $dsize[1]
$x2 $ytxt[1] -$vec_ref[0] $dsize[0]
EOF\n";

    $x1 = ($x1+$x2)/2;
    $y1 = -0.8*$dy;
    print CSH "pstext -JX -R0/1/0/1 -N -K -O -V ${origin_dots} >>$psfile<<EOF
$x1 $ytxt[4] 11 0 $fontno CM pos              neg
$x1 $ytxt[3] $fsize3 0 $fontno CM $sv3
$x1 $ytxt[2] $fsize3 0 $fontno CM $sv2
$x1 $ytxt[1] $fsize3 0 $fontno CM $sv1
$x1 $ytxt[0] 11 0 $fontno CM $unit_rot
EOF\n";
  }
}

if ($idata >= 60 and $idata < 80) {
  $epoles = "${vel_dir0}/misc/${name}_epole_mean_points.dat";
  if (not -f $epoles) { die("Check if $epoles exist or not\n") }
  print CSH "psxy $epoles $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile\n";
}

      # socal box
print CSH "psxy $J1 $R1 -W1.5p -A -O -K -V <<EOF>>$psfile\n"; 
print CSH "$lonmin $latmin\n"; 
print CSH "$lonmax $latmin\n"; 
print CSH "$lonmax $latmax\n"; 
print CSH "$lonmin $latmax\n"; 
print CSH "$lonmin $latmin\n"; 
print CSH "EOF\n"; 

# plot spline centers
#print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";

}
#  END OF FOR LOOP

  #-----
  print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk \nEOF\n"; # FINISH

  print CSH "echo done with $psfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if ($ixv==1) {print CSH "xv $jpgfile &\n"}


#==================================================

close (CSH);
system("csh -f $cshfile");
#print "convert $psfile $jpgfile \n";
#system("convert $psfile $jpgfile");
#if($ipdf==1) {system("ps2pdf $psfile")};
#if($ixv==1) {system("xv $jpgfile &")};

#==================================================
