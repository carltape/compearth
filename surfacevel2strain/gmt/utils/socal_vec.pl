#!/usr/bin/perl -w

#==========================================================
#
#  socal_vec.pl
#  Carl Tape
#  23-Oct-2007
#  
#  This script plots a multiscale residual horizontal velocity field.
#  
#==========================================================

$cshfile = "socal_vec.csh";

$icolor = 0;

# plates and faults
$plate_dir = "/home/carltape/gmt/plates";
$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
if (not -f ${plate_file}) { die("Check if ${plate_file} exist or not\n") }
$kcf_file     = "/home/carltape/gmt/faults/kcf.xy";
$fault_file   = "/home/carltape/gmt/faults/jennings.xy";
if (not -f $fault_file) { die("Check if $fault_file exist or not\n") }

$fault_info_k = "-M -W1.5p,0/0/0";
$fault_info_r = "-M -W1.5p,255/0/0";

$sinfo = "-Sc4p -G255 -W0.5p,0/0/0";
$poleinfo1 = "-Sa16p -G0/255/255 -W0.5p,0/0/0";

# velocity, strain rate, etc
#$vel_dir0   = "${plate_dir}/surface_velocities";
$vel_dir0   = "/home/carltape/SURFACEVEL2STRAIN/matlab_output";

#----------------------------------

# file name, ticks, range for velocities

$iregion = 3;      # region (1=west_us, 2=cal, 3=socal, 6=cascadia)
$idata = 1;        # choose GPS dataset (1 = NASA REASON; 2 = CCMM; 1X = strike-slip; 2X = rotational)
$ndim = 2;         # ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 3;         # min grid order used
$qmax = 8;         # max grid order used
$basis = 1;        # 1 for wavelets; 2 for splines
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$iunrotate = 1;    # 1 to remove rotation; 0 to leave original field

$imask  = 0;       # plot mask
$ieuler = 0;       # plot euler poles
$ifault = 0;       # plot faults
$iplate = 0;       # plot plate boundaries

if ($iregion == 1 ) {
   $itag = "west_us";
   $z_title = 1.9;
   $xtick1 = 2; $xtick2 = 0.25; $wid = 4.5;
   $origin = "-X1.5 -Y1.75";
   $xlegend = 0.25; $ylegend = -1.0;
   $ifault = 0;

} elsif ($iregion == 2) {
   $itag = "cal";
   $z_title = 1.45;
   $xtick1 = 1; $xtick2 = 0.5; $wid = 6.0;
   $origin = "-X1.0 -Y1.75";
   $xlegend = 0.25; $ylegend = -1.0;
   $ifault = 1;

} elsif ($iregion == 3) {
   $itag = "socal";
   $z_title = 1.2;
   $xtick1 = 1; $xtick2 = 0.25; $wid = 6.0;
   #$origin = "-X1.0 -Y1.0";
   #$xlegend = 0.25; $ylegend = 1.3;
   $origin = "-X1.0 -Y2.5";
   $xlegend = 0.25; $ylegend = -1.0;
   $ifault = 1;

} elsif ($iregion == 6) {
   $itag = "cascadia";
   $z_title = 3.55;
   $xtick1 = 2; $xtick2 = 0.5; $wid = 2.5;
   $origin = "-X2.0 -Y1.5";
   $xlegend = 0.25; $ylegend = -1.0;
   $ifault = 1;

}
$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";
$evec_ref = 12;
$vec_conf = 0.95;
$vec_ellipse_frac = 0.1;

# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
$dlab = sprintf("d%2.2i",$idata);
if ($idata==1) {
  $tag = "NASA REASON data set (cGPS)";
  $gps_pts = "${gps_dir}/US/reason_subset_points.dat";
  $gps_vec = "${gps_dir}/US/reason_subset_psvelo.dat";
  $igc = 0; $vec_ref = 50; $vscale = 0.01;
  $cmin = 0; $cmax = 40; $ctick = 10;
  $vumin = -2; $vumax = 2; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;

} elsif ($idata==2) {
  $tag = "California Crustal Motion Map, v1.0  (Shen et al., SCEC-2006)";
  $gps_pts = "${gps_dir}/US/california/socal_vfield_4p0_points.dat";
  $gps_vec = "${gps_dir}/US/california/socal_vfield_4p0_psvelo.dat";
  $igc = 0; $vec_ref = 50; $vscale = 0.01;
  $cmin = 0; $cmax = 40; $ctick = 10;
  $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
}

# synthetic velocity field
if ($idata >= 10) {

  $gps_pts = "${gps_dir}/synthetic/syn_vfield_${itag}_${dlab}_points.dat";
  $gps_vec = "${gps_dir}/synthetic/syn_vfield_${itag}_${dlab}_psvelo.dat"; 

  # strike-slip field
  if ($idata >= 10 and $idata < 20) {
    $igc = 1; $vec_ref = 50; $vscale = 0.01;
    $cmin = 0; $cmax = 40; $ctick = 10;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
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

  # rotational field
  if ($idata >= 20 and $idata < 30) {
    $igc = 0; $vec_ref = 10; $vscale = 0.06;
    $cmin = 0; $cmax = 8; $ctick = 2;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 0.5;
    if ($idata == 20) {
      $tag = "synthetic rotational velocity field (uniform points, no errors)"; $imask = 0;
    } elsif ($idata == 21) {
      $tag = "synthetic rotational velocity field (uniform points, with errors)"; $imask = 0;
    } elsif ($idata == 22) {
      $tag = "synthetic rotational velocity field (REASON points, no errors)";
    } elsif ($idata == 23) {
      $tag = "synthetic rotational velocity field (REASON points, with errors)";
    }
  }

  # 3D strike-slip coseismic field
  if ($idata >= 30 and $idata < 40) {
    $igc = 1; $vec_ref = 1000; $vscale = 0.001;
    $cmin = 0; $cmax = 500; $ctick = 100;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 0.5;
    if ($idata == 30) {
      $tag = "synthetic strike-slip coseismic field (uniform points, no errors)"; $imask = 0;
    } elsif ($idata == 31) {
      $tag = "synthetic strike-slip coseismic field (uniform points, with errors)"; $imask = 0;
    } elsif ($idata == 32) {
      $tag = "synthetic strike-slip coseismic field (REASON points, no errors)";
    } elsif ($idata == 33) {
      $tag = "synthetic strike-slip coseismic field (REASON points, with errors)";
    }
  }

  # 3D thrust coseismic field
  if ($idata >= 50 and $idata < 60) {
    $igc = 1; $vec_ref = 1000; $vscale = 0.0005;
    $cmin = 0; $cmax = 3000; $ctick = 1000;
    $vumin = -2000; $vumax = 2000; $vsmin = -500; $vsmax = 500; $vemin = -3000; $vemax = 3000;
    $evec_ref = 100;
    if ($idata == 50) {
      $tag = "synthetic thrust coseismic field (uniform points, no errors)"; $imask = 0;
    } elsif ($idata == 51) {
      $tag = "synthetic thrust coseismic field (uniform points, with errors)"; $imask = 0;
    } elsif ($idata == 52) {
      $tag = "synthetic thrust coseismic field (REASON points, no errors)";
    } elsif ($idata == 53) {
      $tag = "synthetic thrust coseismic field (REASON points, with errors)";
    }
  }
}
if (not -f $gps_pts) { die("Check if $gps_pts exist or not\n") }
if (not -f $gps_vec) { die("Check if $gps_vec exist or not\n") }

# KEY COMMANDS
$vscale = 0.025;
$vec_ref = 5;
$vec_error = 0.5;
#$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);
$vscale_text = sprintf("%.0f mm/yr", $vec_ref);

#--------------------------------------------

# plotting specifications
$fsize0 = "20";
$fsize1 = "18";
$fsize2 = "14";
$fsize3 = "6";
$fontno = "1";    # 1 or 4
$tick   = "0.2c";
$fpen   = "1p";
$tpen   = "1p";

# plotting specificiations

# A : smallest feature plotted, in km^2; D : resolution
$coast_res     = "-A500 -Df";
$coast_infoK   = "$coast_res -W1.5p,0/0/0 -Na/1.0p";
$coast_infoW   = "$coast_res -W1.0p,255/255/255 -Na/1.0p,255/255/255,t";
$coast_infoR   = "$coast_res -W1.0p,255/0/0 -Na/1.0p,255/0/0,t";

$plate_infoG     = "-M -W1p,0/255/0";
$plate_infoR     = "-M -W1p,255/0/0";
$plate_infoK     = "-M -W1p,0/0/0";
$plate_infoW     = "-M -W1p,255/255/255";
$plate_infoGr    = "-M -W1p,200";

$textinfo = "-G255 -S1p";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

# BOOLEAN: WHICH FIGURES TO PLOT
$ixv    = 1;
$ipdf   = 0;

#==================================================

$glab = sprintf("%2.2i",$igc);
$qlab = sprintf("q%2.2i_q%2.2i",$qmin,$qmax);
$blab = sprintf("b%1.1i",$basis);
$nlab = sprintf("%1iD",$ndim);
$slab = sprintf("s%1i",$iLmat);
$ulab = sprintf("u%1i",$iunrotate);

$name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}_${ulab}";

print "\n $name\n";

if($igc==1){
  $gc_boundary  = "/home/carltape/gmt/gps_data/synthetic/socal_gps_gc_${dlab}.dat";
  #$gc_boundary  = "/home/carltape/gmt/gps_data/synthetic/socal_gps_SAFplanar_gc$glab.dat";
  if (not -f $gc_boundary) { die("Check if $gc_boundary exist or not\n") }
}

# components and magnitudes of v-field
$strain_mag        = "${vel_dir0}/${name}_strain.dat";
$velocity_vec      = "${vel_dir0}/${name}_vec_horz.dat";
$velocity_vec_dat  = "${vel_dir0}/${name}_vec_horz_dat.dat";
$bounds            = "${vel_dir0}/${name}_bounds.dat";
$colors            = "${vel_dir0}/${name}_colors.dat";
$spline_centers    = "${vel_dir0}/${name}_gridpoints.dat";
$euler_poles       = "${vel_dir0}/${name}_euler_vector_psxy.dat";
$euler_anti_poles  = "${vel_dir0}/${name}_euler_anti_vector_psxy.dat";
$euler_poles_scale = "${vel_dir0}/${name}_euler_vector_scale_psxy.dat";

if (not -f $strain_mag)     { die("Check if $strain_mag exist or not\n") }
if (not -f $velocity_vec)   { die("Check if ${velocity_vec} exist or not\n") }
if (not -f $velocity_vec_dat)   { die("Check if ${velocity_vec_dat} exist or not\n") }
if (not -f $bounds)         { die("Check if $bounds exist or not\n") }
if (not -f $colors)         { die("Check if $colors exist or not\n") }
if (not -f ${spline_centers}) { die("Check if ${spline_centers} exist or not\n") }
if (not -f ${euler_poles}) { die("Check if ${euler_poles} exist or not\n") }
if (not -f ${euler_anti_poles}) { die("Check if ${euler_anti_poles} exist or not\n") }
#if (not -f ${euler_poles_scale}) { die("Check if ${euler_pole_scale} exist or not\n") }

# get bounds of the region
open(IN,"$bounds"); @lines = <IN>;
($lonmin,$lonmax,$latmin,$latmax) = split(" ",$lines[0]);

# modify the bounds
if($iregion == 3) {
  $latmin = 32; $latmax = 37; $lonmin = -122; $lonmax = -114;
}

$R1 = "-R$lonmin/$lonmax/$latmin/$latmax";
print "\n bounds are $R1\n\n";

# plot title
$R_title = "-R0/1/0/1";
$x_title = 0.5;

$B0 = "-Ba${xtick1}f${xtick2}d::";

$pname = $name;

#==============================
# from here on out, use portrait plotting with 2-column figures

@flab = ("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p");
$fault_info_k    = "-M -W1.0p,0/0/0";
$fault_info_r    = "-M -W1.0p,255/0/0";

if ($iregion == 3) {
  $wid = 3;
  $origin = "-X1.0 -Y7";
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $xfac = 1.20; $yfac = 1.0;
  $dX1 = $xfac*$wid; $dY1 = 0; 
  $xlegend = 0.25; $ylegend = -1.0;
  $x_stitle = 0.05; $z_stitle = 0.83; $fsize_stitle = 12;  # subtitle

}
$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";

$hwid = $wid/2;
$fsize_title = 10;
$Dscale = "-D$hwid/-0.2/$Dlen/0.10h";
$J1 = "-JM${wid}i";
$J_title = "-JX${wid}";

$dX2 = -$dX1; $dY2 = -$yfac*$wid;
$shift1 = "-X$dX1 -Y$dY1";
$shift2 = "-X$dX2 -Y$dY2";

# open CSH file
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH 0.15c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";

#==============================

  # residuals of velocity field
  $v_file = "${vel_dir0}/${name}_vfield_residual.dat";
  if (not -f $v_file) {die("Check if $vu_file exist or not\n")}
  print "VELOCITY FILE: $v_file\n";

  # additional multiscale files
  $iqs    = "${vel_dir0}/${name}_iqvec.dat";
  $iqs2   = "${vel_dir0}/${name}_iqvec2.dat";
  if (not -f $iqs) {die("Check if $iqs exist or not\n")}
  if (not -f $iqs2) {die("Check if $iqs2 exist or not\n")}

  # get q indexes for the plots
  open(IN,"$iqs"); @qlines1 = <IN>; $nump = @qlines1; print "\n @qlines1";
  open(IN,"$iqs2"); @qlines2 = <IN>; print "\n @qlines2";
  for ($i = 0; $i < $nump; $i = $i+1) {
    ($q1,$q2) = split(" ",$qlines1[$i]);
    $qtags1a[$i] = "q${q1}_q${q2}";
    $qtags1b[$i] = "q=${q1}-${q2}";
  }
  for ($i = 0; $i < $nump-1; $i = $i+1) {
    ($q1,$q2) = split(" ",$qlines2[$i]);
    $qtags2a[$i] = "q${q1}_q${q2}";
    $qtags2b[$i] = "q=${q1}-${q2}";
  }
  print "\n @qtags1a\n @qtags1b\n @qtags2a\n @qtags2b\n";;

#---------

$fname = "${pname}_multi_vec";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = $tag;

# open CSH file
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE 9 ANOT_FONT_SIZE 9 HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 9 FRAME_PEN $fpen TICK_PEN $tpen\n";

@shifts = ($origin,$shift1,$shift2,$shift1);
@Blabs = (12,7,12,7);

#$nump = 2;

for ($i = 1; $i < $nump; $i = $i+1) {

# column indices for velocity field
$vucol = 5 + 3*$i-2;
$vscol = 5 + 3*$i-1;
$vecol = 5 + 3*$i;

$Btag = $Bopts[$Blabs[$i-1]];
$B = "-Ba${xtick1}f${xtick2}d::".$Btag;

if($i==1) {
   print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
 } else {
   print CSH "psbasemap $J1 $R1 $B -K -O -V -P @shifts[$i-1] >> $psfile\n";
}

#plot coast and plate boundary
print CSH "pscoast $J1 $R1 $B $coast_infoR -K -O -V -P >> $psfile\n";
print CSH "psxy $J1 $R1 ${plate_file} $plate_infoR -K -O -V -P >> $psfile\n";

# plot v-field
$vec_info = "-A0.5p/4p/1.5p -Se$vscale/${vec_conf}/0 -W0.5";  # -N to go outside boundary
print CSH "awk '{print \$1,\$2,\$$vecol,-\$$vscol,\$5,\$4}' $v_file | psvelo $J1 $R1 ${vec_info} -K -O -V -P >> $psfile\n";

# plot v-field scale
$origin_arrow = "-Xa0.2 -Ya0.2";
$ytxt = $latmin + 0.35;
print CSH "pstext $J1 $R1 -N -K -O -V -P $origin_arrow >>$psfile<<EOF\n $lonmin $ytxt 10 0 $fontno LM ${vscale_text}\nEOF\n";
print CSH "psvelo $J1 $R1 -N $vec_info $origin_arrow -K -O -V -P >>$psfile<<EOF\n $lonmin $latmin $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

# subtitle
$stitle = "($flab[$i-1])  Data - Est($qtags2b[$i-1])";
print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle\nEOF\n";

}  # end for loop

print CSH "pstext -N $J_title $R_title -O -V -P >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title\nEOF\n"; # FINISH

print CSH "echo convert $psfile to $jpgfile\n";

print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

#==================================================

close (CSH);
system("csh -f $cshfile");

#==================================================
