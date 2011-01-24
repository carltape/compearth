#!/usr/bin/perl -w

#==========================================================
#
#  plot_strain.pl
#  Carl Tape
#  24-Jan-2011
#  
#  This script plots figures using output files from surfacevel2strain/matlab_output/
#  See Tape et al. (GJI 2009) for examples and details.
#  
#==========================================================

$cshfile = "plot_strain.csh";

$icolor = 1;    # ccc

# plates and faults
#$plate_dir = "/home/carltape/gmt/plates";
#$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
$plate_file  = "input/bird_boundaries";
if (not -f ${plate_file}) { die("Check if ${plate_file} exist or not\n") }
$fault_file = "input/jennings_more.xy";
if (not -f $fault_file) { die("Check if $fault_file exist or not\n") }

$fault_info_k = "-m -W1.5p,0/0/0";
$fault_info_r = "-m -W1.5p,255/0/0";

$sinfo = "-Sc4p -G255 -W0.5p,0/0/0";
$poleinfo1 = "-Sa20p -G0/255/255 -W0.75p,0/0/0";
$circleinfo = "-Sc25p -W1.0p,0/0/0,--";

#----------------------------------
# USER PARAMETERS

# base directory
$dir0 = "/home/carltape/compearth/surfacevel2strain";

# velocity, strain rate, etc
$vel_dir0 = "$dir0/matlab_output";

# file name, ticks, range for velocities
# idata
#   1:NASA REASON
#   2: CCMM
#   4: japan
#   10-13: strike-slip
#   20-23: rotational
#   ETC
$iregion = 3;      # region (1=west_us, 2=cal, 3=socal, 6=cascadia, 8=parkfield, 9=japan, 10=wedge)
$idata = 1;        # choose GPS dataset (see above)
$ndim = 2;         # ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 3;         # min grid order used
$qmax = 7;         # max grid order used
$basis = 1;        # 1 for wavelets; 2 for splines
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$iunrotate = 1;    # 1 to remove rotation; 0 to leave original field

$imask  = 1;       # plot mask
$ieuler = 1;       # plot euler poles
$ifault = 0;       # plot faults
$iplate = 0;       # plot plate boundaries

# REGION OF INTEREST
if ($iregion == 1 ) {
   $itag = "west_us";
   #$z_title = 1.9;
   $z_title = 2.2;
   $xtick1 = 2; $xtick2 = 0.25; $wid = 4.5;
   $origin = "-X1.5 -Y1.75";
   $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 2) {
   $itag = "cal";
   $z_title = 1.45;
   $xtick1 = 1; $xtick2 = 0.5; $wid = 6.0;
   $origin = "-X1.0 -Y1.75";
   $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 3) {
   $itag = "socal";
   #$z_title = 1.2;
   $z_title = 1.05;
   $xtick1 = 1; $xtick2 = 0.25; $wid = 6.0;
   #$origin = "-X1.0 -Y1.0";
   #$xlegend = 0.25; $ylegend = 1.3;
   $origin = "-X1.0 -Y2.5";
   $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 6) {
   $itag = "cascadia";
   $z_title = 3.55;
   $xtick1 = 2; $xtick2 = 0.5; $wid = 2.5;
   $origin = "-X2.0 -Y1.5";
   $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 8) {
   $itag = "parkfield";
   $z_title = 0.8;
   $xtick1 = 0.2; $xtick2 = 0.1; $wid = 6.0;
   $origin = "-X1.0 -Y2.5";
   $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 9) {
   $itag = "japan";
   $z_title = 1.2;
   $xtick1 = 2; $xtick2 = 0.5; $wid = 6.0;
   $origin = "-X1.0 -Y2.5";
   $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 10) {
   $itag = "wedge";
   $z_title = 1.85;
   $xtick1 = 20; $xtick2 = 10; $wid = 4.0;
   $origin = "-X2.0 -Y2.0";
   $xlegend = 0.25; $ylegend = -1.0;

}
$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";
$evec_ref = 1.0; $escale = 8;    # default values
$fac_min = 1.0; $fac_max = 1.0;  # default values
$vec_conf0 = 0.95;
$vec_error = 0.5;
$vec_conf = $vec_conf0;

# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
$dlab = sprintf("d%2.2i",$idata);
if ($idata==1) {
  $tag = "NASA REASON data set (cGPS)";
  $gps_pts = "${gps_dir}/US/reason_subset_points.dat";
  $gps_vec = "${gps_dir}/US/reason_subset_psvelo.dat";
  $igc = 0;
  #$cmin = 0; $cmax = 40; $ctick = 10;  # fixed NAM
  if($iunrotate==1) {    # rotation removed
     $cmin = 0; $cmax = 20; $ctick = 5; $vscale = 0.025; $vec_ref = 30;
  } else {                              # original REASoN
     $cmin = 5; $cmax = 50; $ctick = 10; $vscale = 0.01; $vec_ref = 50;
  }
  if($iregion==8) {$cmax = 45; $ctick = 10;}
  $vumin = -2; $vumax = 2; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;
  $ifault = 1; $iplate = 1;
  $fac_min = 0.5;    # 1.0 for uniform color scaling (try 0.5)
  $fac_max = 1.0;    # 1.0 default

} elsif ($idata==2) {
  $tag = "California Crustal Motion Map, v1.0  (Shen et al., SCEC-2006)";
  $gps_pts = "${gps_dir}/US/california/socal_vfield_4p0_points.dat";
  $gps_vec = "${gps_dir}/US/california/socal_vfield_4p0_psvelo.dat";
  $igc = 0; $vec_ref = 50; $vscale = 0.01;
  $cmin = 0; $cmax = 40; $ctick = 10;
  $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;

} elsif ($idata==4) {
  $tag = "Japan data set, 1996--2000";   # from Takeo Ito
  $gps_pts = "${gps_dir}/ASIA/japan/japan_takeo_ito_subset_points.dat";
  $gps_vec = "${gps_dir}/ASIA/japan/japan_takeo_ito_subset_psvelo.dat";
  $igc = 0;
  #$cmin = 0; $cmax = 40; $ctick = 10;  # 
  $cmin = 0; $cmax = 30; $ctick = 5; $vscale = 0.015; $vec_ref = 30;
  $vumin = -2; $vumax = 2; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;
  $ifault = 0; $iplate = 1;
  $fac_min = 0.5;    # 1.0 for uniform color scaling (try 0.5)
  $fac_max = 1.0;    # 1.0 default
}

# SYNTHETIC velocity field
if ($idata >= 10) {

  $gps_pts = "${gps_dir}/synthetic/syn_vfield_${dlab}_points.dat";
  $gps_vec = "${gps_dir}/synthetic/syn_vfield_${dlab}_psvelo.dat"; 

  # strike-slip field
  if ($idata >= 10 and $idata < 20) {
    $igc = 1; $ifault = 0; $iplate = 0;
    $vec_ref = 50; $vscale = 0.01;
    #$cmin = 0; $cmax = 40; $ctick = 10;
    $cmin = 0; $cmax = 20; $ctick = 5;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 2.0; $escale = 5;
    $fac_min = 0.50;   # 1.0 for uniform color scaling
    $fac_max = 1.0;    # 1.0 default
    if($iregion==10) {$cmin = 0; $cmax = 20; $ctick = 5; $vscale = 0.01; $vec_ref = 30; $evec_ref = 10; $escale = 1;}
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
    $igc = 0; $ifault = 0; $iplate = 0;
    $vec_ref = 10; $vscale = 0.06;
    $cmin = 0; $cmax = 8; $ctick = 2;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 1.0; $escale = 8;
    if($iregion==10) {$cmin = 0; $cmax = 80; $ctick = 20; $vscale = 0.006; $vec_ref = 40;}
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
    $igc = 1; $ifault = 0; $iplate = 0;
    $vec_ref = 1000; $vscale = 0.001;
    $cmin = 0; $cmax = 500; $ctick = 100;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 0.80;
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
    $igc = 1; $ifault = 0; $iplate = 0;
    $vec_ref = 1000; $vscale = 0.0005;
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

  # microplate rotation field, I
  if ($idata >= 60 and $idata < 70) {
    $igc = 1; $ifault = 0; $iplate = 0;
    $vec_ref = 10; $vscale = 0.06;
    $cmin = 0; $cmax = 5; $ctick = 1;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 0.2; $escale = 50;
    $fac_min = 1.0;   # 1.0 for uniform color scaling
    $fac_max = 1.0;    # 1.0 default
    if ($idata == 60) {
      $tag = "microplate rotational field (uniform points, no errors)"; $imask = 0;
    } elsif ($idata == 61) {
      $tag = "microplate rotational field (uniform points, with errors)"; $imask = 0;
    } elsif ($idata == 62) {
      $tag = "microplate rotational field (REASON points, no errors)";
    } elsif ($idata == 63) {
      $tag = "microplate rotational field (REASON points, with errors)";
    }
  }

  # microplate rotation field, II
  if ($idata >= 70 and $idata < 80) {
    $igc = 1; $ifault = 0; $iplate = 0;
    $vec_ref = 10; $vscale = 0.06;
    $cmin = 0; $cmax = 5; $ctick = 1;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 0.2; $escale = 50;
    $fac_min = 1.0;   # 1.0 for uniform color scaling
    $fac_max = 1.0;    # 1.0 default
    if ($idata == 70) {
      $tag = "microplate rotational field (uniform points, no errors)"; $imask = 0;
    } elsif ($idata == 71) {
      $tag = "microplate rotational field (uniform points, with errors)"; $imask = 0;
    } elsif ($idata == 72) {
      $tag = "microplate rotational field (REASON points, no errors)";
    } elsif ($idata == 73) {
      $tag = "microplate rotational field (REASON points, with errors)";
    }
  }

  # volcanic dilitation field
  if ($idata >= 80 and $idata < 90) {
    $igc = 0; $ifault = 0; $iplate = 0;
    $vec_ref = 40; $vscale = 0.006;
    $cmin = 0; $cmax = 40; $ctick = 10;
    $vumin = -60; $vumax = 60; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 4; $escale = 1;
    $fac_min = 1.0;   # 1.0 for uniform color scaling
    $fac_max = 1.0;    # 1.0 default
    if ($idata == 80) {
      $tag = "volcanic dilitation field (uniform points, no errors)"; $imask = 0;
    } elsif ($idata == 81) {
      $tag = "volcanic dilitation field (uniform points, with errors)"; $imask = 0;
    } elsif ($idata == 82) {
      $tag = "volcanic dilitation field (REASON points, no errors)";
    } elsif ($idata == 83) {
      $tag = "volcanic dilitation field (REASON points, with errors)";
    }
  }


}
if (not -f $gps_pts) { die("Check if $gps_pts exist or not\n") }
if (not -f $gps_vec) { die("Check if $gps_vec exist or not\n") }

#$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);
$vscale_text = sprintf("%.0f mm/yr", $vec_ref);

#$xtick1 = 1; $xtick2 = 0.25; $wid = 5.5; $cmin = 0; $cmax = 40; $ctick = 10;
#$xtick1 = 1; $xtick2 = 0.25; $wid = 3.25; $cmin = 0; $cmax = 40; $ctick = 10;
#$xtick1 = 1; $xtick2 = 0.25; $wid = 5.5; $cmin = 44; $cmax = 52; $ctick = 2;   # plate field (socal)

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

$phi = "\@~\146\@~";
$theta = "\@~\161\@~";

# plotting specificiations

# A : smallest feature plotted, in km^2; D : resolution
$coast_res     = "-A500 -Df";
$coast_infoK   = "$coast_res -W1.5p,0/0/0 -Na/1.0p";
$coast_infoW   = "$coast_res -W2.0p,255/255/255 -Na/2.0p,255/255/255,t";
$coast_infoR   = "$coast_res -W2.0p,255/0/0 -Na/2.0p,255/0/0,t";

#$coast_infoW   = "$coast_res -W1.0p,255/255/255";

$plate_infoG     = "-m -W2p,0/255/0";
$plate_infoR     = "-m -W2p,255/0/0";
$plate_infoK     = "-m -W2p,0/0/0";
$plate_infoW     = "-m -W2p,255/255/255";
$plate_infoGr    = "-m -W3p,200";

# NOTE: pen attribute (-W) does not work for psvelo
#$vec_scale       = "-Se0.002/0.90/0";
#$vec_scale       = "-Se0.005/0.90/0";
#$vec_infoK       = "-A0.5p/4p/2p -G0/0/0 -W0.5p,0/0/0 $vec_scale";      # add -N
#$vec_info        = "-A0.5p/4p/2p -G255/0/0 -W0.5p,255/0/0 $vec_scale";  # add -N
#$vec_info_scale  = "-N -A1p/8p/4p -G255/0/0 -W0.5p,255/0/0 $vec_scale";
#$vec_info_scaleK = "-N -A1p/8p/4p -G0/0/0 -W0.5p,0/0/0 $vec_scale";

$textinfo = "-G255 -S1p";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

#-----------------------------------------------
# BOOLEAN: WHICH FIGURES TO PLOT
$ixv    = 1;
$ipdf   = 0;

# single subplot figures
$ivel = 1;
$istrain = 0;
$irotate = 0;
$idilate = 0;
$imaxlam = 0;
$imap = 0;
$imap2 = 0;  # see utils/socal_centers.pl

# 2-column figures
$ivel3D = 0;         # estimated velocity field (Vmag, Vr)
#$ivelstr = 0;
#$imaskfig = 0;
#$isetup = 0;

# 3-column figures
$imulti_vel = 0;
$imulti_strain_inc = 0;    # only one at a time!
$imulti_strain_cum = 0;    # only one at a time!
$ivel3Dall = 0;
$ivel_norotate = 0;

if($imulti_vel == 1 && $ndim == 2) {die("\n imulti_vel = 1, so ndim must be 3\n")}

#==================================================

$glab = sprintf("%2.2i",$igc);
$qlab = sprintf("q%2.2i_q%2.2i",$qmin,$qmax);
$blab = sprintf("b%1.1i",$basis);
$nlab = sprintf("%1iD",$ndim);
$slab = sprintf("s%1i",$iLmat);
$ulab = sprintf("u%1i",$iunrotate);
#if ($idata==0) {       $slab = "Lfull";
#} elsif ($idata==1) {  $slab = "Lelastic";
#} elsif ($idata==2) {  $slab = "Lviscous";
#} else {
#  die("\nExit: invalid iLmat option\n\n");
#}

$name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}_${ulab}";

print "\n $name\n";

if($igc==1){
  $gc_boundary  = "/home/carltape/gmt/gps_data/synthetic/gps_gc_${dlab}.dat";
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
#if (not -f ${euler_poles_scale}) { die("Check if ${euler_poles_scale} exist or not\n") }

# get bounds of the region
open(IN,"$bounds"); @lines = <IN>;
($lonmin,$lonmax,$latmin,$latmax) = split(" ",$lines[0]);

# modify the bounds
if($iregion == 1) {
  $latmin = 32; $latmax = 50; $lonmin = -126; $lonmax = -114;
} elsif($iregion == 3) {
  $latmin = 32; $latmax = 38; $lonmin = -122; $lonmax = -114;
} elsif ($iregion == 8) {
  $latmin = 35.45; $latmax = 36.2; $lonmin = -121.2; $lonmax = -119.8;
}
$R1 = "-R$lonmin/$lonmax/$latmin/$latmax";
print "\n bounds are $R1\n\n";
#die("testing");

# get color limits for the plots: dilatation, strain, rotation, maxlam
open(IN,"$colors"); @lines = <IN>;
($cmindilat,$cmaxdilat,$cpwr2) = split(" ",$lines[0]);  # dilatation
($cminstrain,$cmaxstrain,$cpwr3) = split(" ",$lines[1]);  # strain
#$cmaxstrain = 8;
($cminrot,$cmaxrot,$cpwr4) = split(" ",$lines[2]);  # rotation
#$cmaxrot = 0.1;
($cmin5,$cmax5,$cpwr5) = split(" ",$lines[3]);  # maxlam
#$norm2 = "1e$cpwr2"; $norm3 = "1e$cpwr3"; $norm4 = "1e$cpwr4"; $norm5 = "1e$cpwr5";
$norm2 = 10**$cpwr2; $norm3 = 10**$cpwr3; $norm4 = 10**$cpwr4; $norm5 = 10**$cpwr5;
print "\n $cpwr2 $cpwr3 $cpwr4 $cpwr5 $norm2 $norm3 $norm4 $norm5\n\n";

$title_strain  = sprintf("Strain rate, 10\@+%i\@+ yr\@+-1\@+",$cpwr3);   # %2.2i for 07
$title_rot    = sprintf("Rotation rate, |w|, 10\@+%i\@+ rad/yr",$cpwr4);
$title_dilat  = sprintf("Dilatation rate, 10\@+%i\@+ yr\@+-1\@+",$cpwr2);

$title_up    = "Vr (mm/yr)";
$title_south = "V$theta (mm/yr)";
$title_east  = "V$phi (mm/yr)";
#$title_up    = "Estimated V\@-r\@- (mm/yr)";
#$title_south = "Estimated V\@-${theta}\@- (mm/yr)";
#$title_east  = "Estimated V\@-${phi}\@- (mm/yr)";

$Bscale_strain  = sprintf("-B%2.2f:\"Strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);
$Bscale_rot    = sprintf("-B%2.2f:\"Rotation rate, 10\@+%i\@+ rad/yr\": -Ef10p",$cmaxrot*0.25,$cpwr4);
$Bscale_dilat  = sprintf("-B%2.2f:\"Dilatation rate, 10\@+%i\@+ yr\@+-1\@+\": -E10p",$cmaxdilat*0.5,$cpwr2);
$unit_rot      = sprintf("10\@+%i\@+ yr\@+-1\@+",$cpwr4);

#$Bscale_strain2  = sprintf("-B%2.2f:\" \": -E5p",$cmaxstrain*0.5);
#$Bscale_rot2    = sprintf("-B%2.2f:\" \": -E5p",$cmaxrot*0.5);
#$Bscale_dilat2  = sprintf("-B%2.2f:\" \": -E5p",$cmaxdilat);

if(0==1){
  $cmindilat = $cmindilat*$norm2;
  $cmaxdilat = $cmaxdilat*$norm2;
  $cminstrain = $cminstrain*$norm3;
  $cmaxstrain = $cmaxstrain*$norm3;
  $cminrot = $cminrot*$norm4;
  $cmaxrot = $cmaxrot*$norm4;
  $cmin5 = $cmin5*$norm5;
  $cmax5 = $cmax5*$norm5;
}

# subplotting specifications
$hwid = $wid/2;
$xfac = 1.20;
$yfac = 1.65;
$dX = $xfac*$wid; $dY = 0; $shift1 = "-X$dX -Y$dY";
$dX = -$dX; $dY = -$yfac*$wid; $shift2 = "-X$dX -Y$dY";

$J1 = "-JM${wid}i";

# plot title
$J_title = "-JX${wid}";
$R_title = "-R0/1/0/1";
$x_title = 0.5;

$fsize_title = 14;

#-----------------------------------------
# make color files

$cshfile_color = "color.csh";
open(COLOR,">${cshfile_color}");

# surface velocities
$cran = $cmax - $cmin; $dc = $cran/100;
$Tvmag = "-T$cmin/$cmax/$dc";

# velocity components
$dc=($vumax-$vumin)/100; $Tup="-T$vumin/$vumax/$dc";     $vtick[0] = ($vumax-$vumin)/4;
$dc=($vsmax-$vsmin)/100; $Tsouth="-T$vsmin/$vsmax/$dc";  $vtick[1] = ($vsmax-$vsmin)/4;
$dc=($vemax-$vemin)/100; $Teast="-T$vemin/$vemax/$dc";   $vtick[2] = ($vemax-$vemin)/4;

print "\n @vtick\n";

# dilatation
#$cmaxdilat = $cmaxdilat/2; $cmindilat = $cmindilat/2;
$cran = $cmaxdilat - $cmindilat;
$dc = $cran/100;
$Tdilat = "-T$cmindilat/$cmaxdilat/$dc";

# strain
$cran = $cmaxstrain - $cminstrain;
$dc = $cran/100;
$Tstrain = "-T$cminstrain/$cmaxstrain/$dc";

# rotation
$cran = $cmaxrot - $cminrot;
$dc = $cran/100;
$Trot = "-T$cminrot/$cmaxrot/$dc";

# maxlam
$cran = $cmax5 - $cmin5;
$dc = $cran/100;
$T5 = "-T$cmin5/$cmax5/$dc";

# maxlam -- log10 scale
$T6 = "-T0.01/$cmax5/3";

$colorbar = "rainbow";
$colorbar = "seis -I";

# v-field components
$cptup = "color_up.cpt"; print COLOR "makecpt -C$colorbar $Tup -D > $cptup\n";
$cptsouth = "color_south.cpt"; print COLOR "makecpt -C$colorbar $Tsouth -D > $cptsouth\n";
$cpteast = "color_east.cpt"; print COLOR "makecpt -C$colorbar $Teast -D > $cpteast\n";

$cptvmag = "color1.cpt";
print COLOR "makecpt -C$colorbar $Tvmag -D > $cptvmag\n";
#print COLOR "sed 's/^B.*/B       255   255    255  /' temp1  >  temp2\n";
#print COLOR "sed 's/^F.*/F       255     0      0  /' temp2 > $cptvmag\n";

if(1==1) {
  $colorbar = $c_romania;
  $colorbar = "seis -I";
  #$cback = "  0     0    255";
  #$cfor  = "255     0      0";

} else {
  $colorbar = "no_green";
  #$cback = "170     0      0";
  #$cfor  = "  0     0    200";
}

# 170     0      0  to 255   255    255

# -I to invert the color scale

$colorbar = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";

$cptdilat = "color2.cpt";
print COLOR "makecpt -C$colorbar $Tdilat -D > $cptdilat\n";

$cptstrain = "color3.cpt";
print COLOR "makecpt -C$colorbar $Tstrain -D > $cptstrain\n";

$cptrot = "color4.cpt";
print COLOR "makecpt -C$colorbar $Trot -D > $cptrot\n";

$cptfile5 = "color5.cpt";
print COLOR "makecpt -C$colorbar $T5 -D > $cptfile5\n";

#$cptfile6 = "color6.cpt";
#print COLOR "makecpt -C$colorbar $T6 -Qo -D > $cptfile6\n";

close (COLOR);
system("csh -f ${cshfile_color}");

#-----------------------------------------

# color bar
#$Dlen = 0.4*$wid;
$Dlen = 2.5;
$Dx = $xlegend + 2.0 + $Dlen/2;
$Dy = $ylegend + 0.5;
#$Dscale = "-D$hwid/-0.6/$Dlen/0.15h";
$Dscale = "-D$Dx/$Dy/$Dlen/0.15h";

$B0 = "-Ba${xtick1}f${xtick2}d::";

$B = $B0.$Bopts[0];

# KEY: commands for interpolation and masking
$interp = "-I0.1"; $interp_surf = "-I0.04 -T0 -S0";
#if($iregion==1) {$interp_surf = "-I0.2 -T0 -S0";}
if($iregion==8) {$interp = "-I0.025"; $interp_surf = "-I0.01 -T0 -S0";}
if($iregion==9) {$interp = "-I0.2"; $interp_surf = "-I0.08 -T0 -S0";}
if($iregion==10) {$interp = "-I0.1"; $interp_surf = "-I0.4 -T0 -S0";}

$mask_info = "-I0.1 -G200 -S0.2";
if($iregion==8) {$mask_info = "-I0.01 -G200 -S0.02";}
if($iregion==9) {$mask_info = "-I0.2 -G200 -S0.4";}

$grdfile = "temp.grd";
$mask_file = "${vel_dir0}/${name}_masked_pts.dat";
if (not -f ${mask_file}) { die("Check if ${mask_file} exist or not\n") }

# open CSH file
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";

$pname = $name;
#$pname = "${name}_zoom";    # for comparing different functions to plot
#$pname = "${name}_weight";
#$pname = "${name}_weight_3D";

#-----------------------------------------

if ($istrain==1){

$fname = "${pname}_strain";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

$title = "Strain rate from $tag";
#$Bscale  = sprintf("-B%2.2e:\"Strain rate,  yr\@+-1\@+ \": -Ef10p",$cmaxstrain*0.25);
$Bscale  = sprintf("-B%2.2f:\"Strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);

# plot basemap
print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
#print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | pscontour $R1 $J1 -A- -C$cptstrain -I -O -K -V >> $psfile\n";

# make grd file, then plot
if (0==1) {
  #$interp = "-I0.1";
  print CSH "awk '{print \$1,\$2,\$5}' $strain_mag | xyz2grd -G$grdfile $interp $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -Sn -K -O -V >> $psfile\n";

} else {

  #print CSH "awk '{print \$1,\$2,\$5}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -Sn -K -O -V >> $psfile\n";

  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
    print CSH "psmask -C -K -O -V >> $psfile\n";
  }
}

# color bar
#print CSH "gmtset D_FORMAT \%.2e\n";
print CSH "psscale -C$cptstrain $Dscale $Bscale -K -O -V >> $psfile\n";
#print CSH "gmtset D_FORMAT \%lg\n";

print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

# plot observation points
#$dot_info = "-Sc4p -W0.5p";
#print CSH "awk '{ print \$1,\$2}' ${velocity_vec_dat} | psxy $J1 $R1 $dot_info -K -V -O >> $psfile\n";

# plot pole of rotational field
#if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";

#--------------------------------------
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============================

if ($irotate==1) {

$fname = "${pname}_rotation";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = "Rotation rate from $tag";
$Bscale  = sprintf("-B%2.2f:\"Rotation rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxrot*0.25,$cpwr4);

#$R1 = "-R-180/180/-70/70";

@wtags = ("w","wr","wh","ws","wt");

# make several different plots of rotation magnitudes
$wmin = 1; $wmax = 5;      # default
$wmin = 1; $wmax = $wmin;
for ($w = $wmin; $w <= $wmax; $w = $w+1) {

  $wtag = $wtags[$w-1];
  $fname = "${pname}_rotation_${wtag}";
  $title = "Rotation rate, |$wtag|, from $tag";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
  $Bscale  = sprintf("-B%2.2f:\"Rotation rate, |%s|, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxrot*0.25,$wtag,$cpwr4);

  print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
  if ($icolor==1) {

    # plot magnitude of w (3D) or its vertical or horizontal components
    if ($w==1) {
      #print CSH "awk '{print \$1,\$2,\$6/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
      print CSH "awk '{print \$1,\$2,sqrt(\$11*\$11+\$12*\$12+\$13*\$13)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
    } elsif ($w==2) {
      print CSH "awk '{print \$1,\$2,sqrt(\$11*\$11)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
    } elsif ($w==3) {
      print CSH "awk '{print \$1,\$2,sqrt(\$12*\$12+\$13*\$13)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
    } elsif ($w==4) {
      print CSH "awk '{print \$1,\$2,sqrt(\$14*\$14)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
    } elsif ($w==5) {
      print CSH "awk '{print \$1,\$2,sqrt(\$15*\$15)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
    }
    print CSH "grdimage $grdfile $R1 $J1 -C$cptrot -Sn -K -O -V >> $psfile\n";
  }
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
    print CSH "psmask -C -K -O -V >> $psfile\n";
  }
  print CSH "psscale -C$cptrot $Dscale $Bscale -K -O -V >> $psfile\n";
  print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

  # plot faults and/or plate boundary
  if ($ifault==1) {
    #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
    print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  }
  if ($iplate==1) {
    print CSH "psxy ${plate_file} $J1 $R1 $plate_infoW -K -O -V >> $psfile\n";
  }
  if ($igc==1) {
    print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";
  }

  # plot euler poles and anti-poles -- only for the first figure in the loop
  if ($ieuler == 1 && $w == 1) {
    $pinfo      = "-Sc4p -G255/255/255 -W0.5p,0/0/0";
    $pinfo_anti = "-Sc4p -G255/0/0 -W0.5p,0/0/0";
    $seis_info3 = "-Scp -W0.5p/255/255/255";

    # make binary colorpoint file
    $cptfile = "bcolor.cpt"; $cmin = -0.01; $cmax = 0.01; $dc = ($cmax-$cmin)/2;
    $T = "-T$cmin/$cmax/$dc -D -Z "; $colorbar = "-Cpolar";
    print CSH "makecpt $colorbar $T > $cptfile\n";

    if (0==1) {
      print CSH "echo checking euler poles file ${euler_poles}\n";
      print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_poles} > poles_test\n";
    }

    # KEY: plot euler poles
    print CSH "gmtset MEASURE_UNIT point\n";
    #print CSH "psxy ${euler_poles} $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
    #print CSH "psxy ${euler_anti_poles} $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
    print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
    print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_anti_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";

    #if($idata >= 20 && $idata < 30) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n";}

    #-----------------------------
    #$origin_arrow = "-Xa0.5 -Ya-1.5";

    # in inches
    $xdots = $xlegend;
    $ydots = $ylegend - 0.5;
    $origin_dots = "-Xa${xdots}i -Ya${ydots}i";

    # plot euler vector scale
    @vec_ref = ($evec_ref/2, $evec_ref, 1.5*$evec_ref);
    $sv1 = sprintf("%.2f",$vec_ref[0]);
    $sv2 = sprintf("%.2f",$vec_ref[1]);
    $sv3 = sprintf("%.2f",$vec_ref[2]);
    #open(IN,${euler_poles_scale}); $escale = <IN>; chomp($escale);
    #if (not -f ${euler_poles_scale}) {die("Check if ${euler_poles_scale} exist or not\n");}

    # MAKE SURE THAT THIS FORMULA MATCHES WHAT IS USED IN PLOTTING THE POLES ABOVE
    # the factor 72 is to convert inches to dots
    $dsize[0] = $escale*abs($vec_ref[0]);
    $dsize[1] = $escale*abs($vec_ref[1]);
    $dsize[2] = $escale*abs($vec_ref[2]);
    $dy = 0.25*72;
    for ($k = 0; $k < 5; $k = $k+1) {
      $ytxt[$k] = $k*$dy;
    }
    $x1 = 0*72;
    $x2 = $x1 + 1.00*72;
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
$x1 $ytxt[4] $fsize2 0 $fontno CM pos              neg
$x1 $ytxt[3] $fsize1 0 $fontno CM $sv3
$x1 $ytxt[2] $fsize1 0 $fontno CM $sv2
$x1 $ytxt[1] $fsize1 0 $fontno CM $sv1
$x1 $ytxt[0] $fsize2 0 $fontno CM $unit_rot
EOF\n";
    #-----------------------------

    print CSH "gmtset MEASURE_UNIT inch\n";

  }  # plot euler poles as dots

  if ( ($idata >= 60 and $idata < 80) || ($idata >= 20 && $idata < 30) ) {
    $epoles = "${vel_dir0}/${name}_epole_mean_points.dat";
    if (not -f $epoles) {
      die("Check if $epoles exist or not\n");
    }
    print CSH "psxy $epoles $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile\n";
  }

  if ($idata == 11 and $iregion == 10) {
    print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -170.6477 -41.5308\nEOF\n";
  }

  # plot spline centers
  #print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
  print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n"; # FINISH

  print CSH "echo done with $psfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if ($ixv==1) {
    print CSH "gv $psfile &\n";
  }

}				# for $w
}

#==============================

if ($idilate==1) {

$fname = "${pname}_dilatation";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = "Dilatation rate from $tag";
$Bscale  = sprintf("-B%2.2f:\"Dilatation rate, 10\@+%i\@+ yr\@+-1\@+\": -E10p",$cmaxdilat*0.5,$cpwr2);

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
print CSH "awk '{print \$1,\$2,\$4/$norm2}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
print CSH "grdimage $grdfile $R1 $J1 -C$cptdilat -Sn -K -O -V >> $psfile\n";
if ($imask==1) {
  print CSH "echo mask file is ${mask_file}\n";
  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
  print CSH "psmask -C -K -O -V >> $psfile\n";
}
print CSH "psscale -C$cptdilat $Dscale $Bscale -K -O -V >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

#if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============================

if ($imaxlam==1) {

$fname = "${pname}_maxlam";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = "Max strain rate from $tag";
$Bscale  = sprintf("-B%2.2f:\"Max strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmax5*0.25,$cpwr5);

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
print CSH "awk '{print \$1,\$2,\$7/$norm5}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
print CSH "grdimage $grdfile $R1 $J1 -C$cptfile5 -Sn -K -O -V >> $psfile\n";
if ($imask==1) {
  print CSH "echo mask file is ${mask_file}\n";
  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
  print CSH "psmask -C -K -O -V >> $psfile\n";
}
print CSH "psscale -C$cptfile5 $Dscale $Bscale -K -O -V >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============================

if ($ivel==1) {

$fname = "${pname}_vfield";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = $tag;
$Bscale  = "-B$ctick:\" Surface Vel Mag (mm/yr)\": -Ef10p";

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
print CSH "grdimage $grdfile $R1 $J1 -C$cptvmag -Sn -K -O -V >> $psfile\n";
if ($imask==1) {
  print CSH "echo mask file is ${mask_file}\n";
  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
  print CSH "psmask -C -K -O -V >> $psfile\n";
}
print CSH "psscale -C$cptvmag $Dscale $Bscale -K -O -V >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1 && $iregion==8) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
}
if ($iplate==1 && $iregion!=8) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoK -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

# KEY COMMAND: create data file
#$gpstemp = $gps_vec;
$gpstemp = gpstemp1;
print CSH "awk '{ print \$1,\$2,\$3,\$4,\$5,\$6}' ${velocity_vec_dat} > $gpstemp\n";

#$gpstemp = gpstemp1;
#if($igc==0) {
#  #print CSH "awk '{ print \$3,\$2,\$4,\$5,\$6,\$7,\$8,\$1}' $gps_file | tail +2 > $gpstemp\n";
#  $gpstemp = $gps_vec;
#} else {
#  print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec} > $gpstemp\n";
#}

# plot v-field
$vec_info = "-N -A0.5p/4p/2p -Se${vscale}/${vec_conf}/0 -W0.5p";
print CSH "psvelo $gpstemp $J1 $R1 $vec_info -K -V -O >> $psfile\n";

# plot spline centers
#print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

# plot v-field scale
print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0.25 $fsize2 0 $fontno LM ${vscale_text}\nEOF\n";
print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

# plot pole of rotational field
if($idata >= 20 && $idata < 30) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}
if ($idata == 11 and $iregion == 10) {
print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -170.6477 -41.5308\nEOF\n";
}

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
#-----
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============================

if ($imap==1) {

$fname = "${pname}_geometry";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = "$tag: Observations and basis functions";
#$title = "Geodetic observations used in SCEC Crustal Motion Map";

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
if ($imask==1) {
  print CSH "echo mask file is ${mask_file}\n";
  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
  print CSH "psmask -C -K -O -V >> $psfile\n";
}

# plot station locations
$ginfo = "-Si6p -G0/255/255 -W0.5p,0/0/0";
print CSH "psxy $gps_pts $J1 $R1 $ginfo -K -V -O >> $psfile\n";

# plot spline gridpoints
print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_r -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_r -K -V -O >> $psfile\n";
}
print CSH "psxy $J1 $R1 ${plate_file} $plate_infoR -K -O -V >> $psfile\n";

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
#--------
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============================

if ($imap2==1) {

  $fname = "${pname}_geometry2";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
  $title = "  ";

  # color scale for qmax field
  $cptqmax = "color_qmax.cpt";
  $Tq = "-T5.5/9.5/1"; $crad0 = 40;
  #$Tq = "-T5.5/8.5/1"; $crad0 = 35;
  print CSH "makecpt -Cseis $Tq -D -I > $cptqmax\n";

  print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START

  # plot colormap showing maxq at each plotting point
  $qmax_field = "${vel_dir0}/${name}_qmax_plotting.dat";
  if (not -f $qmax_field) {
    die("Check if qmax_field $qmax_field exist or not\n");
  }

  print CSH "awk '{print \$1,\$2,\$3}' $qmax_field | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptqmax -Sn -K -O -V >> $psfile\n";

  # plot station locations
  $ginfo = "-Si6p -G0/255/255 -W0.5p,0/0/0";
  print CSH "psxy $gps_pts $J1 $R1 $ginfo -K -V -O >> $psfile\n";

  # plot spline gridpoints
  $qvec = "${vel_dir0}/${name}_qvec.dat";
  if (not -f $qvec) {
    die("Check if qvec $qvec exist or not\n");
  }
  open(IN,"$qvec"); @qvecs = <IN>; print "\n @qvecs\n";
  $numq = @qvecs;
  $crad = $crad0;
  for ($i = 1; $i <= $numq; $i = $i+1) {
    $q = @qvecs[$i-1]; chomp($q);
    $stq = sprintf("%2.2i",$q);
    $qcen = "${vel_dir0}/${name}_gridpoints_q${stq}.dat";
    if (not -f $qcen) {
      die("Check if qcen $qcen exist or not\n");
    }
    $crad = $crad - 5;
    $ginfo = "-Sc${crad}p -W0.5p,0/0/0";
    print CSH "awk '{print \$1,\$2}' $qcen | psxy $qcen $J1 $R1 $ginfo -K -V -O >> $psfile\n";
  }

  print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
  #print CSH "psxy $J1 $R1 ${plate_file} $plate_infoR -K -O -V >> $psfile\n";

  # plot legend
  $crad = $crad0;
  for ($i = 1; $i <= $numq; $i = $i+1) {
    $q = @qvecs[$i-1]; chomp($q);
    $crad = $crad - 5;
    $ginfo = "-Sc${crad}p -W0.5p,0/0/0 -N";
    $xtxt = $i*0.2;
    $ytxt = 0;
    print CSH "psxy -JX1/1 -R0/1/0/1 $ginfo -K -O -V $origin_arrow >>$psfile<<EOF\n$xtxt $ytxt\nEOF\n";

  }

  print CSH "psscale -C$cptqmax $Dscale $Bscale -K -O -V >> $psfile\n";

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
  #--------
  print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n"; # FINISH

  print CSH "echo done with $psfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if ($ixv==1) {
    print CSH "gv $psfile &\n";
  }

}

#==============================
# from here on out, use portrait plotting with 2-column figures

@labs = ("A","B","C","D");
$fault_info_k    = "-m -W1.0p,0/0/0";
$fault_info_r    = "-m -W1.0p,255/0/0";
$dX1 = 0; $dY1 = 0;
$dX2 = 0; $dY2 = 0;

if ($iregion == 1 ) {
  $wid = 3; $vscale = 0.005;
  $origin = "-X1.0 -Y3";
  $xlab = 0.95*$wid; $ylab = 1.7*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $dX2 = 0; $dY2 = -3.6;
  $xlegend = 0.25; $ylegend = 0.4;

} elsif ($iregion == 2) {
  $wid = 3; $vscale = 0.005;
  $origin = "-X1.0 -Y7";
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 3) {
  $wid = 3; $vscale = 0.005;
  $origin = "-X1.0 -Y7";
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $xlegend = 0.25; $ylegend = -1.0;

} elsif ($iregion == 6) {
  $origin = "-X1.5 -Y1.5"; $vscale = 0.0002;
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $xlegend = 0.25; $ylegend = -0.75;

} elsif ($iregion == 8) {
  $origin = "-X1.5 -Y1.5"; $vscale = 0.0002;
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $xlegend = 0.25; $ylegend = -0.75;

} elsif ($iregion == 9) {
  $wid = 3; $vscale = 0.005;
  $origin = "-X1.0 -Y3";
  $xlab = 0.95*$wid; $ylab = 1.7*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $dX2 = 0; $dY2 = -3.6;
  $xlegend = 0.25; $ylegend = 0.4;
}
$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";

$hwid = $wid/2;
$fsize_title = 10;
$Dscale = "-D$hwid/-0.2/$Dlen/0.10h";
$J1 = "-JM${wid}i";
$J_title = "-JX${wid}";

$shift = "-X$dX1 -Y$dY1";
$shift2 = "-X-$dX2 -Y$dY2";

print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D MEASURE_UNIT inch TICK_LENGTH $tick LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p D_FORMAT \%lg\n";

#==============================

if ($ivel3D==1) {

  $fname = "${pname}_vel3D";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
  $title = "(a)  Horizontal velocity magnitude (mm/yr)";
  $Bscale  = "-B$ctick:\"  \": -Ef10p"; 
  $B = $B0.$Bopts[5];

  print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
  print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptvmag -Sn -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
    print CSH "psmask -C -K -O -V >> $psfile\n";
  }
  print CSH "psscale -C$cptvmag $Dscale $Bscale -K -O -V >> $psfile\n";
  print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  #print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";

  # plot data: horizontal velocity field
  if (0==1) {
    # KEY COMMAND: create data file
    $gpstemp = gpstemp1;
    if ($igc==0) {
      #print CSH "awk '{ print \$3,\$2,\$4,\$5,\$6,\$7,\$8,\$1}' $gps_file | tail +2 > $gpstemp\n";
      $gpstemp = $gps_vec;
    } else {
      print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec} > $gpstemp\n";
    }

    # plot v-field
    $vec_info = "-N -A0.5p/3p/1.5p -Se${vscale}/${vec_conf}/0 -W0.5p";
    print CSH "psvelo $gpstemp $J1 $R1 $vec_info -K -V -O >> $psfile\n";

    # plot v-field scale
    print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 -0.1 14 0 $fontno LM ${vscale_text}\nEOF\n";
    print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";
  }

  # highlight the small volcanic source
  #if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

  if ($ifault==1) {
    #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
    print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  }
  if ($iplate==1) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoW -K -O -V >> $psfile\n";}
  if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
  #print CSH "pstext -N $textinfo -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $xlab $ylab 14 0 $fontno CM $labs[0]\nEOF\n";

  #=======================

  $title = "(b)  Vertical velocity field (mm/yr)";
  $ctick = ($vumax - $vumin)/4;
  $Bscale  = "-B${ctick}:\" \": -E10p";
  $B = $B0.$Bopts[14];

  print CSH "psbasemap $J1 $R1 $B -K -V -O $shift >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$8}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptup -Sn -K -O -V >> $psfile\n";
  print CSH "psscale -C$cptup $Dscale $Bscale -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
    print CSH "psmask -C -K -O -V >> $psfile\n";
  }

  #print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
  #if($igc==0) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n"}
  #else {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n"}

  print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

  # highlight the small volcanic source
  #if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

  # locations of data points
  #$dinfo = "-Sc4p -G255 -W0.5p,0/0/0";
  #print CSH "psxy ${gps_pts} $J1 $R1 $dinfo -K -V -O >> $psfile\n";

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
  #print CSH "pstext -N $textinfo -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $xlab $ylab 14 0 $fontno CM $labs[1]\nEOF\n";

  #-----
  print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n"; # FINISH

  print CSH "echo done with $psfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if ($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============================

#if ($ivelstr==1) {

#$fname = "${pname}_velstr";
#$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
#$title = $tag;
#$Bscale  = "-B$ctick:\" Horizontal Velocity Magnitude (mm/yr)\": -Ef10p";
#$B = $B0.$Bopts[5];

#print CSH "gmtset BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.15c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";

#print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
#print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "grdimage $grdfile $R1 $J1 -C$cptvmag -Sn -K -O -V >> $psfile\n";
#if ($imask==1) {
#  print CSH "echo mask file is ${mask_file}\n";
#  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
#  print CSH "psmask -C -K -O -V >> $psfile\n";
#}
#print CSH "psscale -C$cptvmag $Dscale $Bscale -K -O -V >> $psfile\n";
#print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
##print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
##print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";

## KEY COMMAND: create data file
#$gpstemp = gpstemp1;
#if($igc==0) {
#  #print CSH "awk '{ print \$3,\$2,\$4,\$5,\$6,\$7,\$8,\$1}' $gps_file | tail +2 > $gpstemp\n";
#  $gpstemp = $gps_vec;
#} else {
#  print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec} > $gpstemp\n";
#}

## plot v-field
#$vscale = 0.005;
#$vec_info = "-N -A0.5p/3p/1.5p -Se${vscale}/${vec_conf}/0 -W0.5p";
#print CSH "psvelo $gpstemp $J1 $R1 $vec_info -K -V -O >> $psfile\n";

## plot v-field scale
#print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 -0.1 $fsize2 0 $fontno LM ${vscale_text}\nEOF\n";
#print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

#if($igc==0) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n"}
#else {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n"}

##print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
#print CSH "pstext -N $textinfo -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $xlab $ylab 14 0 $fontno CM $labs[0]\nEOF\n";

##=======================

#$title = "Strain rate from $tag";
##$Bscale  = sprintf("-B%2.2e:\"Strain rate,  yr\@+-1\@+ \": -Ef10p",$cmaxstrain*0.25);
#$Bscale  = sprintf("-B%2.2f:\"Strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);
#$B = $B0.$Bopts[14];

#print CSH "psbasemap $J1 $R1 $B -K -V -O $shift >> $psfile\n";

#print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";

## make grd file, then plot
##print CSH "awk '{print \$1,\$2,\$5}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -Sn -K -O -V >> $psfile\n";

#if ($imask==1) {
#  print CSH "echo mask file is ${mask_file}\n";
#  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
#  print CSH "psmask -C -K -O -V >> $psfile\n";
#}

## color bar
##print CSH "gmtset D_FORMAT \%.2e\n";
#print CSH "psscale -C$cptstrain $Dscale $Bscale -K -O -V >> $psfile\n";
##print CSH "gmtset D_FORMAT \%lg\n";

#print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
##print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";

##if($igc==0) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n"}
##else {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n"}

#print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

#if($igc==0) {
#   print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n";
#   print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
#   print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
#}
#else {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n"}

##print CSH "pstext -N -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
#print CSH "pstext -N $textinfo -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $xlab $ylab 14 0 $fontno CM $labs[1]\nEOF\n";

##-----
#print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

#print CSH "echo done with $psfile\n";
#print CSH "convert $psfile $jpgfile\n";
#if($ixv==1) {print CSH "gv $psfile &\n"}

#}

##==============================

#if ($imaskfig==1) {

#$fname = "${pname}_mask";
#$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
#$title = $tag;
#$title = " ";
#$Bscale  = "-B$ctick:\" Surface Velocity Magnitude (mm/yr)\": -Ef10p";

#print CSH "gmtset BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.15c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";

#for ($k = 0; $k <= 1; $k = $k+1) {

#if($k==0) {
#   $B = $B0.$Bopts[5];
#   print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";   # START
#} else {
#   $B = $B0.$Bopts[6];
#   print CSH "psbasemap $J1 $R1 $B -K -V -O $shift >> $psfile\n";
   
#}
#print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "grdimage $grdfile $R1 $J1 -C$cptvmag -Sn -K -O -V >> $psfile\n";

#if($k==1) {
#print CSH "echo mask file is ${mask_file}\n";
#print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
#print CSH "psmask -C -K -O -V >> $psfile\n";
#}

#print CSH "psscale -C$cptvmag $Dscale $Bscale -K -O -V >> $psfile\n";
#print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

#if(0==1) {

## KEY COMMAND: create data file
#$gpstemp = gpstemp1;
#if($igc==0) {$gpstemp = $gps_vec;
#} else {
#  print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec} > $gpstemp\n";
#}

## plot v-field
#$vec_info = "-N -A0.5p/6p/3p -Se${vscale}/${vec_conf}/0 -W0.5p";
#print CSH "psvelo $gpstemp $J1 $R1 $vec_info -K -V -O >> $psfile\n";

## plot v-field scale
#print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 -0.1 9 0 $fontno LM $vec_ref mm/yr\nEOF\n";
#print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

#}

#if($igc==0) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n"}
#else {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n"}

##print CSH "pstext -N -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
#print CSH "pstext -N $textinfo -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $xlab $ylab 14 0 $fontno CM $labs[$k]\nEOF\n";

#}  # for k

##================================

#$title = "Strain rate from $tag";
##$Bscale  = sprintf("-B%2.2e:\"Strain rate,  yr\@+-1\@+ \": -Ef10p",$cmaxstrain*0.25);
#$Bscale  = sprintf("-B%2.2f:\"Strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);


#for ($k = 0; $k <= 1; $k = $k+1) {

#if($k==0) {
#   $B = $B0.$Bopts[1];
#   print CSH "psbasemap $J1 $R1 $B -K -V -O $shift2 >> $psfile\n";
#} else {
#   $B = $B0.$Bopts[3];
#   print CSH "psbasemap $J1 $R1 $B -K -V -O $shift >> $psfile\n";
   
#}
#print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";

## make grd file, then plot
##print CSH "awk '{print \$1,\$2,\$5}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -Sn -K -O -V >> $psfile\n";

#if($k==1) {
#print CSH "echo mask file is ${mask_file}\n";
#print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
#print CSH "psmask -C -K -O -V >> $psfile\n";
#}

## color bar
##print CSH "gmtset D_FORMAT \%.2e\n";
#print CSH "psscale -C$cptstrain $Dscale $Bscale -K -O -V >> $psfile\n";
##print CSH "gmtset D_FORMAT \%lg\n";

#print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
##print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";

#if($igc==0) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n"}
#else {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n"}

##print CSH "pstext -N -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
#print CSH "pstext -N $textinfo -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $xlab $ylab 14 0 $fontno CM $labs[$k+2]\nEOF\n";

#}  # for k

##-----
#print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

#print CSH "echo done with $psfile\n";
#print CSH "convert $psfile $jpgfile\n";
#if($ixv==1) {print CSH "gv $psfile &\n"}

#}

##==============================

#if ($isetup==1) {

#$fname = "setup";
#$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
#$title = "Topography, bathymetry, faults, and GPS stations";

#$origin = "-X0.75 -Y6.5";

## change for 2-box poster figures
#print CSH "gmtset BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.2c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";

#print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START

## plot station locations
#$ginfo = "-Sc6p -G0/255/255 -W0.5p,0/0/0";
#$gpstemp = $gps_pts;
#print CSH "psxy $gps_pts $J1 $R1 $ginfo -K -V -O >> $psfile\n";

## plot spline gridpoints
##print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

#print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
##print CSH "psxy $kcf_file $J1 $R1 $fault_info_r -K -V -O >> $psfile\n";
#print CSH "psxy $fault_file $J1 $R1 $fault_info_r -K -V -O >> $psfile\n";
##print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";

#print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title 10 0 $fontno CM $title\nEOF\n";

##--------

#$shift = "-X0 -Y-4.5";
#$B = "-Ba${xtick1}f${xtick2}d::".$Bopts[11];
#$Bscale  = "-B$ctick:\" Surface Velocity Magnitude (mm/yr)\": -Ef10p";
#$title = "SCEC Crustal Motion Map, v3.0";

#print CSH "psbasemap $J1 $R1 $B -K -V -O $shift >> $psfile\n";
#print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "grdimage $grdfile $R1 $J1 -C$cptvmag -Sn -K -O -V >> $psfile\n";
#if ($imask==1) {
#  print CSH "echo mask file is ${mask_file}\n";
#  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
#  print CSH "psmask -C -K -O -V >> $psfile\n";
#}
#print CSH "psscale -C$cptvmag $Dscale $Bscale -K -O -V >> $psfile\n";
#print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
##print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
##print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";

## plot v-field
#$vec_info = "-N -A0.5p/6p/3p -Se${vscale}/${vec_conf}/0 -W0.5p";
#print CSH "psvelo $gps_vec $J1 $R1 $vec_info -K -V -O >> $psfile\n";

## plot v-field scale
#print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0.25 $fsize1 0 $fontno LM $vec_ref mm/yr\nEOF\n";
#print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

#print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n";

#print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title 10 0 $fontno CM $title\nEOF\n";

##--------
#print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

#print CSH "echo done with $psfile\n";
#print CSH "convert $psfile $jpgfile\n";
#if($ixv==1) {print CSH "gv $psfile &\n"}

#}

#==============================
# from here on out, use portrait plotting with 3-column figures

@labs = ("A","B","C","D");
$fault_info_k  = "-m -W0.5p,0/0/0";
$fault_info_r  = "-m -W0.5p,255/0/0";
$coast_infoK   = "$coast_res -W1.0p,0/0/0 -Na/1.0p";
$coast_infoW   = "$coast_res -W1.0p,255/255/255 -Na/1.0p,255/255/255,t";
$plate_infoK   = "-m -W1p,0/0/0";
$plate_infoW   = "-m -W1p,255/255/255";
$plate_infoR   = "-m -W1p,255/0/0";

if ($iregion == 1 ) {
  $wid = 2.0; $vscale = 0.005;
  $origin = "-X0.75 -Y3";
  $xfac = 1.20; $yfac = 1.0;
  $dX1 = $xfac*$wid; $dY1 = 0; 
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $xlegend = 0.25; $ylegend = -1.0;
  $xtick1 = 2; $xtick2 = 2;
  $x_title = 0.05; $z_title = 2.2;
  $x_stitle = 0.05; $z_stitle = -0.08; $fsize_stitle = 8;   # subtitle
  $Dscale = "-D1/-0.1/1.5/0.08h";
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH 3p LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 7 HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";

} elsif ($iregion == 2) {
  $wid = 3; $vscale = 0.005;
  $origin = "-X1.0 -Y7";
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $xlegend = 0.25; $ylegend = -1.0;
  $Dscale = "-D1.4/-0.05/1.0/0.08h";

} elsif ($iregion == 3) {
  $wid = 2.0; $vscale = 0.005;
  $origin = "-X0.75 -Y8.75";
  $xfac = 1.20; $yfac = 1.0;
  $dX1 = $xfac*$wid; $dY1 = 0; 
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $xlegend = 0.25; $ylegend = -1.0;
  $xtick1 = 1; $xtick2 = 1;
  $x_title = 0.5; $z_title = 0.87;
  $x_stitle = 0.0; $z_stitle = -0.08; $fsize_stitle = 10;   # subtitle
  $Dscale = "-D1.4/-0.08/1.0/0.08h";
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH 4p LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 8 HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";
  

} elsif ($iregion == 6) {
  $wid = 1.75; $hwid = $wid/2;
  $origin = "-X1.00 -Y3.5";
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $xlegend = 0.25; $ylegend = -0.75;
  $x_title = 0.; $z_title = 3.57;
  $Dscale = "-D$hwid/-0.15/1.5/0.12h";
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH 5p LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 8 HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";

} elsif ($iregion == 8) {
  $wid = 1.75; $hwid = $wid/2;
  $origin = "-X1.00 -Y3.5";
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $dX1 = $wid + 0.5; $dY1 = 0;
  $xlegend = 0.25; $ylegend = -0.75;
  $x_title = 0.; $z_title = 3.57;
  $Dscale = "-D$hwid/-0.15/1.5/0.12h";
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH 5p LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 8 HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";
}
$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";

if ($iregion == 3) {
  $latmin = 32; $latmax = 37; $lonmin = -122; $lonmax = -114;
  $R1 = "-R$lonmin/$lonmax/$latmin/$latmax";
}

$fsize_title = 10;
$J1 = "-JM${wid}i";
$J_title = "-JX${wid}";

$dX2 = -2*$dX1; $dY2 = -$yfac*$wid;
$shift1 = "-X$dX1 -Y$dY1";
$shift2 = "-X$dX2 -Y$dY2";

$B0 = "-Ba${xtick1}f${xtick2}d::";
$B = $B0.$Bopts[5];

@flab = ("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p");

#==============================

if (  $imulti_vel == 1 || ($imulti_strain_inc == 1 || $imulti_strain_cum == 1) ) {

  # components and magnitudes of v-field
  $vu_file = "${vel_dir0}/${name}_vfield_up_plot.dat";
  $vs_file = "${vel_dir0}/${name}_vfield_south_plot.dat";
  $ve_file = "${vel_dir0}/${name}_vfield_east_plot.dat";
  if ($ndim == 3) {if (not -f $vu_file) {die("Check if $vu_file exist or not\n")};}
  if (not -f $vs_file) {die("Check if $vs_file exist or not\n")}
  if (not -f $ve_file) {die("Check if $ve_file exist or not\n")}

  # multiscale strain fields (cumulative and incremental)
  $multi_strain_inc = "${vel_dir0}/${name}_multi_strain_inc.dat";
  $multi_strain_cum = "${vel_dir0}/${name}_multi_strain_cum.dat";
  if (not -f $multi_strain_inc) {die("Check if $multi_strain_inc exist or not\n")}
  if (not -f $multi_strain_cum) {die("Check if $multi_strain_cum exist or not\n")}

  # additional multiscale files
  $cticks = "${vel_dir0}/${name}_tick_colors.dat";
  $colors = "${vel_dir0}/${name}_vel_colors.dat";
  $iqs    = "${vel_dir0}/${name}_iqvec.dat";
  $iqs2   = "${vel_dir0}/${name}_iqvec2.dat";
  if (not -f $colors) {die("Check if $colors exist or not\n")}
  if (not -f $iqs) {die("Check if $iqs exist or not\n")}
  if (not -f $iqs2) {die("Check if $iqs2 exist or not\n")}

  # get color limits for the plots
  open(IN,"$colors"); @clines = <IN>; print "\n @clines\n";

  # get color limits for the plots
  open(IN,"$cticks"); @cticks = <IN>; print "\n @cticks\n";

  # get q indexes for the plots
  open(IN,"$iqs"); @qlines1 = <IN>; print "\n @qlines1";
  open(IN,"$iqs2"); @qlines2 = <IN>; print "\n @qlines2";
  $nump = @qlines1;
  for ($i = 0; $i < $nump; $i = $i+1) {
    ($q1,$q2) = split(" ",$qlines1[$i]);
    $qtags1a[$i] = "q${q1}_q${q2}";
    $qtags1b[$i] = "q = ${q1}-${q2}";
  }
  for ($i = 0; $i < $nump-1; $i = $i+1) {
    ($q1,$q2) = split(" ",$qlines2[$i]);
    $qtags2a[$i] = "q${q1}_q${q2}";
    $qtags2b[$i] = "q = ${q1}-${q2}";
  }
  print "\n -- nump = $nump --\n @qtags1a\n @qtags1b\n @qtags2a\n @qtags2b\n";;
  #die("testing");

  #-----------------------------------------
  # make color files

  #$colorbar = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";
  $colorbar = "seis -I";
  $cshfile_color = "color.csh";
  open(COLOR,">${cshfile_color}");

  # components of multiscale velocity field
  for ($i = 0; $i < $nump; $i = $i+1) {
    ($cmin1,$cmax1,$cmin2,$cmax2,$cmin3,$cmax3) = split(" ",$clines[$i]);
    $dc = ($cmax1 - $cmin1)/100; $T1 = "-T$cmin1/$cmax1/$dc";
    $dc = ($cmax2 - $cmin2)/100; $T2 = "-T$cmin2/$cmax2/$dc";
    $dc = ($cmax3 - $cmin3)/100; $T3 = "-T$cmin3/$cmax3/$dc";
    $cptup = "color_up${i}.cpt";
    $cptsouth = "color_south${i}.cpt";
    $cpteast = "color_east${i}.cpt";
    print COLOR "makecpt -D -C$colorbar $T1 > $cptup\n";
    print COLOR "makecpt -D -C$colorbar $T2 > $cptsouth\n";
    print COLOR "makecpt -D -C$colorbar $T3 > $cpteast\n";
  }
  close (COLOR);
  system("csh -f ${cshfile_color}");
}

#==============

if ($imulti_vel == 1) {

  $fname = "${pname}_multi_vel";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
  $title = $tag;

  $title1 = $title_up;
  $title2 = $title_south;
  $title3 = $title_east;

  #@titles = ("(a)  $title_up","(b)  $title_south","(c)  $title_east");
  @titles = ("$title_up","$title_south","$title_east");

  for ($i = 0; $i < $nump; $i = $i+1) {

    $p = $i + 1;
    $col = $p + 2;
    if($i==0) {$msk = $nump-1} else {$msk = $i}
    $qtag = $qtags1a[$i];

    # flag if the mask file is empty
    if ($imask == 1) {
      $mask_file = "${vel_dir0}/${name}_masked_pts_${msk}.dat";
      if (not -f $mask_file) {die("check if $mask_file exists");}
      #open(IN,"$mask_file"); @test = <IN>;
      ($nline,$nnum,$nchar,$junk) = split(" ",`wc $mask_file`);
    } else {
      $nline = 0;
    }

    ($ctick1,$ctick2,$ctick3) = split(" ",$cticks[$i]);

    $cptfile = "color_up${i}.cpt";
    $Bscale  = "-B$ctick1:\" \": -E5p";
    #$Bscale  = "-B$ctick1:\" Surface Speed (mm/yr)\": -E5p";

    $B = $B0.$Bopts[15];
    if ($i==0) {
      print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
    } else {
      print CSH "psbasemap $J1 $R1 $B -K -O -V $shift2 >> $psfile\n";
    }
    #print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}}' $vu_file | surface -G$grdfile ${interp_surf} $R1\n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -Sn -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile\n";
      print CSH "psmask -C -K -O -V >> $psfile\n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

    if ($ifault == 1 && $iplate == 1) {
      if ($i==0 || $i==$nump-1) { 
	print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
	#print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
      } else {
	print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
      }
    }
    if ($igc==1) {
      print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";
    }
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile\n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {
      print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[0]\nEOF\n";
    }

    # subtitle
    $stitle = "($flab[3*$i])  $qtags1b[$i]";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle\nEOF\n";

    #----------------------------------------
    $B = $B0.$Bopts[9];
    $cptfile = "color_south${i}.cpt";
    $Bscale  = "-B$ctick2:\" \": -E5p";

    print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}}' $vs_file | surface -G$grdfile ${interp_surf} $R1\n";
      #print CSH "surface ${vs_file} -G$grdfile ${interp_surf} $R1\n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -Sn -K -O -V >> $psfile\n";
    }

    #print "\n $cptfile\n $vs_file\n column ${col}\n"; die("testing");

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile\n";
      print CSH "psmask -C -K -O -V >> $psfile\n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
      if ($i==0 || $i==$nump-1) {
	print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
	#print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
      } else {
	print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
      }
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile\n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {
      print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[1]\nEOF\n";
    }

    # subtitle
    $stitle = "($flab[3*$i+1])  $qtags1b[$i]";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle\nEOF\n";

    #----------------------------------------
    $B = $B0.$Bopts[15];
    $cptfile = "color_east${i}.cpt";
    $Bscale  = "-B$ctick3:\" \": -E5p";

    print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}}' $ve_file | surface -G$grdfile ${interp_surf} $R1\n";
      #print CSH "surface ${ve_file} -G$grdfile ${interp_surf} $R1\n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -Sn -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile\n";
      print CSH "psmask -C -K -O -V >> $psfile\n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
      if ($i==0 || $i==$nump-1) {
	print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
	#print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
      } else {
	print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
      }
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile\n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[2]\nEOF\n";}

    # subtitle
    $stitle = "($flab[3*$i+2])  $qtags1b[$i]";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle\nEOF\n";

  }    # end for loop

print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title\nEOF\n"; # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============

if ($imulti_strain_inc == 1 || $imulti_strain_cum == 1) {

  if ($imulti_strain_inc == 1) {$multi_strain = $multi_strain_inc; $xtag = "inc";}
  if ($imulti_strain_cum == 1) {$multi_strain = $multi_strain_cum; $xtag = "cum";}

  $fname = "${pname}_multi_strain_${xtag}";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

  #@titles = ("(a)  $title_strain","(b)  $title_rot","(c)  $title_dilat");
  @titles = ("$title_strain","$title_rot","$title_dilat");
  $B = $B0.$Bopts[15];

  # KEY: adjust limits for multiscale strain plots
  #$fac_min = 0.50;   # 1.0 for uniform color scaling
  #$fac_max = 1.0;    # 1.0 default
  $fac_min = 0.5; $fac_max = 1.0;
  $fac_inc = ($fac_max - $fac_min)/($nump-2);
  print "\n Color scaling for multiscale strain: $fac_min, $fac_inc, $fac_max\n";

  for ($i = 0; $i < $nump-1; $i = $i+1) {

    $p = $i + 1;
    $col = $p + 2;
    if ($imulti_strain_inc == 1) {$qtag = $qtags1b[$i+1]} else {$qtag = $qtags2b[$i];}

    # flag if the mask file is empty
    if ($imask == 1) {
      $mask_file = "${vel_dir0}/${name}_masked_pts_${p}.dat";
      if (not -f $mask_file) {die("check if $mask_file exists");}
      #open(IN,"$mask_file"); @test = <IN>;
      ($nline,$nnum,$nchar,$junk) = split(" ",`wc $mask_file`);
    } else {
      $nline = 0;
    }

    # adjust the plotting range
    # SHOULD WE "SIMPLY" OUTPUT THE RANGES FROM SURFACEVEL2STRAIN.M?
    $colorbar = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";
    $fac = $fac_min + ($p-1)*$fac_inc;
    $cmaxstrain2 = $cmaxstrain*$fac; $cminstrain2 = 0;
    $cmaxrot2   = $cmaxrot*$fac;   $cminrot2 = 0;
    $cmaxdilat2 = $cmaxdilat*$fac; $cmindilat2 = -$cmaxdilat2;

    $c2 = 1.01*$cmaxstrain2;
    $dc = ($cmaxstrain2 - $cminstrain2)/100; $Tstrain2 = "-T$cminstrain2/$c2/$dc";
    $c2 = 1.01*$cmaxrot2;
    $dc = ($cmaxrot2 - $cminrot2)/100; $Trot2 = "-T$cminrot2/$c2/$dc";
    $c1 = 1.01*$cmindilat2; $c2 = 1.01*$cmaxdilat2;
    $dc = ($c2 - $c1)/100; $Tdilat2 = "-T$c1/$c2/$dc";
    $cptstrain2 = "cpt_strain.cpt"; print CSH "makecpt -C$colorbar $Tstrain2 -D > $cptstrain2\n";
    $cptrot2 = "cpt_rot.cpt"; print CSH "makecpt -C$colorbar $Trot2 -D > $cptrot2\n";
    $cptdilat2 = "cpt_dilat.cpt"; print CSH "makecpt -C$colorbar $Tdilat2 -D > $cptdilat2\n";

    $Bscale_strain2  = sprintf("-B%2.2f:\" \": -Ef5p",$cmaxstrain2*0.5);
    $Bscale_rot2    = sprintf("-B%2.2f:\" \": -Ef5p",$cmaxrot2*0.5);
    $Bscale_dilat2  = sprintf("-B%2.2f:\" \": -E5p",$cmaxdilat2);

    print "\n $norm $fac\n";

    if ($i==0) {
      print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
    } else {
      print CSH "psbasemap $J1 $R1 $B -K -O -V $shift2 >> $psfile\n";
    }

    #----------------------------------------
    # dilatation
    $B = $B0.$Bopts[15]; $cptfile = $cptdilat2; $Bscale = $Bscale_dilat2; $col = 2 + 3*$p; $norm = $norm2*$fac;

    #print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}/$norm}' $multi_strain | surface -G$grdfile ${interp_surf} $R1\n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -Sn -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile\n";
      print CSH "psmask -C -K -O -V >> $psfile\n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
      if ($i==$nump-2) {
	print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
	#print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
      } else {
	print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
      }
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile\n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[2]\nEOF\n";}

    # subtitle
    $stitle = "($flab[3*$i])  $qtag";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle\nEOF\n";

    #----------------------------------------
    # strain
    $B = $B0.$Bopts[9]; $cptfile = $cptstrain2; $Bscale  = $Bscale_strain2; $col = 2 + 3*$p - 2; $norm = $norm3*$fac;

    print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}/$norm}' $multi_strain | surface -G$grdfile ${interp_surf} $R1\n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -Sn -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile\n";
      print CSH "psmask -C -K -O -V >> $psfile\n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
      if ($i==$nump-2) {
	print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
	#print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
      } else {
	print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
      }
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile\n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {
      print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[0]\nEOF\n";
    }

    # subtitle
    $stitle = "($flab[3*$i+1])  $qtag";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle\nEOF\n";

    #----------------------------------------
    # rotation
    $B = $B0.$Bopts[15]; $cptfile = $cptrot2; $Bscale = $Bscale_rot2; $col = 2 + 3*$p - 1; $norm = $norm4*$fac;

    print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}/$norm}' $multi_strain | surface -G$grdfile ${interp_surf} $R1\n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -Sn -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile\n";
      print CSH "psmask -C -K -O -V >> $psfile\n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
      if ($i==$nump-2) {
	print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
	#print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
      } else {
	print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
      }
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile\n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {
      print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[1]\nEOF\n";
    }

    # subtitle
    $stitle = "($flab[3*$i+2])  $qtag";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle\nEOF\n";

  }    # end for loop

print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title\nEOF\n"; # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============================

if ($ivel3Dall==1) {

  $fname = "${pname}_vel3Dall";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

  #@titles = ("(a)  Estimated V\@-r\@-","(b)  Estimated V\@-s\@-","(c)  Estimated V\@-e\@-");
  @titles = ("(a)  $title_up","(b)  $title_south","(c)  $title_east");
  @Bins = (6,2,5);
  @cpts = ($cptup,$cptsouth,$cpteast);

  for ($i = 1; $i <= 3; $i = $i+1) {

    $col = $i+7;
    $iB = $Bins[$i-1];
    $B = $B0.$Bopts[$iB];
    $ctick = $vtick[$i-1];
    $Bscale  = "-B$ctick:\"  \": -E5p";

    if ($i==1) {
      print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";
    } else {
      print CSH "psbasemap $J1 $R1 $B -K -V -O $shift1 >> $psfile\n";
    }
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$$col}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cpts[$i-1] -Sn -K -O -V >> $psfile\n";
    }
    print CSH "psscale -C$cpts[$i-1] $Dscale $Bscale -K -O -V >> $psfile\n";
    if ($imask==1) {
      print CSH "echo mask file is ${mask_file}\n";
      print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
      print CSH "psmask -C -K -O -V >> $psfile\n";
    }

    #print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";
    #if($igc==0) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n"}
    #else {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n"}

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

    if ($iplate==1) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n";}
    if ($ifault==1) {
      #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
      print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

    # locations of data points
    #$dinfo = "-Sc4p -G255 -W0.5p,0/0/0";
    #print CSH "psxy ${gps_pts} $J1 $R1 $dinfo -K -V -O >> $psfile\n";

    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno LM $titles[$i-1]\nEOF\n";
    #print CSH "pstext -N $textinfo -JX1/1 -R0/1/0/1 -K -O -V >>$psfile<<EOF\n $xlab $ylab 14 0 $fontno CM $labs[1]\nEOF\n";

  }

  #-----
  print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n"; # FINISH

  print CSH "echo done with $psfile\n";
  print CSH "convert $psfile $jpgfile\n";
  if ($ixv==1) {
    print CSH "gv $psfile &\n";
  }

}

#==============================

if ($ivel_norotate==1 && $iunrotate==1) {

# get the velocity fields
$velocity_vec_before_rotate = "${vel_dir0}/${name}_vec_horz_before_unrotate.dat";
$velocity_vec_rotate = "${vel_dir0}/${name}_vec_horz_rotate.dat";
$euler_vec_unrotate = "${vel_dir0}/${name}_unrotate_euler_vector_psxy.dat";
$euler_vec_unrotate_anti = "${vel_dir0}/${name}_unrotate_euler_anti_vector_psxy.dat";

if (not -f $velocity_vec_before_rotate) {die("Check if $velocity_vec_before_rotate exist or not\n")}
if (not -f $velocity_vec_rotate) {die("Check if $velocity_vec_rotate exist or not\n")}
if (not -f $euler_vec_unrotate) {die("Check if $euler_vec_unrotate exist or not\n")}
if (not -f $euler_vec_unrotate_anti) {die("Check if $euler_vec_unrotate_anti exist or not\n")}

#----------------

$imerge = 0;  # cool plot showing three fields superimposed
$B0 = "-Ba2f1:\" \":/a2f1:\" \":";             # socal
#$B0 = "-Ba0.4f0.2:\" \":/a0.4f0.2:\" \":";    # Parkfield

#@Bins = (6,2,5);
@Bins = (5,14,6,8);
$origin = "-X0.75 -Y7.75";
@titles = ("(a)  Original field","(b)  Uniform rotational field","(c)  Residual field  =  (a) - (b)","(d)  All three fields");
$x_title = 0.05; $z_title = 0.95;
$origin_arrow = "-Xa0.1 -Ya-0.15";

$fname = "${pname}_vfield_norotate";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = $tag;
$Bscale  = "-B$ctick:\" Surface Vel Mag (mm/yr)\": -Ef10p";

# THESE MUST BE CHANGED MANUALLY FOR EACH RUN
# vscale is for the original field and the pure rotation field
# vscaleR is for the residual field
$vscale = 0.008; $vscaleR = $vscale; $vec_refR = $vec_ref;          # REASON, socal
#$vscale = 0.02; $vscaleR = 1.0; $vec_refR = 0.5;                  # uniform rotation field + errors
#$vscale = 0.005; $vscaleR = 0.017;  $vec_refR = $vec_ref;        # REASON, parkfield
$vec_info0 = "-N -L -A0.25p/3p/1.5p";
if ($imerge==1) {
  $vec_info1 = "$vec_info0 -Se${vscale}/${vec_conf}/0 -W0.25p -G0";
  $vec_info2 = "$vec_info0 -Se${vscale}/${vec_conf}/0 -W0.25p,255/0/0,-- -G255/0/0";
  $vec_infoR = "$vec_info0 -Se${vscaleR}/${vec_conf}/0 -W0.25p,0/255/255,-- -G0/255/255";
} else {
  $vec_info1 = "$vec_info0 -Se${vscale}/${vec_conf}/0 -W0.25p -G0";
  $vec_info2 = "$vec_info0 -Se${vscale}/${vec_conf}/0 -W0.25p -G0";
  $vec_infoR = "$vec_info0 -Se${vscaleR}/${vec_conf}/0 -W0.25p -G0";
}
$vscale_text = sprintf("%.0f mm/yr", $vec_ref);
$vscale_textR = sprintf("%.0f mm/yr", $vec_refR);

#----------------

$i = 1;
$iB = $Bins[$i-1];
$B = $B0.$Bopts[$iB];

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
if ($iplate==1 && $iregion!=8) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoR -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoR -K -O -V >> $psfile\n";}

# plot velocity field
$gpstemp = gpstemp1;
#print CSH "awk '{print \$1,\$2,\$3,\$4,\$5,\$6}' ${velocity_vec_before_rotate} > $gpstemp\n";
print CSH "awk '{print \$1,\$2,\$3,\$4}' ${velocity_vec_before_rotate} > $gpstemp\n";
print CSH "psvelo $gpstemp $J1 $R1 $vec_info1 -K -V -O >> $psfile\n";

# plot v-field scale
print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0 10 0 $fontno LM ${vscale_text}\nEOF\n";
print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info1 $origin_arrow -K -O -V >>$psfile<<EOF\n 0.75 0 $vec_ref 0.00\nEOF\n";
#print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info1 $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno LM $titles[$i-1]\nEOF\n";

#--------------------------

$i = 2;
$iB = $Bins[$i-1];
$B = $B0.$Bopts[$iB];

print CSH "psbasemap $J1 $R1 $B -K -V -O $shift1 >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
if ($iplate==1 && $iregion!=8) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoR -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoR -K -O -V >> $psfile\n";}

$gpstemp = gpstemp1;
print CSH "awk '{print \$1,\$2,\$3,\$4}' ${velocity_vec_rotate} > $gpstemp\n";
print CSH "psvelo $gpstemp $J1 $R1 $vec_info2 -K -V -O >> $psfile\n";
print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0 10 0 $fontno LM ${vscale_text}\nEOF\n";
print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info2 $origin_arrow -K -O -V >>$psfile<<EOF\n 0.75 0 $vec_ref 0.00\nEOF\n";

if($idata >= 20 && $idata < 30) {
  print CSH "awk '{print \$1,\$2}' ${euler_vec_unrotate} | psxy $J1 $R1 $poleinfo1 -K -V -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' ${euler_vec_unrotate_anti} | psxy $J1 $R1 $poleinfo1 -K -V -O >> $psfile\n";
  #print CSH "awk '{print \$1,\$2,\$3}' ${euler_vec_unrotate} | psxy $R1 $J1 -O -K -V >> $psfile\n";
}

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno LM $titles[$i-1]\nEOF\n";

#--------------------------

$i = 3;
$iB = $Bins[$i-1];
$B = $B0.$Bopts[$iB];

print CSH "psbasemap $J1 $R1 $B -K -V -O $shift1 >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
if ($iplate==1 && $iregion!=8) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoR -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoR -K -O -V >> $psfile\n";}

$gpstemp = gpstemp1;
print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec_dat} > $gpstemp\n";
print CSH "psvelo $gpstemp $J1 $R1 $vec_infoR -K -V -O >> $psfile\n";
print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0 10 0 $fontno LM ${vscale_textR}\nEOF\n";
print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_infoR $origin_arrow -K -O -V >>$psfile<<EOF\n 0.75 0 $vec_refR 0.00\nEOF\n";

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno LM $titles[$i-1]\nEOF\n";

#--------------------------

if ($imerge==1) {
  $i = 4;
  $iB = $Bins[$i-1];
  $B = "-Ba1f0.25d::$Bopts[$iB]";

  $wid = 6.5;
  $J1 = "-JM${wid}i";
  $origin_arrow = "-Xa0.2 -Ya0.2";

  # THESE MUST BE CHANGED MANUALLY FOR EACH RUN
  # vscale is for the original field and the pure rotation field
  # vscaleR is for the residual field
  $vscale = 2*$vscale; $vscaleR = $vscale; $vec_refR = $vec_ref; # REASON, socal
  #$vscale = 0.02; $vscaleR = 1.0; $vec_refR = 0.5;                  # uniform rotation field + errors
  #$vscale = 0.005; $vscaleR = 0.017;  $vec_refR = $vec_ref;        # REASON, parkfield
  $vec_info1 = "$vec_info0 -Se${vscale}/${vec_conf}/0 -W0.25p -G0";
  $vec_info2 = "$vec_info0 -Se${vscale}/${vec_conf}/0 -W0.25p,255/0/0,-- -G255/0/0";
  $vec_infoR = "$vec_info0 -Se${vscaleR}/${vec_conf}/0 -W0.25p,0/255/255,-- -G0/255/255";
  $vscale_text = sprintf("%.0f mm/yr", $vec_ref);
  $vscale_textR = sprintf("%.0f mm/yr", $vec_refR);

  print CSH "psbasemap $J1 $R1 $B -K -V -O -X-4.5 -Y-6 >> $psfile\n";
  print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
  if ($iplate==1 && $iregion!=8) {print CSH "psxy ${plate_file} $J1 $R1 $plate_infoR -K -O -V >> $psfile\n";}
  if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoR -K -O -V >> $psfile\n";}

  print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec_before_rotate} | psvelo $J1 $R1 ${vec_info1} -K -V -O >> $psfile\n";
  print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec_rotate} | psvelo $J1 $R1 ${vec_info2} -K -V -O >> $psfile\n";
  print CSH "awk '{ print \$1,\$2,\$3,\$4}' ${velocity_vec_dat} | psvelo $J1 $R1 ${vec_infoR} -K -V -O >> $psfile\n";

  print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0 10 0 $fontno LM ${vscale_textR}\nEOF\n";
  print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_infoR -W0.25p,0/0/0 $origin_arrow -K -O -V >>$psfile<<EOF\n 0.75 0 $vec_refR 0.00\nEOF\n";

  #print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno LM $titles[$i-1]\nEOF\n";
}

#--------------------------

#-----
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}


#==================================================

close (CSH);
system("csh -f $cshfile");
#print "convert $psfile $jpgfile\n";
#system("convert $psfile $jpgfile");
#if($ipdf==1) {system("ps2pdf $psfile")};
#if($ixv==1) {system("gv $psfile &")};

#==================================================
