#!/usr/bin/perl -w

#==========================================================
#
#  socal_centers.pl
#  Carl Tape
#  11-June-2009
#  
#  This script inputs a surface velocity field and strain scalar fields
#  computed in test_platemodel2strain.m and outputs a figure.
#  
#==========================================================

$cshfile = "socal_centers.csh";

$icolor = 1;    # ccc

# plates and faults
#$plate_dir = "/home/carltape/gmt/plates";
#$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
$plate_file  = "./bird_boundaries";
if (not -f ${plate_file}) { die("Check if ${plate_file} exist or not\n") }
#$kcf_file     = "/home/carltape/gmt/faults/kcf.xy";
$fault_file   = "/home/carltape/gmt/faults/jennings_more.xy";
if (not -f $fault_file) { die("Check if $fault_file exist or not\n") }

$fault_info_k = "-M -W1.5p,0/0/0";
$fault_info_r = "-M -W1.5p,255/0/0";

$sinfo = "-Sc4p -G255 -W0.5p,0/0/0";
$poleinfo1 = "-Sa20p -G0/255/255 -W0.75p,0/0/0";
$circleinfo = "-Sc25p -W1.0p,0/0/0,--";

# velocity, strain rate, etc
$vel_dir0   = "/home/carltape/SURFACEVEL2STRAIN/matlab_output";
#$vel_dir0   = "${plate_dir}/surface_velocities";

#----------------------------------

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
$qmax = 12;         # max grid order used
$basis = 1;        # 1 for wavelets; 2 for splines
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$iunrotate = 1;    # 1 to remove rotation; 0 to leave original field

$imask  = 1;       # plot mask
$ieuler = 0;       # plot euler poles
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
   $xtick1 = 1; $xtick2 = 0.25; $wid = 5.5;
   $origin = "-X1.0 -Y2.5";
   $xlegend = $wid + 0.5; $ylegend = 0;

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

# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
$dlab = sprintf("d%2.2i",$idata);
if ($idata==1) {
  $tag = "NASA REASON data set (cGPS)";
  $gps_pts = "${gps_dir}/US/reason_subset_points.dat";
  $gps_vec = "${gps_dir}/US/reason_subset_psvelo.dat";
  $igc = 0; $vec_conf = $vec_conf0;
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
  $igc = 0; $vec_ref = 50; $vec_conf = $vec_conf0; $vscale = 0.01;
  $cmin = 0; $cmax = 40; $ctick = 10;
  $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;

} elsif ($idata==4) {
  $tag = "Japan data set, 1996--2000";   # from Takeo Ito
  $gps_pts = "${gps_dir}/ASIA/japan/japan_takeo_ito_subset_points.dat";
  $gps_vec = "${gps_dir}/ASIA/japan/japan_takeo_ito_subset_psvelo.dat";
  $igc = 0; $vec_conf = $vec_conf0;
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
    $vec_ref = 50; $vec_conf = $vec_conf0; $vscale = 0.01;
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
    $vec_ref = 10; $vec_conf = 0.90; $vscale = 0.06;
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
    $vec_ref = 1000; $vec_conf = $vec_conf0; $vscale = 0.001;
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
    $vec_ref = 1000; $vec_conf = $vec_conf0; $vscale = 0.0005;
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
    $vec_ref = 10; $vec_conf = $vec_conf0; $vscale = 0.06;
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
    $vec_ref = 10; $vec_conf = $vec_conf0; $vscale = 0.06;
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
    $vec_ref = 40; $vec_conf = $vec_conf0; $vscale = 0.006;
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

$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);

#$xtick1 = 1; $xtick2 = 0.25; $wid = 5.5; $cmin = 0; $cmax = 40; $ctick = 10;
#$xtick1 = 1; $xtick2 = 0.25; $wid = 3.25; $cmin = 0; $cmax = 40; $ctick = 10;
#$xtick1 = 1; $xtick2 = 0.25; $wid = 5.5; $cmin = 44; $cmax = 52; $ctick = 2;   # plate field (socal)

#--------------------------------------------

# plotting specifications
$fsize0 = "20";
$fsize1 = "14";
$fsize2 = "12";
$fsize3 = "10";
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

$coast_infoW   = "$coast_res -W1.0p,255/255/255";

$plate_infoG     = "-M -W2p,0/255/0";
$plate_infoR     = "-M -W2p,255/0/0";
$plate_infoK     = "-M -W2p,0/0/0";
$plate_infoW     = "-M -W2p,255/255/255";
$plate_infoGr    = "-M -W3p,200";

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

$ivel = 0;
$istrain = 0;
$irotate = 0;
$idilate = 0;
$imaxlam = 0;
$imap = 0;
$imap2 = 1;

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
$title_rot    = sprintf("Rotation rate, 10\@+%i\@+ rad/yr",$cpwr4);
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
$Dlen = 1.5;
$Dx = $xlegend;
$Dy = $ylegend + 3.5;
#$Dscale = "-D$hwid/-0.6/$Dlen/0.15h";
$Dscale = "-D$Dx/$Dy/$Dlen/0.15";

$B = "-Ba${xtick1}f${xtick2}d::WeSn";

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
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize3 ANOT_FONT_SIZE $fsize3  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen CHAR_ENCODING Standard\n";

$pname = $name;
#$pname = "${name}_zoom";    # for comparing different functions to plot
#$pname = "${name}_weight";
#$pname = "${name}_weight_3D";

#-----------------------------------------

$fname = "${pname}_geometry2a";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = "  ";

#$Bscale  = sprintf("-B%2.2f:\"Strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);

# color scale for qmax field
$cptqmax = "color_qmax.cpt";
$Tq = "-T5.5/9.5/1"; $crad0 = 40;
#$Tq = "-T5.5/8.5/1"; $crad0 = 35;
if($qmax == 11) {$Tq = "-T5.5/11.5/1"; $crad0 = 50;}
if($qmax == 12) {$Tq = "-T5.5/12.5/1"; $crad0 = 50;}
print CSH "makecpt -Cseis $Tq -D -I > $cptqmax\n";

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START

# plot colormap showing maxq at each plotting point
$qmax_field = "${vel_dir0}/${name}_qmax_plotting.dat";
if (not -f $qmax_field) { die("Check if qmax_field $qmax_field exist or not\n")}

#print CSH "awk '{print \$1,\$2,\$3}' $qmax_field | surface -G$grdfile ${interp_surf} $R1\n";
#print CSH "grdimage $grdfile $R1 $J1 -C$cptqmax -T -K -O -V >> $psfile\n";
#print CSH "grdcontour $grdfile $R1 $J1 -C$cptqmax -A- -W1p -K -O -V >> $psfile\n";

# plot spline gridpoints
  $qvec = "${vel_dir0}/${name}_qvec.dat";
  if (not -f $qvec) { die("Check if qvec $qvec exist or not\n")}
  open(IN,"$qvec"); @qvecs = <IN>; print "\n @qvecs\n";
  $numq = @qvecs;
  $crad = $crad0;
  for ($i = 1; $i <= $numq; $i = $i+1) {
     $q = @qvecs[$i-1]; chomp($q);
     $stq = sprintf("%2.2i",$q);
     $qcen = "${vel_dir0}/${name}_gridpoints_q${stq}.dat";
     if (not -f $qcen) { die("Check if qcen $qcen exist or not\n")}
     $crad = $crad - 5;
     $ginfo = "-Sc${crad}p -W0.5p,0/0/0";
     print CSH "awk '{print \$1,\$2}' $qcen | psxy $qcen $J1 $R1 $ginfo -K -V -O >> $psfile\n";
  }

print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
#print CSH "psxy $J1 $R1 ${plate_file} $plate_infoR -K -O -V >> $psfile\n";

# plot station locations
$ginfo = "-Si6p -G0/255/255 -W0.5p,0/0/0";
print CSH "psxy $gps_pts $J1 $R1 $ginfo -K -V -O >> $psfile\n";

# plot dots
  $crad = $crad0;
  for ($i = 1; $i <= $numq; $i = $i+1) {
     $q = @qvecs[$i-1]; chomp($q);
     $crad = $crad - 5;
     $ginfo = "-Sc${crad}p -W0.5p,0/0/0 -N";
     $xtxt = 0;
     $ytxt = $i*0.3;
     print CSH "psxy -JX1/1 -R0/1/0/1 $ginfo -K -O -V $origin_arrow >>$psfile<<EOF\n$xtxt $ytxt \nEOF\n";
     $xtxt = 0.5; $stq = sprintf("q = %i",$q);
     print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n$xtxt $ytxt 11 0 $fontno CM $stq\nEOF\n";
  }

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";

print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH
#-----------------------------------------

$fname = "${pname}_geometry2b";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = "  ";

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START

print CSH "awk '{print \$1,\$2,\$3}' $qmax_field | surface -G$grdfile ${interp_surf} $R1\n";
print CSH "grdimage $grdfile $R1 $J1 -C$cptqmax -T -K -O -V >> $psfile\n";
print CSH "grdcontour $grdfile $R1 $J1 -C$cptqmax -A- -W1p -K -O -V >> $psfile\n";

print CSH "pscoast $J1 $R1 $B $coast_infoK -K -O -V >> $psfile\n";
#print CSH "psxy $J1 $R1 ${plate_file} $plate_infoR -K -O -V >> $psfile\n";

# plot station locations
$ginfo = "-Si6p -G0/255/255 -W0.5p,0/0/0";
print CSH "psxy $gps_pts $J1 $R1 $ginfo -K -V -O >> $psfile\n";

# plot color bar
$ginfo = "-Ss32p -W1p,0/0/0 -N";
  for ($i = 4; $i <= $numq; $i = $i+1) {
     $q = @qvecs[$i-1]; chomp($q);
     $xtxt = 0;
     #$ytxt = 2.0 + $i*0.3;
     $ytxt = 0.0 + $i*0.3;
     print CSH "psxy -JX1/1 -R0/1/0/1 -C$cptqmax $ginfo -K -O -V $origin_arrow >>$psfile<<EOF\n$xtxt $ytxt $q\nEOF\n";
     $xtxt = 0.5;
     #$stq = sprintf("q = %i",$q);
     $stq = sprintf("q \@~\243\@~ %i",$q);
     print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n$xtxt $ytxt 11 0 $fontno CM $stq\nEOF\n";
  }

#print CSH "psscale -C$cptqmax $Dscale $Bscale -K -O -V >> $psfile\n";

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";

#--------
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

#==================================================

close (CSH);
system("csh -f $cshfile");

#==================================================
