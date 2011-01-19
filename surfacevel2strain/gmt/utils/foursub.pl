#!/usr/bin/perl -w

#==========================================================
#
#  foursub.pl
#  Carl Tape
#  23-Aug-2007
#  
#  This script inputs a surface velocity field and strain scalar fields
#  computed in test_platemodel2strain.m and outputs a figure.
#  
#==========================================================

$cshfile = "foursub.csh";

$icolor = 1;    # ccc

# plates and faults
$plate_dir = "/home/carltape/gmt/plates";
$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
if (not -f ${plate_file}) { die("Check if ${plate_file} exist or not\n") }
#$kcf_file     = "/home/carltape/gmt/faults/kcf.xy";
$fault_file   = "/home/carltape/gmt/faults/jennings_more.xy";
if (not -f $fault_file) { die("Check if $fault_file exist or not\n") }

$fault_info_k = "-M -W0.75p,0/0/0";
$fault_info_r = "-M -W0.75p,255/0/0";
$fault_info_w = "-M -W0.75p,255/255/255";

$sinfo = "-Sc4p -G255 -W0.5p,0/0/0";
$poleinfo1 = "-Sa20p -G0/255/255 -W0.75p,0/0/0";
$poleinfo2 = "-Sc10p -G255 -W0.75p,0/0/0";
$circleinfo = "-Sc25p -W1.0p,0/0/0,--";

# velocity, strain rate, etc
#$vel_dir0   = "${plate_dir}/surface_velocities";
$vel_dir0   = "/home/carltape/SURFACEVEL2STRAIN/matlab_output";

#----------------------------------

# key parameters
# NOTE: THESE ARE OVER-RIDDEN IN THE LOOP BELOW
$iregion = 3;      # region (1=west_us, 2=cal, 3=socal, 6=cascadia)
$idata = 1;        # choose GPS dataset (1 = NASA REASON; 2 = CCMM; 1X = strike-slip; 2X = rotational)
$ndim = 2;         # ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 3;         # min grid order used
$qmax = 8;         # max grid order used
$basis = 1;        # 1 for wavelets; 2 for splines
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$iunrotate = 1;    # 1 to remove rotation; 0 to leave original field

$imask  = 1;       # plot mask
$ieuler = 0;       # plot euler poles
$ifault = 1;       # plot faults
$iplate = 0;       # plot plate boundaries

$itag = "socal";

$wid = 3.0;
$xtick1 = 1; $xtick2 = 0.25;
$ytick1 = 1; $ytick2 = $xtick2;
$origin = "-X1.0 -Y6.5";
$xlegend = 2; $ylegend = -0.6;

$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";
$evec_ref = 1.0; $escale = 8;    # default values
$fac_min = 1.0; $fac_max = 1.0;  # default values
$vec_conf0 = 0.95;
$vec_error = 0.5;                # ellipse error (mm/yr)

#@tlabs = ("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)");

# location of SAFOD hole
$plat = 35.974;
$plon = -120.552;

#=========================================
# KEY LOOP
#=========================================
#$kmin = 1; $kmax = 6;
$kmin = 6; $kmax = $kmin;

for ($k = $kmin; $k <= $kmax; $k = $k+1) {

if ($k == 1) {          # uniform rotation
   #$idata = 21; $ndim = 2; $imask  = 0; $qmin = 3; $qmax = 7; $iunrotate = 0;
   $idata = 23; $ndim = 2; $imask  = 1; $qmin = 3; $qmax = 7; $iunrotate = 0;
   $ntag = "rotation";

} elsif ($k == 2 ) {    # microplate 1
   $idata = 61; $ndim = 2; $imask  = 0; $qmin = 3; $qmax = 7; $iunrotate = 0; $ieuler = 1;
   $ntag = "micro1";

} elsif ($k == 3) {     # microplate 2
   $idata = 71; $ndim = 2; $imask  = 0; $qmin = 3; $qmax = 7; $iunrotate = 0; $ieuler = 1;
   $ntag = "micro2";

} elsif ($k == 4) {     # dilatation
   $idata = 83; $ndim = 3; $imask  = 1; $qmin = 3; $qmax = 7; $iunrotate = 0;
   $ntag = "dilatation";

} elsif ($k == 5) {     # infinite strike-slip fault
   $idata = 13; $ndim = 2; $imask  = 1; $qmin = 3; $qmax = 7; $iunrotate = 0; $ieuler = 1;
   #$idata = 12; $ndim = 2; $imask  = 1; $qmin = 3; $qmax = 7; $iunrotate = 0;
   $ntag = "strikeslip";

#-----------------------

} elsif ($k == 6) {
   $idata = 1; $ndim = 2; $imask  = 1; $qmin = 3; $qmax = 7; $iunrotate = 1; $ieuler = 0;  # w >= 1e-7
   $ntag = "reason_2D";

} elsif ($k == 7) {
   $idata = 1; $ndim = 3; $imask  = 1; $qmin = 3; $qmax = 7; $iunrotate = 1; $ieuler = 0;  # w >= 1e-7
   $ntag = "reason_3D";

} elsif ($k == 8) {
   $idata = 1; $ndim = 3; $imask  = 1; $qmin = 3; $qmax = 8; $iunrotate = 1;

} elsif ($k == 9) {
   $idata = 1; $ndim = 3; $imask  = 1; $qmin = 3; $qmax = 7; $iunrotate = 0;

} elsif ($k == 10) {
   $idata = 1; $ndim = 2; $imask  = 1; $qmin = 3; $qmax = 8; $iunrotate = 1;

}

#$ntag = sprintf("%2.2i",$k);
$ntag = sprintf("%2.2i_m%1i",$k,$imask);

# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
$dlab = sprintf("d%2.2i",$idata);
$igc = 0; $ifault = 1; $iplate = 0;
if ($idata==1) {
  $tag = "NASA REASON data set (cGPS)";
  $gps_pts = "${gps_dir}/US/reason_subset_points.dat";
  $gps_vec = "${gps_dir}/US/reason_subset_psvelo.dat";
  #$gps_pts = "${gps_dir}/US/reason_fixed_NAM_subset_points.dat";
  #$gps_vec = "${gps_dir}/US/reason_fixed_NAM_subset_psvelo.dat";
  $igc = 0; $vec_conf = $vec_conf0;
  $evec_ref = 2.0; $escale = 4;
  if($iunrotate==1) {    # rotation removed
     $cmin = 0; $cmax = 25; $ctick = 10; $vscale = 0.010; $vec_ref = 30;
  } else {                              # original REASoN
     $cmin = 5; $cmax = 50; $ctick = 10; $vscale = 0.01; $vec_ref = 50;
  }
  $vumin = -2; $vumax = 2; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;

} elsif ($idata==2) {
  $tag = "California Crustal Motion Map, v1.0  (Shen et al., SCEC-2006)";
  $gps_pts = "${gps_dir}/US/california/socal_vfield_4p0_points.dat";
  $gps_vec = "${gps_dir}/US/california/socal_vfield_4p0_psvelo.dat";
  $igc = 0; $vec_ref = 50; $vec_conf = $vec_conf0; $vscale = 0.01;
  $cmin = 0; $cmax = 40; $ctick = 10;
  $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;
}

# SYNTHETIC velocity field
if ($idata >= 10) {

  $gps_pts = "${gps_dir}/synthetic/syn_vfield_${dlab}_points.dat";
  $gps_vec = "${gps_dir}/synthetic/syn_vfield_${dlab}_psvelo.dat"; 

  # strike-slip field
  if ($idata >= 10 and $idata < 20) {
    $igc = 1; $ifault = 0; $iplate = 0;
    $vec_ref = 20; $vec_conf = $vec_conf0; $vscale = 0.01;
    $cmin = 0; $cmax = 22; $ctick = 5;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 2.0; $escale = 2;
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

  # rotational field
  if ($idata >= 20 and $idata < 30) {
    $igc = 0; $ifault = 0; $iplate = 0;
    $vec_ref = 10; $vec_conf = $vec_conf0; $vscale = 0.06;
    $cmin = 0; $cmax = 7; $ctick = 2;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 1.0; $escale = 8;
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
    $vec_ref = 10; $vec_conf = $vec_conf0; $vscale = 0.04;
    $cmin = 0; $cmax = 6; $ctick = 2;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 0.2; $escale = 30;
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
    $vec_ref = 10; $vec_conf = $vec_conf0; $vscale = 0.04;
    $cmin = 0; $cmax = 6; $ctick = 2;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 0.2; $escale = 30;
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
    $cmin = 0; $cmax = 45; $ctick = 10;
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

# text label for reference arrow in legend
#$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);
$vscale_text = sprintf("%.0f mm/yr", $vec_ref);

#--------------------------------------------

# plotting specifications
$fsize0 = "16";
$fsize1 = "14";
$fsize2 = "11";
$fsize3 = "9";
$fontno = "1";    # 1 or 4
$tick   = "0.15c";
$fpen   = "1p";
$tpen   = "1p";

$phi = "\@~\146\@~";
$theta = "\@~\161\@~";

# plotting specificiations

# A : smallest feature plotted, in km^2; D : resolution
$coast_res     = "-A500 -Df";
$coast_infoK   = "$coast_res -W1.5p,0/0/0 -Na/1.0p";
$coast_infoW   = "$coast_res -W1.5p,255/255/255 -Na/1.5p,255/255/255,t";
$coast_infoR   = "$coast_res -W1.5p,255/0/0 -Na/1.5p,255/0/0,t";

$plate_infoG     = "-M -W1.5p,0/255/0";
$plate_infoR     = "-M -W1.5p,255/0/0";
$plate_infoK     = "-M -W1.5p,0/0/0";
$plate_infoW     = "-M -W1.5p,255/255/255";
$plate_infoGr    = "-M -W1.5p,200";

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
$ulab = sprintf("u%1i",$iunrotate);

$name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}_${ulab}";

print "\n $name \n";

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
if (not -f $velocity_vec_dat) { die("Check if ${velocity_vec_dat} exist or not\n") }
if (not -f $bounds)         { die("Check if $bounds exist or not\n") }
if (not -f $colors)         { die("Check if $colors exist or not\n") }
if (not -f ${spline_centers}) { die("Check if ${spline_centers} exist or not\n") }
if (not -f ${euler_poles}) { die("Check if ${euler_poles} exist or not\n") }
if (not -f ${euler_anti_poles}) { die("Check if ${euler_anti_poles} exist or not\n") }
#if (not -f ${euler_poles_scale}) { die("Check if ${euler_poles_scale} exist or not\n") }

# get bounds of the region
open(IN,"$bounds"); @lines = <IN>; close(IN);
($lonmin,$lonmax,$latmin,$latmax) = split(" ",$lines[0]);

# modify the bounds
if($iregion == 1) {
  $latmin = 32; $latmax = 50; $lonmin = -126; $lonmax = -114;
} elsif($iregion == 3) {
  $latmin = 32; $latmax = 37; $lonmin = -121.9; $lonmax = -114;
}

$R1 = "-R$lonmin/$lonmax/$latmin/$latmax";
print "\n bounds are $R1 \n\n";
#die("testing");

# get color limits for the plots: dilatation, strain, rotation, maxlam
open(IN,"$colors"); @lines = <IN>; close(IN);
($cmindilat,$cmaxdilat,$cpwr2) = split(" ",$lines[0]);  # dilatation
($cminstrain,$cmaxstrain,$cpwr3) = split(" ",$lines[1]);  # strain
($cminrot,$cmaxrot,$cpwr4) = split(" ",$lines[2]);  # rotation
($cmin5,$cmax5,$cpwr5) = split(" ",$lines[3]);  # maxlam
#$norm2 = "1e$cpwr2"; $norm3 = "1e$cpwr3"; $norm4 = "1e$cpwr4"; $norm5 = "1e$cpwr5";
$norm2 = 10**$cpwr2; $norm3 = 10**$cpwr3; $norm4 = 10**$cpwr4; $norm5 = 10**$cpwr5;
print "\n $cpwr2 $cpwr3 $cpwr4 $cpwr5 $norm2 $norm3 $norm4 $norm5 \n\n";

$title_strain  = sprintf("Estimated strain, 10\@+%1i\@+ yr\@+-1\@+",$cpwr3);   # %2.2i for 07
$title_rot    = sprintf("Estimated rotation, 10\@+%1i\@+ yr\@+-1\@+",$cpwr4);
$title_dilat  = sprintf("Estimated dilatation, 10\@+%1i\@+ yr\@+-1\@+",$cpwr2);

$Bscale_strain  = sprintf("-B%2.2f:\"Strain rate, 10\@+%2.2i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);
$Bscale_rot    = sprintf("-B%2.2f:\"Rotation rate, 10\@+%2.2i\@+ rad/yr\": -Ef10p",$cmaxrot*0.25,$cpwr4);
$Bscale_dilat  = sprintf("-B%2.2f:\"Dilatation rate, 10\@+%2.2i\@+ yr\@+-1\@+\": -E10p",$cmaxdilat*0.5,$cpwr2);
$unit_rot      = sprintf("10\@+%i\@+ rad/yr",$cpwr4);

# subplotting specifications
$dX = $wid + 0.5;
$dX2 = -$dX;
$shift1 = "-X$dX";
$shift2 = "-X$dX2 -Y-3.4";
$shift3 = $shift1;

$x_title = 0.5; $z_title = 0.9; $z_title2 = 0.85;
$J1 = "-JM${wid}i";

# plot title
$J_title = "-JX${wid}";
$R_title = "-R0/1/0/1";

$fsize_title = 12;
$fsize_title0 = 12;

# open CSH file
if($k==$kmin) {
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize3  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen D_FORMAT \%g\n";
}

#-----------------------------------------
# make color files

# surface velocities
$cran = $cmax - $cmin; $dc = $cran/100;
$Tvmag = "-T$cmin/$cmax/$dc";

$colorbar = "seis";
$colorbar = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";

$cptvmag = "color1.cpt";
print CSH "makecpt -Crainbow $Tvmag -D > $cptvmag\n";

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

$cptdilat = "color2.cpt";
print CSH "makecpt -C$colorbar $Tdilat -D > $cptdilat\n";

$cptstrain = "color3.cpt";
print CSH "makecpt -C$colorbar $Tstrain -D > $cptstrain\n";

$cptrot = "color4.cpt";
print CSH "makecpt -C$colorbar $Trot -D > $cptrot\n";

#-----------------------------------------

# color bar
#$Dlen = 0.6*$wid; $Dx = $wid/2;
$Dlen = 0.5*$wid; $Dx = $Dlen/2;
$Dy = -0.35;
$Dscale = "-D$Dx/$Dy/$Dlen/0.1h";

#$B0 = "-Ba${xtick1}f${xtick2}d::";
$B0 = sprintf("-Ba%3.3ff%3.3fd:\" \":/a%3.3ff%3.3fd:\" \"::.\" \":",$xtick1,$xtick2,$ytick1,$ytick2);

# KEY: commands for interpolation and masking
$interp = "-I0.1";
$interp_surf = "-I0.04 -T0 -S0";
#if($iregion==1) {$interp_surf = "-I0.2 -T0 -S0";}
$mask_info = "-I0.1 -G200 -S0.2";
$grdfile = "temp.grd";
$mask_file = "${vel_dir0}/${name}_masked_pts.dat";
if (not -f ${mask_file}) { die("Check if ${mask_file} exist or not\n") }

$pname = $name;

#-----------------------------------------

$fname = "foursub_${ntag}";
$psfile = "${fname}.ps"; $jpgfile = "${fname}.jpg";

$title0  = "(a)  Surface velocity, mm/yr";
$Bscale  = "-B$ctick:\"  \": -Ef10p";
#$Dlen = 1; $Dx = 1; $Dscale0 = "-D1.2/$Dy/$Dlen/0.1h";
$Dscale0 = $Dscale;
$B = $B0.$Bopts[7];

print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
if ($icolor==1) {print CSH "grdimage $grdfile $R1 $J1 -C$cptvmag -T -K -O -V >> $psfile\n";}
if ($imask==1) {
  print CSH "echo MASK DATA FILE: ${mask_file}\n";
  print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
  print CSH "psmask -C -K -O -V >> $psfile\n";
}
print CSH "psscale -C$cptvmag $Dscale0 $Bscale -K -O -V >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_w -K -V -O -P >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_w -K -V -O -P >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

print "VELOCITY DATA FILE: ${velocity_vec_dat}\n";

#$gpstemp = $gps_vec;
$gpstemp = gpstemp1;
print CSH "awk '{print \$1,\$2,\$3,\$4,\$5,\$6}' ${velocity_vec_dat} > $gpstemp\n";

# plot black arrows on top of white arrows
$vec_info0 = "-N -L -Se${vscale}/${vec_conf}/0";
#$vec_info = "$vec_info0 -A1.5p/3p/2p -W1.5p,255/255/255";
#print CSH "psvelo $gpstemp $J1 $R1 $vec_info -K -V -O >> $psfile\n";
$vec_info = "$vec_info0 -A0.25p/3p/1.5p -W0.25p,0/0/0";
print CSH "psvelo $gpstemp $J1 $R1 $vec_info -K -V -O >> $psfile\n";

# plot spline centers
#print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

# plot v-field scale
print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0.25 $fsize2 0 $fontno LM ${vscale_text}\nEOF\n";
print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

# plot pole of rotational field
#if($idata >= 20 && $idata < 30) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}

#if($idata==1) {    # SAFOD hole
#print CSH "psxy $J1 $R1 $poleinfo2 -K -O -V >>$psfile<<EOF\n$plon $plat\nEOF\n";
#}

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title2 $fsize_title0 0 $fontno CM $title0\nEOF\n";

#-----------------------------------------

$title0  = sprintf("(b)  Dilatation rate, 10\@+%i\@+ yr\@+-1\@+",$cpwr2);
$Bscale  = sprintf("-B%2.2f:\" \": -E10p",$cmaxdilat*0.5);
$B = $B0.$Bopts[4];

print CSH "psbasemap $J1 $R1 $B -K -V -O $shift1 >> $psfile\n";
if ($icolor == 1) {
  print CSH "echo STRAIN DATA FILE: ${strain_mag}\n";
  print CSH "awk '{print \$1,\$2,\$4}' $strain_mag > ${name}_dilat.xyz\n";
  print CSH "awk '{print \$1,\$2,\$4/$norm2}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptdilat -T -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo MASK DATA FILE: ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
    print CSH "psmask -C -K -O -V >> $psfile\n";
  }
}
print CSH "psscale -C$cptdilat $Dscale $Bscale -K -O -V >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

# highlight the small volcanic source
#if($k==3) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n";}

#if($idata==1) {    # SAFOD hole
#print CSH "psxy $J1 $R1 $poleinfo2 -K -O -V >>$psfile<<EOF\n$plon $plat\nEOF\n";
#}

#if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}
#$print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title2 $fsize_title0 0 $fontno CM $title0\nEOF\n";

#-----------------------------

$title0  = sprintf("(c)  Strain rate, 10\@+%i\@+ yr\@+-1\@+",$cpwr3);
$Bscale  = sprintf("-B%2.2f:\"  \": -Ef10p",$cmaxstrain*0.25);
$B = $B0.$Bopts[7];

# plot basemap
print CSH "psbasemap $J1 $R1 $B -K -V -O $shift2 >> $psfile\n";
#print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | pscontour $R1 $J1 -A- -C$cptstrain -I -O -K -V >> $psfile\n";

if($icolor==1) {
  print CSH "echo STRAIN DATA FILE: ${strain_mag}\n";
  print CSH "awk '{print \$1,\$2,\$5}' $strain_mag > ${name}_strain.xyz\n";
  print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -T -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo MASK DATA FILE: ${mask_file}\n";
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
if ($iplate==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

# plot pole of rotational field
#if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}

# highlight the small volcanic source
#if($k==3) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n";}

if ($idata==1) {		# SAFOD hole
  #print CSH "psxy $J1 $R1 $poleinfo2 -K -O -V >>$psfile<<EOF\n$plon $plat\nEOF\n";
  $platmin = 35.5; $platmax = 36.1; $plonmin = -121.2; $plonmax = -119.8;
  print CSH "psxy $J1 $R1 -W2p,0/0/0 -A -O -K -V <<EOF>>$psfile\n"; 
  print CSH "$plonmin $platmin\n $plonmax $platmin\n $plonmax $platmax\n $plonmin $platmax\n $plonmin $platmin\nEOF\n"; 
}

#print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title2 $fsize_title0 0 $fontno CM $title0\nEOF\n";

#==============================

$title0 = sprintf("(d)  Rotation rate, |w|, 10\@+%i\@+ rad/yr",$cpwr4);
$Bscale  = sprintf("-B%2.2f:\"  \": -Ef10p",$cmaxrot*0.25);
$B = $B0.$Bopts[4];

print CSH "psbasemap $J1 $R1 $B -K -V -O $shift3 >> $psfile\n";
if ($icolor==1) {

  # plot magnitude of w (3D) or its vertical or horizontal components
  print CSH "echo STRAIN DATA FILE: ${strain_mag}\n";
  print CSH "awk '{print \$1,\$2,\$6/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  #print CSH "awk '{print \$1,\$2,sqrt(\$11*\$11+\$12*\$12+\$13*\$13)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  #print CSH "awk '{print \$1,\$2,sqrt(\$11*\$11)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  #print CSH "awk '{print \$1,\$2,sqrt(\$12*\$12+\$13*\$13)/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";

  print CSH "grdimage $grdfile $R1 $J1 -C$cptrot -T -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo MASK DATA FILE: ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
    print CSH "psmask -C -K -O -V >> $psfile\n";
  }
}
print CSH "psscale -C$cptrot $Dscale $Bscale -K -O -V >> $psfile\n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

#-----------------------------
#if ( ($idata >= 20 && $idata < 30) || ($idata >= 60 && $idata < 80) ) {
if($ieuler==1) {
  $pinfo      = "-Sc4p -G255/255/255 -W0.5p,0/0/0";
  $pinfo_anti = "-Sc4p -G255/0/0 -W0.5p,0/0/0";
  $seis_info3 = "-Scp -W0.5p/150/150/150";

  # make colorpoint file
  $cptfile = "color.cpt"; $cmin = -0.01; $cmax = 0.01; $dc = ($cmax-$cmin)/2;
  $T = "-T$cmin/$cmax/$dc -D -Z";
  #$colorbar = "-Cpolar";
  $colorbar = "-Cgray -I";
  print CSH "makecpt $colorbar $T > $cptfile\n";

  # KEY: plot euler poles
  print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_anti_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";

  if ( ($idata >= 20 && $idata < 30) || ($idata >= 60 && $idata < 80) ) {
    $epoles = "${vel_dir0}/${name}_epole_mean_points.dat";
    if (not -f $epoles) { die("Check if $epoles exist or not\n") }
    print CSH "psxy $epoles $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile\n";
  }

  # plot key for euler vectors
  #-----------------------------
  #$origin_arrow = "-Xa0.5 -Ya-1.5";

  $xdots = 0.6*$wid;
  $ydots = -1.0;
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
  $dy = 0.2;
  for ($k = 0; $k < 5; $k = $k+1) {
    $ytxt[$k] = $k*$dy;
  }
  $x1a = 0;
  $x2a = $x1a + 1.20;
  $x1m = ($x1a+$x2a)/2;
  $y1 = -0.8*$dy;

  # text labels
  $textinfo = "-N -Dj -C8p -W255/255/255,O1p,255/255/255";
  print CSH "pstext -JX1 -R0/1/0/1 $textinfo -K -O -V ${origin_dots} >>$psfile<<EOF
$x1m $ytxt[4] 10 0 $fontno CM pos $unit_rot neg
EOF\n";
  print CSH "pstext -JX1 -R0/1/0/1 -N -K -O -V ${origin_dots} >>$psfile<<EOF
$x1m $ytxt[3] 10 0 $fontno CM $sv3
$x1m $ytxt[2] 10 0 $fontno CM $sv2
$x1m $ytxt[1] 10 0 $fontno CM $sv1
EOF\n";

  # dots
  print CSH "psxy -JX1 -R0/1/0/1 -C$cptfile -N ${seis_info3} ${origin_dots} -K -O -V >>$psfile<<EOF
$x1a $ytxt[3] $vec_ref[2] $dsize[2]
$x1a $ytxt[2] $vec_ref[1] $dsize[1]
$x1a $ytxt[1] $vec_ref[0] $dsize[0]
$x2a $ytxt[3] -$vec_ref[2] $dsize[2]
$x2a $ytxt[2] -$vec_ref[1] $dsize[1]
$x2a $ytxt[1] -$vec_ref[0] $dsize[0]
EOF\n";


  #-----------------------------

}
#-----------------------------

# highlight the small volcanic source
#if($k==3) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n";}

#if($idata==1) {    # SAFOD hole
#print CSH "psxy $J1 $R1 $poleinfo2 -K -O -V >>$psfile<<EOF\n$plon $plat\nEOF\n";
#}

#print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title2 $fsize_title0 0 $fontno CM $title0\nEOF\n";

#--------------------------------

# horizontal line
print CSH "psxy -JX10 -R0/1/0/1 -W1p -K -O -V -Xa-3.6 -Ya2.65 >>$psfile<<EOF\n0 0\n0.67 0\nEOF\n";

print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}   # END LOOP

#==================================================

close (CSH);
system("csh -f $cshfile");
#print "convert $psfile $jpgfile\n";
#system("convert $psfile $jpgfile");
#if($ipdf==1) {system("ps2pdf $psfile")};
#if($ixv==1) {system("xv $jpgfile &")};

#==================================================
