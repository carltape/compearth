#!/usr/bin/perl -w

#==========================================================
#
#  syn_four_strain.pl
#  Carl Tape
#  23-Aug-2007
#  
#  This script makes a 4 x 3 figure containing dilatation, strain, and rotation
#  for four synthetic examples.  Instead of using this figure, we now use
#  syn_one.pl, which makes four separate figures.
#  
#==========================================================

$cshfile = "syn_four_strain.csh";

$icolor = 1;    # ccc

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
$poleinfo1 = "-Sa20p -G0/255/255 -W0.75p,0/0/0";
$circleinfo = "-Sc25p -W1.0p,0/0/0,--";

# velocity, strain rate, etc
$vel_dir0   = "${plate_dir}/surface_velocities";

#----------------------------------

# file name, ticks, range for velocities

$iregion = 3;      # region (1=west_us, 2=cal, 3=socal, 6=cascadia)
$idata = 71;        # choose GPS dataset (1 = NASA REASON; 2 = CCMM; 1X = strike-slip; 2X = rotational)
$ndim = 2;         # ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 1;         # min grid order used
$qmax = 7;         # max grid order used
$basis = 2;        # 1 for splines; 2 for wavelets
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$imask  = 1;       # plot mask
$ieuler = 1;       # plot euler poles

$ifault = 0;       # plot faults
$iplate = 0;       # plot plate boundaries

   $itag = "socal";

   $x_title = 0.1; $z_title = -0.1;
   $xtick1 = 1; $xtick2 = 0.25; $wid = 2.0;
   $origin = "-X0.75 -Y8.5";
   $xlegend = 0.25; $ylegend = -1.0;

$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";
$evec_ref = 1.0; $escale = 8;    # default values
$fac_min = 1.0; $fac_max = 1.0;  # default values

$fname = "syn_four";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

$t = -1;
@tlabs = ("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)");

#=========================================
# KEY LOOP
#=========================================
for ($k = 1; $k <= 4; $k = $k+1) {

if ($k == 1 ) {
   $idata = 61; $ndim = 2; $imask  = 0; $qmin = 1; $qmax = 7;
   $ib = 6;

} elsif ($k == 2) {
   $idata = 71; $ndim = 2; $imask  = 0; $qmin = 1; $qmax = 7;
   $ib = 2;

} elsif ($k == 3) {
   $idata = 83; $ndim = 3; $imask  = 1; $qmin = 1; $qmax = 7;
   $ib = 6;

} elsif ($k == 4) {
   $idata = 13; $ndim = 2; $imask  = 1; $qmin = 1; $qmax = 7;
   $ib = 2;
}


# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
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

# SYNTHETIC velocity field
if ($idata >= 10) {

  $gps_pts = "${gps_dir}/synthetic/syn_vfield_${dlab}_points.dat";
  $gps_vec = "${gps_dir}/synthetic/syn_vfield_${dlab}_psvelo.dat"; 

  # strike-slip field
  if ($idata >= 10 and $idata < 20) {
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

  # rotational field
  if ($idata >= 20 and $idata < 30) {
    $igc = 0; $ifault = 0; $iplate = 0;
    $vec_ref = 10; $vec_conf = 0.9; $vscale = 0.06;
    $cmin = 0; $cmax = 8; $ctick = 2;
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
    $vec_ref = 1000; $vec_conf = 0.9; $vscale = 0.001;
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
    $vec_ref = 1000; $vec_conf = 0.9; $vscale = 0.0005;
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
    $vec_ref = 10; $vec_conf = 0.9; $vscale = 0.06;
    $cmin = 0; $cmax = 5; $ctick = 1;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 2.0; $escale = 5;
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
    $vec_ref = 10; $vec_conf = 0.9; $vscale = 0.06;
    $cmin = 0; $cmax = 5; $ctick = 1;
    $vumin = -1000; $vumax = 0; $vsmin = -35; $vsmax = 0; $vemin = -35; $vemax = -$vemin;
    $evec_ref = 2.0; $escale = 5;
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
    $vec_ref = 40; $vec_conf = 0.9; $vscale = 0.006;
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

#--------------------------------------------

# plotting specifications
$fsize0 = "16";
$fsize1 = "14";
$fsize2 = "11";
$fsize3 = "8";
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

$name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}";

print "\n $name \n";

if($igc==1){
  $gc_boundary  = "/home/carltape/gmt/gps_data/synthetic/gps_gc_${dlab}.dat";
  #$gc_boundary  = "/home/carltape/gmt/gps_data/synthetic/socal_gps_SAFplanar_gc$glab.dat";
  if (not -f $gc_boundary) { die("Check if $gc_boundary exist or not\n") }
}

# components and magnitudes of v-field
$strain_mag        = "${vel_dir0}/misc/${name}_spline_strain.dat";
$velocity_vec      = "${vel_dir0}/misc/${name}_spline_vec.dat";
$bounds            = "${vel_dir0}/misc/${name}_spline_bounds.dat";
$colors            = "${vel_dir0}/misc/${name}_spline_colors.dat";
$spline_centers    = "${vel_dir0}/misc/${name}_spline_gridpoints.dat";
$euler_poles       = "${vel_dir0}/misc/${name}_euler_vector_psxy.dat";
$euler_anti_poles  = "${vel_dir0}/misc/${name}_euler_anti_vector_psxy.dat";
$euler_poles_scale = "${vel_dir0}/misc/${name}_euler_vector_scale_psxy.dat";

if (not -f $strain_mag)     { die("Check if $strain_mag exist or not\n") }
if (not -f $velocity_vec)   { die("Check if ${velocity_vec} exist or not\n") }
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
  $latmin = 32; $latmax = 38; $lonmin = -122; $lonmax = -114;
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
$Bscale_rot    = sprintf("-B%2.2f:\"Rotation rate, 10\@+%2.2i\@+ yr\@+-1\@+\": -Ef10p",$cmaxrot*0.25,$cpwr4);
$Bscale_dilat  = sprintf("-B%2.2f:\"Dilatation rate, 10\@+%2.2i\@+ yr\@+-1\@+\": -E10p",$cmaxdilat*0.5,$cpwr2);
$unit_rot      = sprintf("10\@+%2.2i\@+ yr\@+-1\@+",$cpwr4);

# subplotting specifications
$hwid = $wid/2;
$xfac = 1.20;
$yfac = 1.65;
$dX = $xfac*$wid; $dY = 0; $shift1 = "-X$dX -Y$dY";
$dX = -$dX; $dY = -$yfac*$wid; $shift2 = "-X$dX -Y$dY";

$dX = $wid + 0.3;
$dX2 = -2*$dX;
$shift = "-X$dX";
$shift0 = "-X$dX2 -Y-2.35";

$J1 = "-JM${wid}i";

# plot title
$J_title = "-JX${wid}";
$R_title = "-R0/1/0/1";

$fsize_title = 12;
$fsize_title0 = 12;

# open CSH file
if ($k == 1) {
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain MEASURE_UNIT inch PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize3  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";
}

#-----------------------------------------
# make color files

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

$colorbar = "rainbow";
$c_romania = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";

$colorbar = $c_romania;
$colorbar = "seis -I";
# -I to invert the color scale

$colorbar = "seis";
$cptdilat = "color2.cpt";
print CSH "makecpt -C$colorbar $Tdilat -D > $cptdilat\n";

$colorbar = "seis";
$cptstrain = "color3.cpt";
print CSH "makecpt -C$colorbar $Tstrain -D > $cptstrain\n";

$colorbar = "seis";
$cptrot = "color4.cpt";
print CSH "makecpt -C$colorbar $Trot -D > $cptrot\n";

#-----------------------------------------

# color bar
$Dlen = 0.6*$wid;;
$Dx = $Dlen/2;
$Dy = -0.1;
$Dscale = "-D1.2/$Dy/$Dlen/0.1h";

$B0 = "-Ba${xtick1}f${xtick2}d::";

$B = $B0.$Bopts[15];

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

$fname = "syn_four_strain";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

$t = $t+1;
$title0  = sprintf("Dilatation rate, 10\@+%i\@+ yr\@+-1\@+",$cpwr2);
$title  = $tlabs[$t];
$Bscale  = sprintf("-B%2.2f:\" \": -E10p",$cmaxdilat*0.5);

if ($k == 1) {
  print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
} else {
  print CSH "psbasemap $J1 $R1 $B -K -O -V $shift0 >> $psfile\n";
}

if ($icolor == 1) {
  print CSH "awk '{print \$1,\$2,\$4/$norm2}' $strain_mag | surface -G$grdfile ${interp_surf} $R1 \n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptdilat -T -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile\n";
    print CSH "psmask -C -K -O -V >> $psfile\n";
  }
}
print CSH "psscale -C$cptdilat $Dscale $Bscale -K -O -V >> $psfile \n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

# highlight the small volcanic source
#if($k==3) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n";}

#if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}
print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
if($k==1) {print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n 0.5 1.03 $fsize_title0 0 $fontno CM $title0 \nEOF\n";}

#-----------------------------

$t = $t+1;
$title  = $tlabs[$t];
$title0  = sprintf("Strain rate, 10\@+%i\@+ yr\@+-1\@+",$cpwr3);
$Bscale  = sprintf("-B%2.2f:\"  \": -Ef10p",$cmaxstrain*0.25);

# plot basemap
print CSH "psbasemap $J1 $R1 $B -K -V -O $shift >> $psfile\n";
#print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | pscontour $R1 $J1 -A- -C$cptstrain -I -O -K -V >> $psfile\n";

if($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1 \n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -T -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile \n";
    print CSH "psmask -C -K -O -V >> $psfile \n";
  }
}

# color bar
#print CSH "gmtset D_FORMAT \%.2e \n";
print CSH "psscale -C$cptstrain $Dscale $Bscale -K -O -V >> $psfile \n";
#print CSH "gmtset D_FORMAT \%lg \n";

print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

# plot pole of rotational field
#if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}

# highlight the small volcanic source
#if($k==3) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n";}

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
if($k==1) {print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n 0.5 1.03 $fsize_title0 0 $fontno CM $title0 \nEOF\n";}

#==============================

$t = $t+1;
$title  = $tlabs[$t];
$title0 = sprintf("Rotation rate, 10\@+%i\@+ yr\@+-1\@+",$cpwr4);
$Bscale  = sprintf("-B%2.2f:\"  \": -Ef10p",$cmaxrot*0.25);

print CSH "psbasemap $J1 $R1 $B -K -V -O $shift >> $psfile\n";
if ($icolor==1) {
  print CSH "awk '{print \$1,\$2,\$6/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1 \n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptrot -T -K -O -V >> $psfile\n";
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile \n";
    print CSH "psmask -C -K -O -V >> $psfile \n";
  }
}
print CSH "psscale -C$cptrot $Dscale $Bscale -K -O -V >> $psfile \n";
print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

# plot spline centers
#print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

# highlight the small volcanic source
#if($k==3) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n";}

print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title \nEOF\n";
if($k==1) {print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n 0.5 1.03 $fsize_title0 0 $fontno CM $title0 \nEOF\n";}

}   # END LOOP

#--------------------------------

print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk \nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "xv $jpgfile &\n"}

#==================================================

close (CSH);
system("csh -f $cshfile");
#print "convert $psfile $jpgfile \n";
#system("convert $psfile $jpgfile");
#if($ipdf==1) {system("ps2pdf $psfile")};
#if($ixv==1) {system("xv $jpgfile &")};

#==================================================
