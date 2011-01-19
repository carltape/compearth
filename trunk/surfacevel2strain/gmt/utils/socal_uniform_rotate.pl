#!/usr/bin/perl -w

#==========================================================
#
#  socal_uniform_rotate.pl
#  Carl Tape
#  23-Aug-2007
#  
#  This script inputs a surface velocity field and strain scalar fields
#  computed in test_platemodel2strain.m and outputs a figure.
#  
#==========================================================

$cshfile = "socal_uniform_rotate.csh";

$icolor = 1;

# plates and faults
$plate_dir = "/home/carltape/gmt/plates";
$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
if (not -f ${plate_file}) {die("Check if ${plate_file} exist or not\n")}
$kcf_file     = "/home/carltape/gmt/faults/kcf.xy";
$fault_file   = "/home/carltape/gmt/faults/jennings.xy";
if (not -f $fault_file) {die("Check if $fault_file exist or not\n")}

$fault_info_k = "-M -W1.5p,0/0/0";
$fault_info_r = "-M -W1.5p,255/0/0";

# plot euler poles and anti-poles
$pinfo      = "-Sc4p -G255/255/255 -W0.5p,0/0/0";
$pinfo_anti = "-Sc4p -G255/0/0 -W0.5p,0/0/0";
$seis_info3 = "-Scp -W0.5p/255/255/255";

$sinfo = "-Sc4p -G255 -W0.5p,0/0/0";              # spline centers
$poleinfo1 = "-Sa16p -G0/255/255 -W0.5p,0/0/0";   # euler pole star
$poleinfo2 = "-Sa16p -G255 -W0.5p,0/0/0";   # euler pole star

# velocity, strain rate, etc
$vel_dir0   = "${plate_dir}/surface_velocities";

#----------------------------------

# file name, ticks, range for velocities

$iregion = 3;			# region (1=west_us, 2=cal, 3=socal, 6=cascadia)
$idata = 21;			# choose GPS dataset (1 = NASA REASON; 2 = CCMM; 1X = strike-slip; 2X = rotational)
$ndim = 2;			# ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 0;			# min grid order used
$qmax = 7;			# max grid order used
$basis = 2;			# 1 for splines; 2 for wavelets
$iLmat = 1;			# sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$imask  = 0;			# plot mask
$ieuler = 1;			# plot euler poles

if ($iregion == 3) {
  $itag = "socal";
  $xtick1 = 1; $xtick2 = 0.25; $wid = 6.0;
  #$origin = "-X1.0 -Y1.0";
  #$xlegend = 0.25; $ylegend = 1.3;
  $origin = "-X1.0 -Y2.5";
  $xlegend = 0.25; $ylegend = -1.0;
  $ifault = 1;
}
$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";
$evec_ref = 1.0;

# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
$dlab = sprintf("d%2.2i",$idata);

$tag = "NASA REASON data set (cGPS)";
$gps_pts = "${gps_dir}/US/reason_fixed_NAM_subset_points.dat";
$gps_vec = "${gps_dir}/US/reason_fixed_NAM_subset_psvelo.dat";
$igc = 0; $vec_ref = 50; $vec_conf = 0.9; $vscale = 0.01;
$cmin = 0; $cmax = 40; $ctick = 10;
$vumin = -2; $vumax = 2; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;

# synthetic velocity field

#$gps_pts = "${gps_dir}/synthetic/syn_vfield_${itag}_${dlab}_points.dat";
#$gps_vec = "${gps_dir}/synthetic/syn_vfield_${itag}_${dlab}_psvelo.dat"; 

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

$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);
$vscale_text1 = sprintf("%.0f mm/yr", $vec_ref);
$vscale_text2 = sprintf("(%.0f%%)",$vec_conf*100);

#--------------------------------------------

# plotting specifications
$fsize0 = "20";
$fsize1 = "18";
$fsize2 = "14";
$fsize3 = "6";
$fontno = "1";			# 1 or 4
$tick   = "0.2c";
$fpen   = "1p";
$tpen   = "1p";

$phi = "\@~\146\@~";
$theta = "\@~\161\@~";

# plotting specificiations

# A : smallest feature plotted, in km^2; D : resolution
$coast_res     = "-A500 -Df";
$coast_infoK   = "$coast_res -W1.0p,0/0/0 -Na/1.0p";
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

$ivel = 0;
$ishear = 0;
$irotate = 0;
$idilate = 0;
$imaxlam = 0;
$imap = 0;

#==================================================

$glab = sprintf("%2.2i",$igc);
$qlab = sprintf("q%2.2i_q%2.2i",$qmin,$qmax);
$blab = sprintf("b%1.1i",$basis);
$nlab = sprintf("%1iD",$ndim);
$slab = sprintf("s%1i",$iLmat);
#if ($idata==0) {       $slab = "Lfull";
#} elsif ($idata==1) {  $slab = "Lelastic";
#} elsif ($idata==2) {  $slab = "Lviscous";
#} else {
#  die("\nExit: invalid iLmat option\n\n");
#}

$name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}";

print "\n $name \n";

# bounds and color limits
$bounds            = "${vel_dir0}/misc/${name}_spline_bounds.dat";
$colors            = "${vel_dir0}/misc/${name}_spline_colors.dat";
if (not -f $bounds) {die("Check if $bounds exist or not\n")}
if (not -f $colors) {die("Check if $colors exist or not\n")}

# get bounds of the region
open(IN,"$bounds"); @lines = <IN>;
($lonmin,$lonmax,$latmin,$latmax) = split(" ",$lines[0]);

# modify the bounds
$latmin = 32;
#if($iregion == 3) {
#  $latmin = 31; $latmax = 38; $lonmin = -122; $lonmax = -114;
#}

$R1 = "-R$lonmin/$lonmax/$latmin/$latmax";
print "\n bounds are $R1 \n\n";
#die("testing");

# get color limits for the plots: dilatation, shear, rotation, maxlam
open(IN,"$colors"); @lines = <IN>;
($cmin2,$cmaxdilat,$cpwr2) = split(" ",$lines[0]); # dilatation
($cmin3,$cmaxshear,$cpwr3) = split(" ",$lines[1]); # shear
($cmin4,$cmaxrot,$cpwr4) = split(" ",$lines[2]); # rotation
($cmin5,$cmax5,$cpwr5) = split(" ",$lines[3]); # maxlam
#$norm2 = "1e$cpwr2"; $norm3 = "1e$cpwr3"; $norm4 = "1e$cpwr4"; $norm5 = "1e$cpwr5";
$norm2 = 10**$cpwr2; $norm3 = 10**$cpwr3; $norm4 = 10**$cpwr4; $norm5 = 10**$cpwr5;
print "\n $cpwr2 $cpwr3 $cpwr4 $cpwr5 $norm2 $norm3 $norm4 $norm5 \n\n";

$title_shear  = sprintf("Estimated shear, 10\@+%1i\@+ yr\@+-1\@+",$cpwr3); # %2.2i for 07
$title_rot    = sprintf("Estimated rotation, 10\@+%1i\@+ yr\@+-1\@+",$cpwr4);
$title_dilat  = sprintf("Estimated dilatation, 10\@+%1i\@+ yr\@+-1\@+",$cpwr2);

$title_up    = "Estimated Vr (mm/yr)";
$title_south = "Estimated V$theta (mm/yr)";
$title_east  = "Estimated V$phi (mm/yr)";

$Bscale_shear  = sprintf("-B%2.2f:\"Shear strain rate, 10\@+%i\@+ yr\@+-1\@+\": -E10p",$cmaxshear*0.25,$cpwr3);
$Bscale_rot    = sprintf("-B%2.2f:\"Rotation rate, 10\@+%i\@+ yr\@+-1\@+\": -E10p",$cmaxrot*0.25,$cpwr4);
$Bscale_dilat  = sprintf("-B%2.2f:\"Dilatation rate, 10\@+%i\@+ yr\@+-1\@+\": -E10p",$cmaxdilat*0.5,$cpwr2);
$unit_rot      = sprintf("10\@+%i\@+ yr\@+-1\@+",$cpwr4);

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

print "\n @vtick \n";

# dilatation
#$cmaxdilat = $cmaxdilat/2; $cmin2 = $cmin2/2;
$cran = ($cmaxdilat - $cmin2);
$dc = $cran/100;
$Tdilat = "-T$cmin2/$cmaxdilat/$dc";

# shear
$cran = $cmaxshear - $cmin3;
$dc = $cran/100;
$Tshear = "-T$cmin3/$cmaxshear/$dc";

# rotation
$cran = $cmaxrot - $cmin4;
$dc = $cran/100;
$Trot = "-T$cmin4/$cmaxrot/$dc";

# maxlam
$cran = $cmax5 - $cmin5;
$dc = $cran/100;
$T5 = "-T$cmin5/$cmax5/$dc";

# maxlam -- log10 scale
$T6 = "-T0.01/$cmax5/3";

$colorbar = "rainbow";
$c_romania = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";

# v-field components
$cptup = "color_up.cpt"; print COLOR "makecpt -C${c_romania} $Tup -D > $cptup \n";
$cptsouth = "color_south.cpt"; print COLOR "makecpt -C${c_romania} $Tsouth -D > $cptsouth \n";
$cpteast = "color_east.cpt"; print COLOR "makecpt -C${c_romania} $Teast -D > $cpteast \n";

$cptvmag = "color1.cpt";
print COLOR "makecpt -C$colorbar $Tvmag -D > $cptvmag \n";
#print COLOR "sed 's/^B.*/B       255   255    255  /' temp1  >  temp2\n";
#print COLOR "sed 's/^F.*/F       255     0      0  /' temp2 > $cptvmag\n";

if (1==1) {
  $colorbar = $c_romania;
  #$cback = "  0     0    255";
  #$cfor  = "255     0      0";

} else {
  $colorbar = "no_green";
  #$cback = "170     0      0";
  #$cfor  = "  0     0    200";
}

# 170     0      0  to 255   255    255

$cptdilat = "color2.cpt";
print COLOR "makecpt -C$colorbar $Tdilat -D > $cptdilat\n";
#print COLOR "makecpt -C$colorbar $Tdilat -D > temp1\n";
#print COLOR "sed 's/^B.*/B       $cback  /' temp1  >  temp2\n";
#print COLOR "sed 's/^F.*/F       $cfor  /' temp2 > $cptdilat\n";

$cptshear = "color3.cpt";
print COLOR "makecpt -C$colorbar $Tshear -D > $cptshear\n";

$cptrot = "color4.cpt";
print COLOR "makecpt -C$colorbar $Trot -D > $cptrot\n";

$cptfile5 = "color5.cpt";
print COLOR "makecpt -C$colorbar $T5 -D > $cptfile5\n";

#$cptfile6 = "color6.cpt";
#print COLOR "makecpt -C$colorbar $T6 -Qo -D > $cptfile6\n";

close (COLOR);
system("csh -f ${cshfile_color}");

#-----------------------------------------

# KEY: commands for interpolation and masking
$interp = "-I0.1";
$interp_surf = "-I0.04 -T0 -S0";
#if($iregion==1) {$interp_surf = "-I0.2 -T0 -S0";}
$mask_info = "-I0.1 -G200 -S0.2";
$grdfile = "temp.grd";

#==============================
# estimated velocity fields

@labs = ("A","B","C","D");
$fault_info_k    = "-M -W1.0p,0/0/0";
$fault_info_r    = "-M -W1.0p,255/0/0";

if ($iregion == 3) {
  $wid = 2.5; $vscale = 0.01;
  $origin = "-X1.0 -Y7.25";
  $xfac = 1.20; $yfac = 1.1;
  $dX1 = $xfac*$wid; $dY1 = 0; 
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $xlegend = 0.25; $ylegend = -1.0;
  $xtick1 = 2; $xtick2 = 1;
  $x_title = 0.; $z_title = 0.92;

  $Dlen = 1;
  $Dx = $wid + 0.25;
  $Dy = $Dlen/2;
  $Dscale = "-D$Dx/$Dy/$Dlen/0.12";

  $origin_arrow = "-Xa$Dx -Ya2";
  $origin_dots = "-Xa$Dx -Ya1.5";

}
$B0 = "-Ba${xtick1}f${xtick2}d::";

$hwid = $wid/2;
$fsize_title = 10;
$J1 = "-JM${wid}i";
$J_title = "-JX${wid}";

$dX2 = -$dX1; $dY2 = -$yfac*$wid;
$shift1 = "-X$dX1 -Y$dY1";
$shift2 = "-X$dX2 -Y$dY2";

$fname = "socal_uniform_rotate";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";
$title = "Rotation rate from $tag";

# SUB-PLOTTING FORMAT (idata)
#
#    21   23
#    21   23
#    20   22
#

# open CSH file
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter PLOT_DEGREE_FORMAT D TICK_LENGTH 0.15c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 10 FRAME_PEN 1p TICK_PEN 1p\n";

@idata_vel = (21,23);
@Bins = (7,4);
#$Bscale  = "-B$ctick:\"  \": -Ef5p";
$Bscale  = "-B$ctick:\"  \":/:\"mm/yr\": -Ef5p";
@titles = ("(a)  Synthetic velocity field (mm/yr)","(c)  Synthetic velocity field (mm/yr)");
#@titles = ("(a)","(b)");

# plot the velocity field
for ($i = 1; $i <= 2; $i = $i+1) {

  $idata = $idata_vel[$i-1];
  $dlab = sprintf("d%2.2i",$idata);
  $gps_pts = "${gps_dir}/synthetic/syn_vfield_${dlab}_points.dat";
  $gps_vec = "${gps_dir}/synthetic/syn_vfield_${dlab}_psvelo.dat"; 
  if (not -f $gps_pts) {die("Check if $gps_pts exist or not\n")}
  if (not -f $gps_vec) {die("Check if $gps_vec exist or not\n")}
  $name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}";

  $mask_file = "${vel_dir0}/misc/${name}_masked_pts.dat";
  if (not -f ${mask_file}) {die("Check if ${mask_file} exist or not\n")}
  if ($i == 2) {$imask = 1}

  # components and magnitudes of v-field
  $strain_mag        = "${vel_dir0}/misc/${name}_spline_strain.dat";
  $velocity_vec      = "${vel_dir0}/misc/${name}_spline_vec.dat";
  $spline_centers    = "${vel_dir0}/misc/${name}_spline_gridpoints.dat";
  if (not -f ${strain_mag}) {die("Check if ${strain_mag} exist or not\n")}
  if (not -f ${velocity_vec}) {die("Check if ${velocity_vec} exist or not\n")}
  if (not -f ${spline_centers}) {die("Check if ${spline_centers} exist or not\n")}

  $iB = $Bins[$i-1];
  $B = $B0.$Bopts[$iB];

  if ($i == 1) {
    print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
  } else {
    print CSH "psbasemap $J1 $R1 $B -K -O -V $shift1 >> $psfile\n";
  }
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1 \n";
    print CSH "grdimage $grdfile $R1 $J1 -C$cptvmag -T -K -O -V >> $psfile\n";
  }
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile \n";
    print CSH "psmask -C -K -O -V >> $psfile \n";
  }
  print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

  # plot v-field
  $vec_info = "-N -A0.5p/3p/1.5p -Se${vscale}/${vec_conf}/0 -W0.25p";
  print CSH "psvelo $gps_vec $J1 $R1 $vec_info -K -V -O >> $psfile \n";

  # plot spline centers
  #print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

  # plot color scale and v-field scale
  if ($i==2) {
    print CSH "psscale -C$cptvmag $Dscale $Bscale -K -O -V >> $psfile \n";

    print CSH "psvelo $J_title $R_title $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_conf $vec_conf\nEOF\n";
    print CSH "pstext $J_title $R_title -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 -0.08 10 0 $fontno LM ${vscale_text1} \nEOF\n";
    print CSH "pstext $J_title $R_title -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 -0.16 10 0 $fontno LM ${vscale_text2} \nEOF\n";
  }

  # plot faults and/or plate boundary
  #print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n";

  # plot pole of rotational field
  print CSH "psxy $J1 $R1 $B $poleinfo2 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n";

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno LM $titles[$i-1] \nEOF\n";
}

#---------------------------------------

@idata_vel = (21,23);
@Bins = (3,15);
$Bscale = sprintf("-B%2.2f:\"  \":/:\"10\@+%i\@+ yr\@+-1\@+\": -Ef5p",$cmaxrot*0.25,$cpwr4);
@titles = ("(b)  Rotation rate  (regular grid)","(d)  Rotation rate  (irregular grid)");
@shifts = ($shift2,$shift1);

# plot the rotation-rate field
for ($i = 1; $i <= 2; $i = $i+1) {

  $idata = $idata_vel[$i-1];
  $dlab = sprintf("d%2.2i",$idata);
  $gps_pts = "${gps_dir}/synthetic/syn_vfield_${dlab}_points.dat";
  $gps_vec = "${gps_dir}/synthetic/syn_vfield_${dlab}_psvelo.dat"; 
  if (not -f $gps_pts) {die("Check if $gps_pts exist or not\n")}
  if (not -f $gps_vec) {die("Check if $gps_vec exist or not\n")}
  $name = "${itag}_${dlab}_${qlab}_${blab}_${nlab}_${slab}";

    if ($idata == 20) {$imask = 0;
    } elsif ($idata == 21) {$imask = 0;
    } elsif ($idata == 22) {$imask = 1;
    } elsif ($idata == 23) {$imask = 1;
    }

  $mask_file = "${vel_dir0}/misc/${name}_masked_pts.dat";
  if (not -f ${mask_file}) {die("Check if ${mask_file} exist or not\n")}

  # v-field files
  $strain_mag        = "${vel_dir0}/misc/${name}_spline_strain.dat";
  $spline_centers    = "${vel_dir0}/misc/${name}_spline_gridpoints.dat";
  $euler_poles       = "${vel_dir0}/misc/${name}_euler_vector_psxy.dat";
  $euler_anti_poles  = "${vel_dir0}/misc/${name}_euler_anti_vector_psxy.dat";
  $euler_poles_scale = "${vel_dir0}/misc/${name}_euler_vector_scale_psxy.dat";
  if (not -f $strain_mag) {die("Check if $strain_mag exist or not\n")} 
  if (not -f ${spline_centers}) {die("Check if ${spline_centers} exist or not\n")}
  if (not -f ${euler_poles}) {die("Check if ${euler_poles} exist or not\n")}
  if (not -f ${euler_anti_poles}) {die("Check if ${euler_anti_poles} exist or not\n")}
  #if (not -f ${euler_poles_scale}) {die("Check if ${euler_poles_scale} exist or not\n")}

  $iB = $Bins[$i-1];
  $B = $B0.$Bopts[$iB];

  print CSH "psbasemap $J1 $R1 $B -K -O -V $shifts[$i-1] >> $psfile\n";
  if ($icolor==1) {
    print CSH "awk '{print \$1,\$2,\$6/$norm4}' $strain_mag | surface -G$grdfile ${interp_surf} $R1 \n";
    print CSH "grdimage $grdfile $R1 $J1 -C$cptrot -T -K -O -V >> $psfile\n";
  }
  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V >> $psfile \n";
    print CSH "psmask -C -K -O -V >> $psfile \n";
  }
  print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

  # plot faults and/or plate boundary
  #print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";

  #-----------------------------

  if ($i == 2) {

    # colorbar for scalar field
    print CSH "psscale -C$cptrot $Dscale $Bscale -K -O -V >> $psfile \n";

    # plot euler vector scale
    @vec_ref = ($evec_ref/2, $evec_ref, 1.5*$evec_ref);
    $sv1 = sprintf("%.2f",$vec_ref[0]);
    $sv2 = sprintf("%.2f",$vec_ref[1]);
    $sv3 = sprintf("%.2f",$vec_ref[2]);
    #open(IN,${euler_poles_scale}); $fscale = <IN>; chomp($fscale);
    #if (not -f ${euler_poles_scale}) {die("Check if ${euler_poles_scale} exist or not\n");}

  $dsize[0] = $escale*abs($vec_ref[0]);
  $dsize[1] = $escale*abs($vec_ref[1]);
  $dsize[2] = $escale*abs($vec_ref[2]);

    $dy = 0.25;
    for ($k = 0; $k < 5; $k = $k+1) {
      $ytxt[$k] = $k*$dy;
    }
    $x1 = 0;
    $x2 = $x1 + 0.7;
    print CSH "psxy -JX1 $R_title -C$cptfile -N ${seis_info3} ${origin_dots} -K -O -V >>$psfile<<EOF
$x1 $ytxt[3] $vec_ref[2] $dsize[2]
$x1 $ytxt[2] $vec_ref[1] $dsize[1]
$x1 $ytxt[1] $vec_ref[0] $dsize[0]
$x2 $ytxt[3] -$vec_ref[2] $dsize[2]
$x2 $ytxt[2] -$vec_ref[1] $dsize[1]
$x2 $ytxt[1] -$vec_ref[0] $dsize[0]
EOF\n";

    $x1 = ($x1+$x2)/2;
    $y1 = -0.8*$dy;
  print CSH "pstext -JX1 $R_title -N -K -O -V ${origin_dots} >>$psfile<<EOF
$x1 $ytxt[4] 10 0 $fontno CM pos              neg
$x1 $ytxt[3] 10 0 $fontno CM $sv3
$x1 $ytxt[2] 10 0 $fontno CM $sv2
$x1 $ytxt[1] 10 0 $fontno CM $sv1
$x1 $ytxt[0] 10 0 $fontno CM $unit_rot
EOF\n";
  }

  # plot faults and/or plate boundary
  #print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n";

  # make colorpoint file
  $cptfile = "color.cpt"; $cmin = -0.01; $cmax = 0.01; $dc = ($cmax-$cmin)/2;
  $T = "-T$cmin/$cmax/$dc -D -Z "; $colorbar = "-Cpolar";
  print CSH "makecpt $colorbar $T > $cptfile \n";

  # KEY: plot euler poles
  print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$3/$norm4,$escale*sqrt(\$3*\$3)/$norm4}' ${euler_anti_poles} | psxy $J1 $R1 -C$cptfile ${seis_info3} -K -V -O >> $psfile\n";

  #print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n";  # star

  # plot mean pole
  $epoles = "${vel_dir0}/misc/${name}_epole_mean_points.dat";
  if (not -f $epoles) { die("Check if $epoles exist or not\n") }
  print CSH "psxy $epoles $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile\n";

  print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno LM $titles[$i-1] \nEOF\n";

}

#-----
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk \nEOF\n"; # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if ($ixv==1) {
  print CSH "xv $jpgfile &\n";
}

#==================================================

close (CSH);
system("csh -f $cshfile");
#print "convert $psfile $jpgfile \n";
#system("convert $psfile $jpgfile");
#if($ipdf==1) {system("ps2pdf $psfile")};
#if($ixv==1) {system("xv $jpgfile &")};

#==================================================
