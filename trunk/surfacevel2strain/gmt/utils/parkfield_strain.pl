#!/usr/bin/perl -w

#==========================================================
#
#  parkfield_strain.pl
#  Carl Tape
#  23-Aug-2007
#  
#  This script inputs a surface velocity field and strain scalar fields
#  computed in test_platemodel2strain.m and outputs a figure.
#  
#==========================================================

$cshfile = "parkfield_strain.csh";

$icolor = 1;    # ccc

# plates and faults
$plate_dir = "/home/carltape/gmt/plates";
$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
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
#$vel_dir0   = "${plate_dir}/surface_velocities";
$vel_dir0   = "/home/carltape/SURFACEVEL2STRAIN/matlab_output";

#----------------------------------

# file name, ticks, range for velocities

$iregion = 8;      # region (1=west_us, 2=cal, 3=socal, 6=cascadia, 8=parkfield)
$idata = 1;        # choose GPS dataset (1 = NASA REASON; 2 = CCMM; 10-13 = strike-slip; 20-23 = rotational)
$ndim = 3;         # ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 6;         # min grid order used
$qmax = 9;         # max grid order used
$basis = 1;        # 1 for wavelets; 2 for splines
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$iunrotate = 1;    # 1 to remove rotation; 0 to leave original field

$imask  = 1;       # plot mask
$ieuler = 0;       # plot euler poles
$ifault = 0;       # plot faults
$iplate = 0;       # plot plate boundaries

   $itag = "parkfield";
   $z_title = 0.8;
   $xtick1 = 0.2; $xtick2 = 0.1; $wid = 6.0;
   $origin = "-X1.0 -Y2.5";
   $xlegend = 0.25; $ylegend = -1.0;

$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";
$evec_ref = 1.0; $escale = 8;    # default values
$fac_min = 1.0; $fac_max = 1.0;  # default values

# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
$dlab = sprintf("d%2.2i",$idata);

  $tag = "NASA REASON data set (cGPS)";
  #$gps_pts = "${gps_dir}/US/reason_fixed_NAM_subset_points.dat";
  #$gps_vec = "${gps_dir}/US/reason_fixed_NAM_subset_psvelo.dat";
  $gps_pts = "${gps_dir}/US/reason_subset_points.dat";
  $gps_vec = "${gps_dir}/US/reason_subset_psvelo.dat";
  $igc = 0; $vec_conf = 0.95; $vec_error = 0.5;
  $cmin = 5; $cmax = 50; $ctick = 10; $vscale = 0.04; $vec_ref = 10;
  $cmax = 45; $ctick = 10;
  $vumin = -2; $vumax = 2; $vsmin = -35; $vsmax = 5; $vemin = -35; $vemax = 5;
  $ifault = 1; $iplate = 0;
  $fac_min = 0.5;    # 1.0 for uniform color scaling (try 0.5)
  $fac_max = 1.0;    # 1.0 default

if (not -f $gps_pts) { die("Check if $gps_pts exist or not\n") }
if (not -f $gps_vec) { die("Check if $gps_vec exist or not\n") }

#$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);
$vscale_text = sprintf("%.0f mm/yr", $vec_ref);

#--------------------------------------------

# plotting specifications
$fsize0 = "20";
$fsize1 = "12";
$fsize2 = "12";
$fsize3 = "6";
$fontno = "1";    # 1 or 4
$tick   = "0.3c";
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

$plate_infoG     = "-M -W2p,0/255/0";
$plate_infoR     = "-M -W2p,255/0/0";
$plate_infoK     = "-M -W2p,0/0/0";
$plate_infoW     = "-M -W2p,255/255/255";
$plate_infoGr    = "-M -W3p,200";

$textinfo = "-G255 -S1p";

# which borders to plot the lat-lon
# 1 four sides, 4 single sides, 6 two sides, 4 three sides, 1 zero sides
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");

#-----------------------------------------------
# BOOLEAN: WHICH FIGURES TO PLOT
$ixv    = 1;
$ipdf   = 0;

$ivel = 0;
$istrain = 1;
$irotate = 0;
$idilate = 0;
$imaxlam = 0;
$imap = 0;

# 2-column figures
$ivel3D = 0;         # estimated velocity field (Vmag, Vr)
#$ivelstr = 0;
#$imaskfig = 0;
#$isetup = 0;

# 3-column figures
$imulti_vel = 0;
$imulti_strain_inc = 0;
$imulti_strain_cum = 0;
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
$latmin = 35.5; $latmax = 36.1; $lonmin = -121.2; $lonmax = -119.8;
$R1 = "-R$lonmin/$lonmax/$latmin/$latmax";
print "\n bounds are $R1\n\n";
#die("testing");

# get color limits for the plots: dilatation, strain, rotation, maxlam
open(IN,"$colors"); @lines = <IN>;
($cmindilat,$cmaxdilat,$cpwr2) = split(" ",$lines[0]);  # dilatation
($cminstrain,$cmaxstrain,$cpwr3) = split(" ",$lines[1]);  # strain
$cmaxstrain = 16;
($cminrot,$cmaxrot,$cpwr4) = split(" ",$lines[2]);  # rotation
($cmin5,$cmax5,$cpwr5) = split(" ",$lines[3]);  # maxlam
#$norm2 = "1e$cpwr2"; $norm3 = "1e$cpwr3"; $norm4 = "1e$cpwr4"; $norm5 = "1e$cpwr5";
$norm2 = 10**$cpwr2; $norm3 = 10**$cpwr3; $norm4 = 10**$cpwr4; $norm5 = 10**$cpwr5;
print "\n $cpwr2 $cpwr3 $cpwr4 $cpwr5 $norm2 $norm3 $norm4 $norm5\n\n";

$title_strain  = sprintf("Estimated strain, 10\@+%i\@+ yr\@+-1\@+",$cpwr3);   # %2.2i for 07
$title_rot    = sprintf("Estimated rotation, 10\@+%i\@+ yr\@+-1\@+",$cpwr4);
$title_dilat  = sprintf("Estimated dilatation, 10\@+%i\@+ yr\@+-1\@+",$cpwr2);

$title_up    = "Estimated Vr (mm/yr)";
$title_south = "Estimated V$theta (mm/yr)";
$title_east  = "Estimated V$phi (mm/yr)";

$Bscale_strain  = sprintf("-B%2.2f:\"Strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);
$Bscale_rot    = sprintf("-B%2.2f:\"Rotation rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxrot*0.25,$cpwr4);
$Bscale_dilat  = sprintf("-B%2.2f:\"Dilatation rate, 10\@+%i\@+ yr\@+-1\@+\": -E10p",$cmaxdilat*0.5,$cpwr2);
$unit_rot      = sprintf("10\@+%i\@+ yr\@+-1\@+",$cpwr4);

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

$colorbar = "seis -I";

# v-field components
$cptup = "color_up.cpt"; print COLOR "makecpt -C$colorbar $Tup -D > $cptup\n";
$cptsouth = "color_south.cpt"; print COLOR "makecpt -C$colorbar $Tsouth -D > $cptsouth\n";
$cpteast = "color_east.cpt"; print COLOR "makecpt -C$colorbar $Teast -D > $cpteast\n";

$colorbar = "rainbow";

$cptvmag = "color1.cpt";
print COLOR "makecpt -C$colorbar $Tvmag -D > $cptvmag\n";
#print COLOR "sed 's/^B.*/B       255   255    255  /' temp1  >  temp2\n";
#print COLOR "sed 's/^F.*/F       255     0      0  /' temp2 > $cptvmag\n";

# if(1==1) {
#   $colorbar = $c_romania;
#   $colorbar = "seis -I";
#   #$cback = "  0     0    255";
#   #$cfor  = "255     0      0";

# } else {
#   $colorbar = "no_green";
#   #$cback = "170     0      0";
#   #$cfor  = "  0     0    200";
# }

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

$B = $B0.$Bopts[12];

# KEY: commands for interpolation and masking
$interp = "-I0.1"; $interp_surf = "-I0.04 -T0 -S0";
#if($iregion==1) {$interp_surf = "-I0.2 -T0 -S0";}
if($iregion==8) {$interp = "-I0.025"; $interp_surf = "-I0.01 -T0 -S0";}

$mask_info = "-I0.1 -G200 -S0.2";
if($iregion==8) {$mask_info = "-I0.01 -G200 -S0.02";}

$grdfile = "temp.grd";
$mask_file = "${vel_dir0}/${name}_masked_pts.dat";
if (not -f ${mask_file}) { die("Check if ${mask_file} exist or not\n") }

# open CSH file
open(CSH,">$cshfile");
print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter MEASURE_UNIT inch PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen D_FORMAT \%.1f\n";

$pname = $name;

#-----------------------------------------

if ($istrain==1){

$fname = "parkfield_strain_${qlab}_${nlab}";
$psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

$title = "Strain rate from $tag";
$title = " ";
$Bscale  = sprintf("-B%2.2f:\"Strain rate, 10\@+%i\@+ yr\@+-1\@+\": -Ef10p",$cmaxstrain*0.25,$cpwr3);

# plot basemap
print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n";  # START
#print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | pscontour $R1 $J1 -A- -C$cptstrain -I -P -O -K -V >> $psfile\n";

# make grd file, then plot
if (0==1) {
  #$interp = "-I0.1";
  print CSH "awk '{print \$1,\$2,\$5}' $strain_mag | xyz2grd -G$grdfile $interp $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -T -K -O -V -P >> $psfile\n";

} else {

  #print CSH "awk '{print \$1,\$2,\$5}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "awk '{print \$1,\$2,\$5/$norm3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1\n";
  print CSH "grdimage $grdfile $R1 $J1 -C$cptstrain -T -K -O -V -P >> $psfile\n";

  if ($imask==1) {
    print CSH "echo mask file is ${mask_file}\n";
    print CSH "psmask ${mask_file} $R1 $J1 $mask_info -K -O -V -P >> $psfile\n";
    print CSH "psmask -C -K -O -V -P >> $psfile\n";
  }

}

# KEY COMMAND: create data file
#$gpstemp = $gps_vec;
$gpstemp = gpstemp1;
print CSH "awk '{ print \$1,\$2,\$3,\$4,\$5,\$6}' ${velocity_vec_dat} > $gpstemp\n";

# SAFOD hole
$plat = 35.974;
$plon = -120.552;
$poleinfo2 = "-Sc15p -G255 -W0.75p,0/0/0";
#print CSH "psxy $J1 $R1 $poleinfo2 -K -O -V >>$psfile<<EOF\n$plon $plat\nEOF\n";

# plot v-field
$vec_info = "-N -A1p/6p/3p -Se${vscale}/${vec_conf}/0 -W0.5p";
print CSH "psvelo $gpstemp $J1 $R1 $vec_info -K -V -O >> $psfile\n";

# plot spline centers
#print CSH "psxy ${spline_centers} $J1 $R1 $sinfo -K -O -V >> $psfile\n";

# plot v-field scale
print CSH "pstext -JX1/1 -R0/1/0/1 -N -K -O -V $origin_arrow >>$psfile<<EOF\n 0 0.25 13 0 $fontno LM ${vscale_text}\nEOF\n";
print CSH "psvelo -JX1/1 -R0/1/0/1 $vec_info $origin_arrow -K -O -V >>$psfile<<EOF\n 0 0 $vec_ref 0.00 $vec_error $vec_error\nEOF\n";

# color bar
print CSH "gmtset D_FORMAT \%.0f\n";
print CSH "psscale -C$cptstrain $Dscale $Bscale -K -O -V -P >> $psfile\n";
print CSH "gmtset D_FORMAT \%.1f\n";

#print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V -P >> $psfile\n";

# plot faults and/or plate boundary
if ($ifault==1) {
  #print CSH "psxy $kcf_file $J1 $R1 $fault_info_k -K -V -O -P >> $psfile\n";
  print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O -P >> $psfile\n";
}
if ($iplate==1) {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoW -K -O -V >> $psfile\n";}
if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoW -K -O -V >> $psfile\n";}

print CSH "psbasemap $J1 $R1 $B -K -V -O -P >> $psfile\n";

# plot pole of rotational field
#if($ieuler == 1) {print CSH "psxy $J1 $R1 $B $poleinfo1 -K -O -V >>$psfile<<EOF\n -116.0 35.0\nEOF\n"}

print CSH "pstext -N $J_title $R_title -K -O -V -P >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $title\nEOF\n";

#--------------------------------------
print CSH "pstext -N -JX10 -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 10 0 1 CM junk\nEOF\n";  # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==================================================

close (CSH);
system("csh -f $cshfile");

#==================================================
