#!/usr/bin/perl -w

#==========================================================
#
#  mogi_dilat.pl
#  Carl Tape
#  05-Aug-2008
#  
#  This script inputs a surface velocity field and strain scalar fields
#  computed in test_platemodel2strain.m and outputs a figure.
#  
#==========================================================

$cshfile = "mogi_dilat.csh";

$icolor = 1;    # ccc

# plates and faults
$plate_dir = "/home/carltape/gmt/plates";
$plate_file  = "${plate_dir}/plate_boundaries/bird_boundaries";
if (not -f ${plate_file}) { die("Check if plate_file ${plate_file} exist or not\n") }
$kcf_file     = "/home/carltape/gmt/faults/kcf.xy";
$fault_file   = "/home/carltape/gmt/faults/jennings.xy";
if (not -f $fault_file) { die("Check if fault_file $fault_file exist or not\n") }

$fault_info_k = "-M -W1.5p,0/0/0";
$fault_info_r = "-M -W1.5p,255/0/0";

$sinfo = "-Sc4p -G255 -W0.5p,0/0/0";
$poleinfo1 = "-Sa20p -G0/255/255 -W0.75p,0/0/0";
$circleinfo = "-Sc30p -W1.0p,0/0/0,--";

# velocity, strain rate, etc
#$vel_dir0   = "${plate_dir}/surface_velocities";
$vel_dir0   = "/home/carltape/SURFACEVEL2STRAIN/matlab_output";

#----------------------------------

# file name, ticks, range for velocities

$iregion = 3;      # region (1=west_us, 2=cal, 3=socal, 6=cascadia, 8=parkfield)
$idata = 83;        # choose GPS dataset (1 = NASA REASON; 2 = CCMM; 10-13 = strike-slip; 20-23 = rotational)
$ndim = 3;         # ndim = 3 for verticals; ndim = 2 for horizontals
$qmin = 3;         # min grid order used
$qmax = 8;         # max grid order used
$basis = 1;        # 1 for wavelets; 2 for splines
$iLmat = 1;        # sopt : 0 (full L), 1 (elasticity), 2 (viscosity)
$iunrotate = 0;    # 1 to remove rotation; 0 to leave original field

$imask  = 1;       # plot mask
$ieuler = 0;       # plot euler poles
$ifault = 0;       # plot faults
$iplate = 0;       # plot plate boundaries

   $itag = "socal";
   #$z_title = 1.2;
   $z_title = 1.05;
   $xtick1 = 1; $xtick2 = 0.25; $wid = 6.0;
   #$origin = "-X1.0 -Y1.0";
   #$xlegend = 0.25; $ylegend = 1.3;
   $origin = "-X1.0 -Y2.5";
   $xlegend = 0.25; $ylegend = -1.0;

$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";
$evec_ref = 1.0; $escale = 8;    # default values
$fac_min = 1.0; $fac_max = 1.0;  # default values

# VELOCITY FIELD DATASET
$gps_dir  = "/home/carltape/gmt/gps_data";
$dlab = sprintf("d%2.2i",$idata);

# SYNTHETIC velocity field
  # volcanic dilitation field
 $gps_pts = "${gps_dir}/synthetic/syn_vfield_${dlab}_points.dat";
  $gps_vec = "${gps_dir}/synthetic/syn_vfield_${dlab}_psvelo.dat"; 

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

if (not -f $gps_pts) { die("Check if gps_pts $gps_pts exist or not\n") }
if (not -f $gps_vec) { die("Check if gps_vec $gps_vec exist or not\n") }

$vscale_text = sprintf("%.0f mm/yr  (%.0f%%)", $vec_ref, $vec_conf*100);

#$xtick1 = 1; $xtick2 = 0.25; $wid = 5.5; $cmin = 0; $cmax = 40; $ctick = 10;
#$xtick1 = 1; $xtick2 = 0.25; $wid = 3.25; $cmin = 0; $cmax = 40; $ctick = 10;
#$xtick1 = 1; $xtick2 = 0.25; $wid = 5.5; $cmin = 44; $cmax = 52; $ctick = 2;   # plate field (socal)

#--------------------------------------------

# plotting specifications
$fsize0 = "20";
$fsize1 = "14";
$fsize2 = "12";
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

# 3-column figures
$imulti_vel = 0;
$imulti_strain_inc = 1;
$imulti_strain_cum = 0;

if($imulti_vel == 1 && $ndim == 2) {die("\n imulti_vel = 1, so ndim must be 3\n")}

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
print "\n bounds are $R1 \n\n";
#die("testing");

# get color limits for the plots: dilatation, strain, rotation, maxlam
open(IN,"$colors"); @lines = <IN>;
($cmindilat,$cmaxdilat,$cpwr2) = split(" ",$lines[0]);  # dilatation
$cmaxdilat = 3;
$cmindilat = -$cmaxdilat;

($cminstrain,$cmaxstrain,$cpwr3) = split(" ",$lines[1]);  # strain
($cminrot,$cmaxrot,$cpwr4) = split(" ",$lines[2]);  # rotation
($cmin5,$cmax5,$cpwr5) = split(" ",$lines[3]);  # maxlam
#$norm2 = "1e$cpwr2"; $norm3 = "1e$cpwr3"; $norm4 = "1e$cpwr4"; $norm5 = "1e$cpwr5";
$norm2 = 10**$cpwr2; $norm3 = 10**$cpwr3; $norm4 = 10**$cpwr4; $norm5 = 10**$cpwr5;
print "\n $cpwr2 $cpwr3 $cpwr4 $cpwr5 $norm2 $norm3 $norm4 $norm5 \n\n";

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

print "\n @vtick \n";

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

#$colorbar = "rainbow";
$colorbar = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";

# v-field components
$cptup = "color_up.cpt"; print COLOR "makecpt -C$colorbar $Tup -D > $cptup \n";
$cptsouth = "color_south.cpt"; print COLOR "makecpt -C$colorbar $Tsouth -D > $cptsouth \n";
$cpteast = "color_east.cpt"; print COLOR "makecpt -C$colorbar $Teast -D > $cpteast \n";

$cptvmag = "color1.cpt";
print COLOR "makecpt -C$colorbar $Tvmag -D > $cptvmag \n";
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
$interp = "-I0.1";
$interp_surf = "-I0.04 -T0 -S0";
#if($iregion==1) {$interp_surf = "-I0.2 -T0 -S0";}

$mask_info = "-I0.1 -G200 -S0.2";
if($iregion==8) {$mask_info = "-I0.01 -G200 -S0.02";}

$grdfile = "temp.grd";
$mask_file = "${vel_dir0}/${name}_masked_pts.dat";
if (not -f ${mask_file}) { die("Check if ${mask_file} exist or not\n") }

# open CSH file
open(CSH,">$cshfile");

$pname = $name;

#==============================
# from here on out, use portrait plotting with 3-column figures

@labs = ("A","B","C","D");
$fault_info_k  = "-M -W0.5p,0/0/0";
$fault_info_r  = "-M -W0.5p,255/0/0";
$coast_infoK   = "$coast_res -W1.0p,0/0/0 -Na/1.0p";
$coast_infoW   = "$coast_res -W1.0p,255/255/255 -Na/1.0p,255/255/255,t";
$plate_infoK   = "-M -W1p,0/0/0";
$plate_infoW   = "-M -W1p,255/255/255";
$plate_infoR   = "-M -W1p,255/0/0";

  $wid = 2.0; $vscale = 0.005;
  $origin = "-X0.75 -Y8.75";
  $xfac = 1.20; $yfac = 1.0;
  $dX1 = $xfac*$wid; $dY1 = 0; 
  $xlab = 0.95*$wid; $ylab = 1.0*$wid;
  $xlegend = 0.25; $ylegend = -1.0;
  $xtick1 = 1; $xtick2 = 1;
  $x_title = 0.5; $z_title = 0.87;
  $x_stitle = 0.05; $z_stitle = -0.08; $fsize_stitle = 8;   # subtitle
  $Dscale = "-D1.4/-0.1/1.0/0.08h";
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter MEASURE_UNIT inch PLOT_DEGREE_FORMAT D TICK_LENGTH 4p LABEL_FONT_SIZE 8 ANOT_FONT_SIZE 9 HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE 8 FRAME_PEN 1p TICK_PEN 1p\n";

$origin_arrow = "-Xa${xlegend} -Ya${ylegend}";

  $latmin = 32; $latmax = 37; $lonmin = -122; $lonmax = -114;
  $R1 = "-R$lonmin/$lonmax/$latmin/$latmax";

$fsize_title = 8;
$J1 = "-JM${wid}i";
$J_title = "-JX${wid}";

$dX2 = -2*$dX1; $dY2 = -$yfac*$wid;
$shift1 = "-X$dX1 -Y$dY1";
$shift2 = "-X$dX2 -Y$dY2";

$B0 = "-Ba${xtick1}f${xtick2}d::";
$B = $B0.$Bopts[5];

@flab = ("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p");

#==============================

if ($imulti_vel == 1 || ($imulti_strain_inc == 1 || $imulti_strain_cum == 1) ) {

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
  open(IN,"$colors"); @clines = <IN>; print "\n @clines \n";

  # get color limits for the plots
  open(IN,"$cticks"); @cticks = <IN>; print "\n @cticks \n";

  # get q indexes for the plots
  open(IN,"$iqs"); @qlines1 = <IN>; print "\n @qlines1";
  open(IN,"$iqs2"); @qlines2 = <IN>; print "\n @qlines2";
  $nump = @qlines1;
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
  print "\n -- nump = $nump -- \n @qtags1a \n @qtags1b \n @qtags2a \n @qtags2b \n";;
  #die("testing");

  #-----------------------------------------
  # make color files

  $colorbar = "/home/carltape/gmt/color_maps/Romanian_flag_smooth.cpt";
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
  $B = $B0.$Bopts[15];

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

    if ($i==0) {
      print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
    } else {
      print CSH "psbasemap $J1 $R1 $B -K -O -V $shift2 >> $psfile\n";
    }
    #print CSH "awk '{print \$1,\$2,\$3}' $strain_mag | surface -G$grdfile ${interp_surf} $R1 \n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}}' $vu_file | surface -G$grdfile ${interp_surf} $R1 \n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -T -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile \n";
      print CSH "psmask -C -K -O -V >> $psfile \n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";

    if ($ifault == 1 && $iplate == 1) {
       if ($i==0 || $i==$nump-1) {print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
       } else {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";}
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile \n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {
      print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[0] \nEOF\n";
    }

    # subtitle
    $stitle = "($flab[3*$i])  $qtags1b[$i]";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle \nEOF\n";

    #----------------------------------------
    $cptfile = "color_south${i}.cpt";
    $Bscale  = "-B$ctick2:\" \": -E5p";

    print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}}' $vs_file | surface -G$grdfile ${interp_surf} $R1 \n";
      #print CSH "surface ${vs_file} -G$grdfile ${interp_surf} $R1 \n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -T -K -O -V >> $psfile\n";
    }

    #print "\n $cptfile \n $vs_file \n column ${col} \n"; die("testing");

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile \n";
      print CSH "psmask -C -K -O -V >> $psfile \n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
       if ($i==0 || $i==$nump-1) {print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
       } else {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";}
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile \n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {
      print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[1] \nEOF\n";
    }

    # subtitle
    $stitle = "($flab[3*$i+1])  $qtags1b[$i]";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle \nEOF\n";

    #----------------------------------------

    $cptfile = "color_east${i}.cpt";
    $Bscale  = "-B$ctick3:\" \": -E5p";

    print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}}' $ve_file | surface -G$grdfile ${interp_surf} $R1 \n";
      #print CSH "surface ${ve_file} -G$grdfile ${interp_surf} $R1 \n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -T -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile \n";
      print CSH "psmask -C -K -O -V >> $psfile \n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
       if ($i==0 || $i==$nump-1) {print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
       } else {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";}
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile \n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    if ($i == 0) {print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title $fsize_title 0 $fontno CM $titles[2] \nEOF\n";}

    # subtitle
    $stitle = "($flab[3*$i+2])  $qtags1b[$i]";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle $fsize_stitle 0 $fontno LM $stitle \nEOF\n";

  }    # end for loop

print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==============

if ($imulti_strain_inc == 1 || $imulti_strain_cum == 1) {

  if ($imulti_strain_inc == 1) {$multi_strain = $multi_strain_inc; $xtag = "inc";}
  if ($imulti_strain_cum == 1) {$multi_strain = $multi_strain_cum; $xtag = "cum";}

  $fname = "mogi_dilat";
  $psfile = "$fname.ps"; $jpgfile = "$fname.jpg";

  $title1    = sprintf("(a)  Dilatation rate, %s",$qtags1b[1]);
  $subtitle1 = sprintf("10\@+%i\@+ yr\@+-1\@+",$cpwr2);
  $title2    = sprintf("(b)  Dilatation rate, %s",$qtags1b[2]);
  $subtitle2 = sprintf("10\@+%i\@+ yr\@+-1\@+",$cpwr2);
  @titles = ($title1,$title2);
  @stitles = ($subtitle1,$subtitle2);

  $B = $B0.$Bopts[15];

  # KEY: adjust limits for multiscale strain plots
  #$fac_min = 0.50;   # 1.0 for uniform color scaling
  #$fac_max = 1.0;    # 1.0 default
  $fac_inc = ($fac_max - $fac_min)/($nump-2);
  print "\n Color scaling for multiscale strain: $fac_min, $fac_inc, $fac_max\n";

  for ($i = 0; $i < $nump-1; $i = $i+1) {

    $p = $i + 1;
    $col = $p + 2;

    # flag if the mask file is empty
    if ($imask == 1) {
      $mask_file = "${vel_dir0}/${name}_masked_pts_${p}.dat";
      if (not -f $mask_file) {die("check if $mask_file exists");}
      #open(IN,"$mask_file"); @test = <IN>;
      ($nline,$nnum,$nchar,$junk) = split(" ",`wc $mask_file`);
    } else {
      $nline = 0;
    }

    #print "\n $norm $fac\n";

    if ($i==0) {
      $B = $B0.$Bopts[9];
      print CSH "psbasemap $J1 $R1 $B -K -V -P $origin > $psfile\n"; # START
      $cmaxdilat2 = 10;
    } else {
      $B = $B0.$Bopts[3];
      print CSH "psbasemap $J1 $R1 $B -K -O -V $shift1 >> $psfile\n";
      $cmaxdilat2 = 4;
    }

    $fac = $fac_min + ($p-1)*$fac_inc;
    #$cmaxdilat2 = $cmaxdilat*$fac;
    $cmindilat2 = -$cmaxdilat2;
    $dc = ($cmaxdilat2 - $cmindilat2)/100; $Tdilat2 = "-T$cmindilat2/$cmaxdilat2/$dc";
    $cptdilat2 = "cpt_dilat.cpt"; print CSH "makecpt -C$colorbar $Tdilat2 -D > $cptdilat2\n";
    $Bscale_dilat2  = sprintf("-B%2.2f:\" \": -E5p",$cmaxdilat2);

    #----------------------------------------
    # dilatation
    $cptfile = $cptdilat2; $Bscale = $Bscale_dilat2; $col = 2 + 3*$p; $norm = $norm2*$fac;

    #print CSH "psbasemap $J1 $R1 $B $shift1 -K -O -V >> $psfile\n";
    if ($icolor==1) {
      print CSH "awk '{print \$1,\$2,\$${col}/$norm}' $multi_strain | surface -G$grdfile ${interp_surf} $R1 \n";
      print CSH "grdimage $grdfile $R1 $J1 -C$cptfile -T -K -O -V >> $psfile\n";
    }

    if ($nline > 0) {
      print CSH "echo mask file is $mask_file\n";
      print CSH "psmask $mask_file $R1 $J1 $mask_info -K -O -V >> $psfile \n";
      print CSH "psmask -C -K -O -V >> $psfile \n";
    }

    print CSH "pscoast $J1 $R1 $B $coast_infoW -K -O -V >> $psfile\n";
    if ($ifault == 1 && $iplate == 1) {
       if ($i==$nump-2) {print CSH "psxy $fault_file $J1 $R1 $fault_info_k -K -V -O >> $psfile \n";
       } else {print CSH "psxy $J1 $R1 ${plate_file} $plate_infoK -K -O -V >> $psfile\n";}
    }
    if ($igc==1) {print CSH "psxy $J1 $R1 ${gc_boundary} $plate_infoK -K -O -V >> $psfile\n";}
    print CSH "psscale -C$cptfile $Dscale $Bscale -K -O -V >> $psfile \n";

    # highlight the small volcanic source
    if($idata >= 80 and $idata < 90) {print CSH "psxy $J1 $R1 $B $circleinfo -K -O -V >>$psfile<<EOF\n -118.0 34.0\nEOF\n"}

    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_title $z_title 10 0 $fontno CM $titles[$i] \nEOF\n";
    print CSH "pstext -N $J_title $R_title -K -O -V >>$psfile<<EOF\n $x_stitle $z_stitle 10 0 $fontno LM $stitles[$i] \nEOF\n";

  }    # end for loop

print CSH "pstext -N $J_title $R_title -O -V >>$psfile<<EOF\n 10 10 $fsize_title 0 $fontno CM $title \nEOF\n"; # FINISH

print CSH "echo done with $psfile\n";
print CSH "convert $psfile $jpgfile\n";
if($ixv==1) {print CSH "gv $psfile &\n"}

}

#==================================================

close (CSH);
system("csh -f $cshfile");

#==================================================
