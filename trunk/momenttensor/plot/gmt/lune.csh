gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.3c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p
psbasemap -JH0/3i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -V -P -X2 -Y1 > lune_hammer_iplot1.ps
psxy ./dfiles//beach_patch_01.lonlat -G120 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_patch_02.lonlat -G255 -J -R -K -O -V >>lune_hammer_iplot1.ps
psbasemap -JH0/3i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_01.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_02.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_03.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_04.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_05.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_06.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps
pstext -N -J -R -K -O -V -D0.00p/-15.00p >>lune_hammer_iplot1.ps<<EOF
0.0000 -90.0000 14 0 1 CM ISO
EOF
pstext -N -J -R -K -O -V -D15.00p/15.00p >>lune_hammer_iplot1.ps<<EOF
0.0000 0.0000 14 0 1 CM DC
EOF
pstext -N -J -R -K -O -V -D0.00p/15.00p >>lune_hammer_iplot1.ps<<EOF
0.0000 90.0000 14 0 1 CM ISO
EOF
pstext -N -J -R -K -O -V -D-18.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 60.5038 14 0 1 CM C
EOF
pstext -N -J -R -K -O -V -D-22.50p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 35.2644 14 0 1 CM LVD
EOF
pstext -N -J -R -K -O -V -D-30.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 0.0000 14 0 1 CM CLVD
EOF
pstext -N -J -R -K -O -V -D-15.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 -54.7356 14 0 1 CM -
EOF
pstext -N -J -R -K -O -V -D18.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 -60.5038 14 0 1 CM C
EOF
pstext -N -J -R -K -O -V -D22.50p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 -35.2644 14 0 1 CM LVD
EOF
pstext -N -J -R -K -O -V -D30.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 0.0000 14 0 1 CM CLVD
EOF
pstext -N -J -R -K -O -V -D15.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 54.7356 14 0 1 CM -
EOF
psxy ./dfiles//beachpts_Ford2009_points.dat -N -Sc8p -W0.5p,0/0/0 -G255/0/0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Foulger2004_points.dat -N -Sc8p -W0.5p,0/0/0 -G255/165/0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Minson2007_points.dat -N -Sc8p -W0.5p,0/0/0 -G50/205/50 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Minson2008_points.dat -N -Sc8p -W0.5p,0/0/0 -G255/255/255 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Walter2009_points.dat -N -Sc8p -W0.5p,0/0/0 -G0/255/255 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Walter2010_points.dat -N -Sc8p -W0.5p,0/0/0 -G160/32/240 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Pesicek2012_points.dat -N -Sc8p -W0.5p,0/0/0 -G0/0/0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy -N -Sc8p -W1p,0/0/0 -G255/0/0 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 1.2
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 1.2 12 0 1 LM Ford 2009 (n=32)
EOF
psxy -N -Sc8p -W1p,0/0/0 -G255/165/0 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0.9
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0.9 12 0 1 LM Foulger 2004 (n=26)
EOF
psxy -N -Sc8p -W1p,0/0/0 -G50/205/50 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0.6
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0.6 12 0 1 LM Minson 2007 (n=18)
EOF
psxy -N -Sc8p -W1p,0/0/0 -G255/255/255 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0.3
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0.3 12 0 1 LM Minson 2008 (n=7)
EOF
psxy -N -Sc8p -W1p,0/0/0 -G0/255/255 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0 12 0 1 LM Walter 2009 (n=13)
EOF
psxy -N -Sc8p -W1p,0/0/0 -G160/32/240 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 -0.3
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 -0.3 12 0 1 LM Walter 2010 (n=14)
EOF
psxy -N -Sc8p -W1p,0/0/0 -G0/0/0 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 -0.6
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 -0.6 12 0 1 LM Pesicek 2012 (n=7)
EOF
pstext -R0/1/0/1 -JX1i -Xa-1 -Ya8.7 -O -V >>lune_hammer_iplot1.ps<<EOF
 -1 -1 11 0 1 LM TEST
EOF
