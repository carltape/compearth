gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.0c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p
psbasemap -JH0/3i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G255/255/255 -K -V -P -X2 -Y1 > lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_03.lonlat -W2p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_04.lonlat -W2p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_05.lonlat -W2p,0 -W2p,50/205/50 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beach_arc_06.lonlat -W2p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
0.0000 -90.0000
EOF
pstext -N -J -R -K -O -V -D0.00p/-10.00p >>lune_hammer_iplot1.ps<<EOF
0.0000 -90.0000 14 0 1 CT (-1,-1,-1)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
-30.0000 -54.7356
EOF
pstext -N -J -R -K -O -V -D-10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 -54.7356 14 0 1 RM (0,-1,-1)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
-30.0000 0.0000
EOF
pstext -N -J -R -K -O -V -D-10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 0.0000 14 0 1 RM (2,-1,-1)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
-30.0000 35.2644
EOF
pstext -N -J -R -K -O -V -D-10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 35.2644 14 0 1 RM (1,0,0)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
0.0000 90.0000
EOF
pstext -N -J -R -K -O -V -D0.00p/10.00p >>lune_hammer_iplot1.ps<<EOF
0.0000 90.0000 14 0 1 CB (1,1,1)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
30.0000 54.7356
EOF
pstext -N -J -R -K -O -V -D10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 54.7356 14 0 1 LM (1,1,0)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
30.0000 0.0000
EOF
pstext -N -J -R -K -O -V -D10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 0.0000 14 0 1 LM (1,1,-2)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
30.0000 -35.2644
EOF
pstext -N -J -R -K -O -V -D10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 -35.2644 14 0 1 LM (0,0,-1)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
0.0000 0.0000
EOF
pstext -N -J -R -K -O -V -D10.00p/10.00p >>lune_hammer_iplot1.ps<<EOF
0.0000 0.0000 14 0 1 LB (1,0,-1)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
-30.0000 60.5038
EOF
pstext -N -J -R -K -O -V -D-10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
-30.0000 60.5038 14 0 1 RM (3,1,1)
EOF
psxy -N -Sp8p -W1p,0 -J -R -K -O -V >>lune_hammer_iplot1.ps<<EOF
30.0000 -60.5038
EOF
pstext -N -J -R -K -O -V -D10.00p/0.00p >>lune_hammer_iplot1.ps<<EOF
30.0000 -60.5038 14 0 1 LM (-1,-1,-3)
EOF
psxy ./dfiles//beachpts_Walter2010_points.dat -N -Sc8p -W0.5p,0/0/0 -G160/32/240 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Walter2009_points.dat -N -Sc8p -W0.5p,0/0/0 -G50/205/50 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Ford2009_points.dat -N -Sc8p -W0.5p,0/0/0 -G200 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Minson2007_points.dat -N -Sc8p -W0.5p,0/0/0 -G255/0/0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy ./dfiles//beachpts_Foulger2004_points.dat -N -Sc8p -W0.5p,0/0/0 -G255/165/0 -J -R -K -O -V >>lune_hammer_iplot1.ps
psxy -N -Sc8p -W0.5p,0/0/0 -G160/32/240 -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 1.2
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 1.2 12 0 1 LM Walter2010 (n=14)
EOF
psxy -N -Sc8p -W0.5p,0/0/0 -G50/205/50 -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0.9
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0.9 12 0 1 LM Walter2009 (n=13)
EOF
psxy -N -Sc8p -W0.5p,0/0/0 -G200 -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0.6
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0.6 12 0 1 LM Ford2009 (n=32)
EOF
psxy -N -Sc8p -W0.5p,0/0/0 -G255/0/0 -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0.3
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0.3 12 0 1 LM Minson2007 (n=18)
EOF
psxy -N -Sc8p -W0.5p,0/0/0 -G255/165/0 -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
0 0
EOF
pstext -N -R0/1/0/1 -JX1i -Xa3.0 -Ya7.2 -K -O -V >>lune_hammer_iplot1.ps<<EOF
 0.2 0 12 0 1 LM Foulger2004 (n=26)
EOF
pstext -R0/1/0/1 -JX1i -Xa2 -Ya8.7 -O -V >>lune_hammer_iplot1.ps<<EOF
 -1 -1 12 0 1 LM TEST
EOF
