gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.3c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p
psbasemap -JH0/3i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -V -P -X2 -Y1 > lune_hammer_iplot2_kplot4.ps
psxy ./dfiles//beach_patch_01.lonlat -G120 -J -R -K -O -V >>lune_hammer_iplot2_kplot4.ps
psxy ./dfiles//beach_patch_02.lonlat -G255 -J -R -K -O -V >>lune_hammer_iplot2_kplot4.ps
psbasemap -JH0/3i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot2_kplot4.ps
psxy ./dfiles//beach_arc_01.lonlat -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_kplot4.ps
psxy ./dfiles//beach_arc_02.lonlat -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_kplot4.ps
psxy ./dfiles//beach_arc_03.lonlat -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_kplot4.ps
psxy ./dfiles//beach_arc_04.lonlat -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_kplot4.ps
psxy ./dfiles//beach_arc_05.lonlat -W3p,0/0/0 -J -R -K -O -V >>lune_hammer_iplot2_kplot4.ps
psmeca ./dfiles//beachballs_4_psmeca -JH0/3i -R-30/30/-90/90 -Sm0.5 -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_kplot4.ps
pstext -N -R0/1/0/1 -JX1i -Xa-1 -Ya9.0 -K -O -V >>lune_hammer_iplot2_kplot4.ps<<EOF
 0 0 14 0 1 LM 
EOF
pstext -N -R0/1/0/1 -JX1i -Xa-1 -Ya8.7 -O -V >>lune_hammer_iplot2_kplot4.ps<<EOF
 0 0 11 0 1 LM 
EOF
