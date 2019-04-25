gmtset PAPER_MEDIA Custom_9.5ix21.5i PAGE_ORIENTATION landscape MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.0c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10 HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -V -X0.5 -Y0.5 > lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | psxy -G120 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | psxy -G255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psmeca ../dfiles//beachballs_ipts1_iref1_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sm0.45/8p -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
pstext -N -R0/1/0/1 -JX1i -Xa0.5 -Ya8.6 -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps<<EOF
 0 0 16 0 1 LM Reference sets of moment tensors (input files available in carltape compearth github repository)
EOF
pstext -N -R0/1/0/1 -JX1i -Xa0.5 -Ya8.3 -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps<<EOF
 0 0 12 0 1 LM Plotted using GMT 4.5.3 with a modified version of psmeca
EOF
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -O -V -X3.5 >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | psxy -G120 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | psxy -G255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psmeca ../dfiles//beachballs_ipts1_iref2_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sm0.45/8p -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -O -V -X3.5 >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | psxy -G120 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | psxy -G255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psmeca ../dfiles//beachballs_ipts1_iref3_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sm0.45/8p -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -O -V -X3.5 >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | psxy -G120 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | psxy -G255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psmeca ../dfiles//beachballs_ipts1_iref4_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sm0.45/8p -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -O -V -X3.5 >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | psxy -G120 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | psxy -G255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psmeca ../dfiles//beachballs_ipts1_iref5_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sm0.45/8p -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -O -V -X3.5 >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | psxy -G120 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | psxy -G255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psbasemap -JH0/2.8i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1.ps
psmeca ../dfiles//beachballs_ipts2_iref3_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sm0.45/8p -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1.ps
psxy -R -J -O -T -V >> lune_hammer_iplot2_lplot1_kplot1.ps
