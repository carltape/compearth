gmt gmtset PS_MEDIA Custom_9.5ix11i PS_PAGE_ORIENTATION landscape PROJ_LENGTH_UNIT inch MAP_FRAME_TYPE plain MAP_TICK_LENGTH 0.0c MAP_FRAME_PEN 2p MAP_TICK_PEN 2p FONT_ANNOT_PRIMARY 12p,Helevetica,black FONT_HEADING 18p,Helvetiva,black FONT_LABEL 10p,Helvetiva,black
gmt psbasemap -JH0/2.8i -R-30/30/-90/90 -Bxa10f5g10 -Bya10f5g10 -Bwesn+g200 -K -V -X0.5 -Y0.5 > lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | gmt psxy -G120 -L+x0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | gmt psxy -G255 -L+x0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt psbasemap -JH0/2.8i -R-30/30/-90/90 -Bxa10f5g10 -Bya10f5g10 -Bwesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | gmt psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | gmt psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt psmeca ../dfiles//beachballs_ipts1_iref1_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sm0.45/8p -L0.5p,0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt pstext -N -R0/1/0/1 -JX1i -Xa0.5 -Ya8.6 -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps<<EOF
 0 0 16 0 1 LM Examples using the psmeca flags -Sm (full), -Sz (deviatoric), and -Sd (double couple)
EOF
gmt pstext -N -R0/1/0/1 -JX1i -Xa0.5 -Ya8.3 -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps<<EOF
 0 0 12 0 1 LM Plotted using gmt600pr2067_psmeca_flags
EOF
gmt psbasemap -JH0/2.8i -R-30/30/-90/90 -Bxa10f5g10 -Bya10f5g10 -Bwesn+g200 -K -O -V -X3.5 >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | gmt psxy -G120 -L+x0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | gmt psxy -G255 -L+x0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt psbasemap -JH0/2.8i -R-30/30/-90/90 -Bxa10f5g10 -Bya10f5g10 -Bwesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | gmt psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | gmt psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt psmeca ../dfiles//beachballs_ipts1_iref1_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sz0.45/8p -L0.5p,0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt psbasemap -JH0/2.8i -R-30/30/-90/90 -Bxa10f5g10 -Bya10f5g10 -Bwesn+g200 -K -O -V -X3.5 >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_01.dat | gmt psxy -G120 -L+x0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_patch_02.dat | gmt psxy -G255 -L+x0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt psbasemap -JH0/2.8i -R-30/30/-90/90 -Bxa10f5g10 -Bya10f5g10 -Bwesn -K -O -V >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_01.dat | gmt psxy -W3p,160/32/240 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_02.dat | gmt psxy -W3p,255/0/0 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_03.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_04.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
awk '{print $1,$2}' ../dfiles//sourcetype_arc_06.dat | gmt psxy -W3p,30/144/255 -J -R -K -O -V >>lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
gmt psmeca ../dfiles//beachballs_ipts1_iref1_lune_psmeca -JH0/2.8i -R-30/30/-90/90 -Sd0.45/8p -L0.5p,0/0/0 -G255/0/0 -N -K -O -V >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
psxy -R -J -O -T -V >> lune_hammer_iplot2_lplot1_kplot1_gmt600pr2067_psmeca_flags.ps
