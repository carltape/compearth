gmtset PAPER_MEDIA Custom_7.7ix3.2i PAGE_ORIENTATION landscape MEASURE_UNIT inch BASEMAP_TYPE plain TICK_LENGTH 0.0c FRAME_PEN 2p TICK_PEN 2p
psbasemap -JXi/6.72i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -V -X0.5 -Y0.5 > lune_ngeo_iplot2_lplot1_kplot1_min.ps
psmeca ./dfiles//beachballs_ipts1_iref1_lune_psmeca -JXi/6.72i -R-30/30/-90/90 -Sm0.45/8p -L0.5p/0/0/0 -G255/0/0 -N -K -O -V >> lune_ngeo_iplot2_lplot1_kplot1_min.ps
psxy -R -J -O -T -V >> lune_ngeo_iplot2_lplot1_kplot1_min.ps
