gmt gmtset PS_MEDIA Custom_3.2ix7.7i PS_PAGE_ORIENTATION portrait PROJ_LENGTH_UNIT inch MAP_FRAME_TYPE plain MAP_TICK_LENGTH 0.0c MAP_FRAME_PEN 2p MAP_TICK_PEN 2p
gmt psbasemap -JXi/6.72i -R-30/30/-90/90 -Ba10f5g10:" ":/a10f5g10:" ":wesn -G200 -K -V -X0.5 -Y0.5 > lune_gmt600.ps
gmt psmeca beachballs_ipts1_iref1_lune_psmeca -JXi/6.72i -R-30/30/-90/90 -Sm0.45/8p -L0.5p,0/0/0 -W0.5p,0/0/0 -G255/0/0 -N -O -V >> lune_gmt600.ps
