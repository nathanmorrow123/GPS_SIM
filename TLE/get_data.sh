#!/bin/bash
#Download current GNSS TLE from celestrak.org
wget https://celestrak.org/NORAD/elements/gp.php?GROUP=gnss&FORMAT=tle
mv gp.php?GROUP=gnss gnss_tle.txt
sudo rm 'gp.php?GROUP=gnss'


