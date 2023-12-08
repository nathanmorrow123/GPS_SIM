#!/bin/bash
#Download current yuma and opsadvisories off of NAVCENT
wget https://www.navcen.uscg.gov/sites/default/files/gps/opsadvisory/current_opsadvisory.txt -O current_opsadvisory.txt
wget https://www.navcen.uscg.gov/sites/default/files/gps/almanac/current_yuma.alm -O current_yuma.alm 

