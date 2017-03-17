#!/bin/bash

# Extract coordinates for the current meters (which only have names like "Site 1").

for d in *; do
    if [ ! -d $d ]; then
        continue
    fi
    cd $d
    for i in *.m21fm; do
        # Do a more sophisticated version which cuts all the time series
        # outputs I specified.
        # OUTPUT_2 = BODC tides
        # OUTPUT_3 = Current meters
        # OUTPUT_4 = North Sea Tides
        names=("tide" "points" "north_sea_points")
        outputs=("OUTPUT_2" "OUTPUT_3" "OUTPUT_4")
        for ((n=0; n<${#names[@]}; n++)); do
            awk '/\['${outputs[n]}'\]/{flag=1;next}/EndSect  \/\/ '${outputs[n]}'/{flag=0}flag' $i | \
                tr -d '\015' | \
                grep --group-separator="" -A4 POINT | \
                grep . | \
                grep -v POINT | \
                sed '/\[LINE\]/,$d' | \
                cut -f2 -d= | \
                sed 's/Site /Site_/g' | \
                xargs -n4 | \
                tr " " "," > "${i%.*}_${names[n]}.txt"
                if [ $(grep . "${i%.*}_${names[n]}.txt" | wc -l) -eq 0 ]; then
                    rm "${i%.*}_${names[n]}.txt"
                fi
            done
        done


        #echo "name,lon,lat,depth" > "${i%.*}_current_meter_sites.txt"
        ## Holy horrible pipeline Batman!
        #grep --group-separator="" Site\ .* -A3 "$i" | \
        #    tr -d '\015' | \
        #    sed 's/Site /Site_/g' | \
        #    grep . | \
        #    cut -f2 -d= | \
        #    xargs -n4 | \
        #    tr " " "," >> "${i%.*}_current_meter_sites.txt"
        #if [ $(grep . "${i%.*}_current_meter_sites.txt" | wc -l) -eq 1 ]; then
        #    rm "${i%.*}_current_meter_sites.txt"
        #fi
        #done
    cd ~-
done

