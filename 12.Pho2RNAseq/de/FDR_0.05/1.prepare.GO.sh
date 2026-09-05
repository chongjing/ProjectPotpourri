#!/bin/bash

cd /home/cx264/rds/hpc-work/project/2.Jeongmin/01.Saskia/4.DEanalysis.FDR_0.05/
for folder in */; do
    folder_name="${folder%/}"
    cd $folder_name
#    awk -F"," '$8 == "UP" {print $2}' 1.5.QLF.DE.results.*.FINAL.csv | sed 's/AGIS_//g' > 2.1.$folder_name\.UP.list
#    awk -F"," '$8 == "DOWN" {print $2}' 1.5.QLF.DE.results.*.FINAL.csv | sed 's/AGIS_//g' > 2.1.$folder_name\.DOWN.list
#    awk -F"," '$8 == "Not_significant" {print $2}' 1.5.QLF.DE.results.*.FINAL.csv | sed 's/AGIS_//g' > 2.1.$folder_name\.Non-Sig.list

#    while read line; do awk -v line="$line" '$1 ~ line {print $0}' /rds/user/cx264/hpc-work/project/0.ref/1.Rice_Nipponbare/NIP-T2T_aa.GO.list >> 2.2.$folder_name\.UP-GO.list; done < 2.1.$folder_name\.UP.list
#    while read line; do awk -v line="$line" '$1 ~ line {print $0}' /rds/user/cx264/hpc-work/project/0.ref/1.Rice_Nipponbare/NIP-T2T_aa.GO.list >> 2.2.$folder_name\.DOWN-GO.list; done < 2.1.$folder_name\.DOWN.list
#    while read line; do awk -v line="$line" '$1 ~ line {print $0}' /rds/user/cx264/hpc-work/project/0.ref/1.Rice_Nipponbare/NIP-T2T_aa.GO.list >> 2.2.$folder_name\.Non-Sig-GO.list; done < 2.1.$folder_name\.Non-Sig.list

    ## only need gene names for GO enrichment
    awk '{print $1}' 2.2.$folder_name\.UP-GO.list > 2.3.$folder_name\.UP-GO.list
    awk '{print $1}' 2.2.$folder_name\.DOWN-GO.list > 2.3.$folder_name\.DOWN-GO.list
    awk '{print $1}' 2.2.$folder_name\.Non-Sig-GO.list > 2.3.$folder_name\.Non-Sig-GO.list
    cd ../
done

