DIRBAM=/home/projects/dp_00007/people/hmon/Shucking/06_realigned

    for POP in AGAB
    do
        for IND in `echo -n 01 02 03 04 05 06 08 10 11 12 13 14 15 18 19 20 23 24`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done
