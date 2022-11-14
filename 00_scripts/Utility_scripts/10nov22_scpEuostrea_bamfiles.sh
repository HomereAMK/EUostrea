DIRBAM=/home/projects/dp_00007/people/hmon/Shucking/06_realigned
##AGAB
    for POP in AGAB
    do
        for IND in `echo -n 01 02 03 04 05 06 08 10 11 12 13 14 15 18 19 20 23 24`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##BARR
    for POP in BARR
    do
        for IND in `echo -n 01 02 03 04 06 07 08 09 10 11 12 13 14 15`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##BUNN
    for POP in BUNN
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##CLEW
    for POP in CLEW
    do
        for IND in `echo -n 02 03 04 05 06 07 08 09 10`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##CORS
    for POP in CORS
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##COLN
    for POP in COLN
    do
        for IND in `echo -n 01 02 03 04 06 08 09 10 12 13 15 17 18 19 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##CRES
    for POP in CRES
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##DOLV
    for POP in DOLV
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##GREV
    for POP in GREV
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##HAFR
    for POP in HAFR
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08 09 10 11  13 14 15 16 17`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `echo -n 18 19 23 25`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##HALS
    for POP in HALS
    do
        for IND in `seq -w 01 15`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##HAUG
    for POP in HAUG
    do
        for IND in `seq -w 01 14`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done
    for POP in HAUG
    do
        for IND in `echo -n 16 17 18 19 20 21 24`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##HYPP
    for POP in HYPP
    do
        for IND in `seq -w 01 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##INNE
    for POP in INNE
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 17 19 20 21 22 25`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


#KALV
    for POP in KALV
    do
        for IND in `echo -n 01 02 03 04 05 06 07 08 09 10 12 13 16 18 19 20 21 23`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

    scp $DIRBAM/${POP}_11_EKDL210004* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea



##LOGS
    for POP in LOGS
    do
        for IND in `seq -w 01 11`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done
    scp $DIRBAM/${POP}_99* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
    for POP in LOGS
    do
        for IND in `seq -w 13 19`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##MOLU
    for POP in MOLU
    do
        for IND in `seq -w 01 19`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done
    scp $DIRBAM/${POP}_25* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea


##NISS
    for POP in NISS
    do
        for IND in `seq -w 02 10`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `echo -n 12 13 16 21 26 27 31 32 38 39 42 49`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##ORIS
    for POP in ORIS
    do
        for IND in `seq -w 01 15`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##OSTR
    for POP in OSTR
    do
        for IND in `seq -w 01 05`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done
    for POP in OSTR
    do
        for IND in `seq -w 07 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##PONT
    for POP in PONT
    do
        for IND in `seq -w 01 10`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `echo -n 13 15 16 17 19 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##RIAE
    for POP in RIAE
    do
        for IND in `seq -w 05 22`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `echo -n 03`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##RYAN
    for POP in RYAN
    do
        for IND in `seq -w 01 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##THIS
    for POP in THIS
    do
        for IND in `seq -w 01 19`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `echo -n 20 21 23 24 26 `
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##TOLL
    for POP in TOLL
    do
        for IND in `seq -w 01 09`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `seq -w 12 18`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done
    

##TRAL
    for POP in TRAL
    do
        for IND in `seq -w 01 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##VAGS
    for POP in VAGS
    do
        for IND in `seq -w 01 06`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `seq -w 08 22`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done


##VENO
    for POP in VENO
    do
        for IND in `seq -w 01 05`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `echo -n 07 08 09 10 11 12 14 15 16 17 18`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##WADD
    for POP in WADD
    do
        for IND in `seq -w 01 19`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done

##ZECE
    for POP in ZECE
    do
        for IND in `seq -w 01 08`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `seq -w 10 20`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
        for IND in `seq -w 22 24`
        do
        scp $DIRBAM/${POP}_${IND}* /home/projects/dp_00007/people/hmon/Bamfile_EUostrea
        done
    done