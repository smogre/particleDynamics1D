
    echo param.nfusedist_pf$1; 
    # copy template file to all files
    cp param.fixtotlen2x_2mt param.fixtotlen${j}x_${i}mt;
    # set number of microtubules
    sed -i "s/NMT 2/NMT ${i}/g" param.fixtotlen${j}x_${i}*;
    # set filenames to read MT lengths and nucleation points
    sed -i "s/40um/$((${j}*20))um/g" param.fixtotlen${j}x_${i}*;
    # set tip capture property
    sed -i "s/TIPCAPTURE T/TIPCAPTURE F/g" param.fixtotlen${j}x_${i}*;
    # set stop on capture property
    sed -i "s/STOPONCAPTURE T/STOPONCAPTURE F/g" param.fixtotlen${j}x_${i}*;
    # set random nucleation property
    #sed -i "s/RANDNUC F/RANDNUC T/g" param.fixtotlen${j}x_${i}*;
    # set READNUC
    #sed -i "s/READNUC F/READNUC T/g" param.fixtotlen${j}x_${i}*;
