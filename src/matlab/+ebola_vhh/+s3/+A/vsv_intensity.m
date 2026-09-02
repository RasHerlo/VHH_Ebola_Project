function fig = vsv_intensity()
%VSV_INTENSITY  MeGFP intensity of VSV particles on glass (SDCM).
    fig = ebola_vhh.render('s3a_vsv_intensity', 's3', 'A', 'vsv_intensity');
end
