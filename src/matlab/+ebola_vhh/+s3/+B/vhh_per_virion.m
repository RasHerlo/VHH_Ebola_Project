function fig = vhh_per_virion()
%VHH_PER_VIRION  AF647-VHH copies per MeGFP-VSV at 44 nM (G10, G84, G68).
    fig = ebola_vhh.render('s3b_vhh_per_virion', 's3', 'B', 'vhh_per_virion');
end
