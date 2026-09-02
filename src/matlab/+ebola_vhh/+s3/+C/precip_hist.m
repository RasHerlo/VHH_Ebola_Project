function fig = precip_hist()
%PRECIP_HIST  Off-cell precipitation with no VHH vs G10 vs G84.
    fig = ebola_vhh.render('s3c_precip_hist', 's3', 'C', 'precip_hist');
end
