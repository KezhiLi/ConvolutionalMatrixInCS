function snr_val=snr2(x0,x1)
err=norm(x0-x1);
sig=norm(x0);
snr_val=-20*log10(err/sig);
