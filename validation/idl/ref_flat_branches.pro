; The flat branches ref_flat.pro does not reach: a single lamps-on frame
; (no scaling, no combining), a single dark, and no darks.
@setup
on5 = raw+'SPEC_DFLAT_'+string(indgen(5)+2,FORMAT='(I4.4)')+'.fits'
on1 = raw+'SPEC_DFLAT_0006.fits'
off1 = raw+'SPEC_DFLAT_0060.fits'
flatcase, on1, '', 'flatb_one_nodark', out, pkg, keywords
flatcase, on5, off1, 'flatb_five_onedark', out, pkg, keywords
flatcase, on5, '', 'flatb_five_nodark', out, pkg, keywords
flatcase, on1, off1, 'flatb_one_onedark', out, pkg, keywords
print, 'DONE'
exit
