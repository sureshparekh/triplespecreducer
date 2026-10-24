!path=expand_path('+'+getenv('SPEXTOOL_DIR'))+':'+!path
arc = readfits(getenv('TSR_ROOT')+'/validation/idl_out/wc_arc.fits',/SILENT)
x = findgen(2048) & guessx = 402.110306423524d & midslit = 427.132572279766d & k = 6 & ystep = 3 & ysum = 2 & linewin = 5.0
xsub = x[((guessx-linewin) > 0):((guessx+linewin) < (2048-1))]
ysub = median(arc[((guessx-linewin) > 0):((guessx+linewin) < (2048-1)), (midslit+k*ystep-ysum):midslit+k*ystep+ysum], DIMEN=2)
tabinv,xsub,guessx,zidx & guessz = ysub[zidx]
yfit = mpfitpeak(double(xsub),double(ysub),a,NTERMS=4,ESTIMATES=double([guessz,guessx,2.4/2.354,min(ysub)]))
print, 'IDL', a, FORMAT='(A,4F16.9)'
print, 'YS', ysub
print, 'SIZE', n_elements(xsub), xsub[0]
exit
