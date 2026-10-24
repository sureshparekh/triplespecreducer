!path=expand_path('+'+getenv('SPEXTOOL_DIR'))+':'+!path
c = readfits('../idl_out/mp_cases.fits',/SILENT)
for i=0,2 do begin & v = c[*,i] & xs = double(v[0:10]) & ys = double(v[11:21]) & e = double(v[22:25]) & yf = mpfitpeak(xs,ys,a,NTERMS=4,ESTIMATES=e) & print, 'IDL', a, FORMAT='(A,4F16.9)' & endfor
exit
