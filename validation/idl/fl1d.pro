@setup
calinfo = mc_readlinecalfile(pkg+'data/OHXD_LineCal.fits')
readcol,pkg+'data/OHll.dat',lines,FORMAT='D',COMMENT='#',/SILENT
s = readfits(out+'wc_spec1d_6.fits',/SILENT)
wspec = reform(s[*,0]) & fspec = reform(s[*,1])
idx = mc_nantrim(wspec,2) & wspec=wspec[idx] & fspec=fspec[idx]
xguess = poly(lines,calinfo.w2pcoeffs[*,3]) - 0.656294
z = where(xguess gt min(wspec,MAX=xmax) and xguess lt xmax,cnt)
xguess = xguess[z] & plines = lines[z]
xpos = mc_findlines1d(wspec,fspec,xguess,2.4/2.354,WLINES=plines,GOODBAD=goodbad)
for i=0,n_elements(xguess)-1 do print, plines[i], xguess[i], xpos[i], goodbad[i], FORMAT='(F12.8,2F12.4,I3)'
exit
