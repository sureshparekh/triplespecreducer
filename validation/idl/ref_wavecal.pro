; Reference wavecal, reproducing ts41cals_mkwavecal (kna_ts41cals2d.pro), Sky mode,
; with Spextool's default autoxcorr=1, without the display.
@setup
oname = 'wavecal_NGC134'
arcfiles = raw+'SPEC_SO2025D-006-'+string(indgen(12)+120,FORMAT='(I4.4)')+'.fits'
mc_readflat,out+'flat.fits',flat,ncols,nrows,obsmode,ps,slith_pix,slith_arc, $
            slitw_pix,slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms, $
            rotation,edgedeg,OMASK=omask,CANCEL=cancel
modeinfo = mc_readmodefile(pkg+'data/XD.dat')
kna_readts4fits,arcfiles,data,hdrinfo,var,ROTATE=rotation,KEYWORDS=keywords,NIMAGES=nimages
if nimages ne 1 then arc = median(data,/EVEN,DIMENSION=3) else arc = reform(data)
arc = temporary(arc)/flat
writefits, out+'wc_arc.fits', arc
fake = mc_simwavecal2d(ncols,nrows,edgecoeffs,xranges,slith_arc,CANCEL=cancel)
spec = kna_medextspec1d(data,var,omask,orders,fake[*,*,0],fake[*,*,1],[1],$
                        APRANGES=rebin([slith_arc/2.-5,slith_arc/2+5],2,n_elements(orders)),CANCEL=cancel)
for i=0,n_elements(orders)-1 do writefits, out+'wc_spec1d_'+strtrim(orders[i],2)+'.fits', spec.(i)
calinfo = mc_readlinecalfile(pkg+'data/OHXD_LineCal.fits',CANCEL=cancel)
whome = reform(calinfo.xcorspec[*,0]) & fhome = reform(calinfo.xcorspec[*,1])
idx = mc_nantrim(whome,2,CANCEL=cancel) & whome=whome[idx] & fhome=fhome[idx]
z = where(orders eq calinfo.xcororder)
wobs = reform((spec.(z))[*,0]) & fobs = reform((spec.(z))[*,1])
linterp,wobs,fobs,whome,ifobs,MISSING=0
lag = findgen(ncols)-fix(ncols/2.)
corr = c_correlate(fhome,ifobs,lag)
max = max(corr,idx)
fit = gaussfit(lag[idx-20:idx+20],corr[idx-20:idx+20],a,NTERMS=4)
offset = a[1]
print, 'OFFSET', offset
readcol,pkg+'data/OHll.dat',lines,FORMAT='D',COMMENT='#',/SILENT
linetypes = replicate('OH',n_elements(lines))
coeffs = mc_xdwavecal1d(spec,orders,calinfo.w2pcoeffs,offset,lines, $
                        linetypes,calinfo.ordertype,slitw_pix, $
                        calinfo.homeorder,modeinfo.dispdeg, $
                        modeinfo.ordrdeg,3.5,out+oname+'_Wavecal',$
                        NLINES=nlines,WRANGES=wranges,RMS=rms1d,CANCEL=cancel)
print, 'COEFFS1D', coeffs
print, 'NLINES', nlines
print, 'RMS1D', rms1d
print, 'WRANGES', wranges
wcinfo1d = {coeffs:coeffs,dispdeg:modeinfo.dispdeg,ordrdeg:modeinfo.ordrdeg,homeorder:modeinfo.homeorder}
mc_rdohll,pkg+'data/OHll.dat',plines,RES=rp,THRESH=3,CANCEL=cancel
print, 'NPLINES', n_elements(plines)
wcinfo2d = mc_findlines2d_knahack(arc,edgecoeffs,xranges,orders,plines, $
                        calinfo.w2pcoeffs,offset,wranges,slitw_pix, $
                        slith_arc,modeinfo.linedeg,modeinfo.linereg/2., $
                        modeinfo.ystep,modeinfo.ysum,CANCEL=cancel)
for i=0,norders-1 do writefits, out+'wc_lines2d_'+strtrim(orders[i],2)+'.fits', [[wcinfo2d.(i).sx],[wcinfo2d.(i).sy],[wcinfo2d.(i).ss],[wcinfo2d.(i).sl],[wcinfo2d.(i).sw]]
mc_wavecal2d_knahack,wcinfo1d,wcinfo2d,edgecoeffs,omask,orders,xranges,wranges,rotation, $
             modeinfo.wxdeg,modeinfo.wydeg,modeinfo.sxdeg,modeinfo.sydeg, $
             out+oname+'.fits',modeinfo.ps,out+oname+'_Distortion.ps',/SILENT,CANCEL=cancel
print, 'DONE'
exit
