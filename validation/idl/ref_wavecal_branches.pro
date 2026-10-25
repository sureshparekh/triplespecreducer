; The wavecal branches ref_wavecal.pro does not reach: a single sky frame
; (no median), and a manual offset (Auto X-Correlate off: XGETOFFSET's FLOAT).
@setup
mc_readflat,out+'flat.fits',flat,ncols,nrows,obsmode,ps,slith_pix,slith_arc, $
            slitw_pix,slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms, $
            rotation,edgedeg,OMASK=omask,CANCEL=cancel
modeinfo = mc_readmodefile(pkg+'data/XD.dat')
for run = 0, 1 do begin & $
  if run eq 0 then arcfiles = raw+'SPEC_SO2025D-006-0120.fits' else arcfiles = raw+'SPEC_SO2025D-006-'+string(indgen(12)+120,FORMAT='(I4.4)')+'.fits' & $
  oname = (run eq 0) ? 'wcb_single' : 'wcb_manual' & $
  kna_readts4fits,arcfiles,data,hdrinfo,var,ROTATE=rotation,KEYWORDS=keywords,NIMAGES=nimages & $
  if nimages ne 1 then arc = median(data,/EVEN,DIMENSION=3) else arc = reform(data) & $
  arc = temporary(arc)/flat & $
  fake = mc_simwavecal2d(ncols,nrows,edgecoeffs,xranges,slith_arc,CANCEL=cancel) & $
  if run eq 0 then writefits, out+'wcb_simwave.fits', fake & $
  spec = kna_medextspec1d(data,var,omask,orders,fake[*,*,0],fake[*,*,1],[1],APRANGES=rebin([slith_arc/2.-5,slith_arc/2+5],2,n_elements(orders)),CANCEL=cancel) & $
  calinfo = mc_readlinecalfile(pkg+'data/OHXD_LineCal.fits',CANCEL=cancel) & $
  whome = reform(calinfo.xcorspec[*,0]) & fhome = reform(calinfo.xcorspec[*,1]) & $
  idx = mc_nantrim(whome,2,CANCEL=cancel) & whome=whome[idx] & fhome=fhome[idx] & $
  z = where(orders eq calinfo.xcororder) & $
  wobs = reform((spec.(z))[*,0]) & fobs = reform((spec.(z))[*,1]) & $
  linterp,wobs,fobs,whome,ifobs,MISSING=0 & $
  if run eq 0 then begin & $
     lag = findgen(ncols)-fix(ncols/2.) & corr = c_correlate(fhome,ifobs,lag) & $
     max = max(corr,idx) & fit = gaussfit(lag[idx-20:idx+20],corr[idx-20:idx+20],a,NTERMS=4) & offset = a[1] & $
  endif else offset = 1.25 - 1.9 & $
  print, 'OFFSET', offset & help, offset & $
  readcol,pkg+'data/OHll.dat',lines,FORMAT='D',COMMENT='#',/SILENT & $
  linetypes = replicate('OH',n_elements(lines)) & $
  coeffs = mc_xdwavecal1d(spec,orders,calinfo.w2pcoeffs,offset,lines,linetypes,calinfo.ordertype,slitw_pix,calinfo.homeorder,modeinfo.dispdeg,modeinfo.ordrdeg,3.5,out+oname+'_Wavecal',NLINES=nlines,WRANGES=wranges,RMS=rms1d,CANCEL=cancel) & $
  print, 'NLINES', nlines & $
  wcinfo1d = {coeffs:coeffs,dispdeg:modeinfo.dispdeg,ordrdeg:modeinfo.ordrdeg,homeorder:modeinfo.homeorder} & $
  mc_rdohll,pkg+'data/OHll.dat',plines,RES=rp,THRESH=3,CANCEL=cancel & $
  wcinfo2d = mc_findlines2d_knahack(arc,edgecoeffs,xranges,orders,plines,calinfo.w2pcoeffs,offset,wranges,slitw_pix,slith_arc,modeinfo.linedeg,modeinfo.linereg/2.,modeinfo.ystep,modeinfo.ysum,CANCEL=cancel) & $
  for i=0,norders-1 do writefits, out+oname+'_lines2d_'+strtrim(orders[i],2)+'.fits', [[wcinfo2d.(i).sx],[wcinfo2d.(i).sy],[wcinfo2d.(i).ss],[wcinfo2d.(i).sl],[wcinfo2d.(i).sw]] & $
  mc_wavecal2d_knahack,wcinfo1d,wcinfo2d,edgecoeffs,omask,orders,xranges,wranges,rotation,modeinfo.wxdeg,modeinfo.wydeg,modeinfo.sxdeg,modeinfo.sydeg,out+oname+'.fits',modeinfo.ps,out+oname+'_Distortion.ps',CANCEL=cancel & $
endfor
print, 'DONE'
exit
