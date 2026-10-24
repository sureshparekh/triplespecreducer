; Reference extractions reproducing xspextool (2D wavecal mode) without the display:
; point source (HD225200 116-117, A-B) and extended source (NGC134 120-121, A-B).
@setup
mc_readflat, out+'flat.fits', flat,ncols,nrows,modename,ps,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms,rotation,OMASK=omask
mc_rdwavecal2d, out+'wavecal_NGC134.fits', wavecal, spatcal, wctype, xo2w, wdeg, odeg, homeorder, xy2w, xy2s, wxdeg, wydeg, sxdeg, sydeg, disp, ROTATE=rotation
atr = readfits(sp+'/data/atran4000.fits',/SILENT) & awave = reform(atr[*,0]) & atrans = reform(atr[*,1])
print, 'PS', ps, ' SLITH', slith_arc, slith_pix
; ---------------- point source
files = raw+['SPEC_SO2025D-006-0116.fits','SPEC_SO2025D-006-0117.fits']
kna_readts4fits, files, data, hdrinfo, var, KEYWORD=keywords, ROTATE=rotation, /PAIR, ITOT=itot
data = data/flat & var = var/flat^2
writefits, out+'ps_work.fits', data
prof = mc_mkspatprof2d(data, omask, wavecal, spatcal, slith_arc/slith_pix, 1, awave, atrans, 0.4, SPATCOEFFS=spatcoeffs, MEDSPATCOEFFS=medspatcoeffs)
for i=0,norders-1 do writefits, out+'ps_prof_'+strtrim(orders[i],2)+'.fits', prof.(i)
doorders = intarr(norders)+1
mc_findpeaks, prof, 2, doorders, positions, apsign, /AUTO
print, 'PSPOS', positions & print, 'PSSIGN', apsign
tracecoeffs = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs, xranges, positions, 5, 5, 3.0, 0.1, 1)
print, 'TRACE', tracecoeffs
apsigndat = fltarr(2)
for i=0,norders-1 do for j=0,1 do begin & profile = prof.(i) & med = median(profile[*,1]) & tabinv,profile[*,0],positions[j,i],idx & apsigndat[j] = (profile[idx,1] gt med) ? 1:-1 & endfor
apsign = fix(apsigndat)
print, 'APSIGN', apsign
apradii = replicate(1.0,2,norders)
s = mc_extspec2d(data,var,omask,orders,wavecal,spatcal,ps,apsign,tracecoeffs,apradii,PSBGINFO=[2.5,5.0],BGDEG=1)
for k=0,n_tags(s)-1 do writefits, out+'ps_spec_'+strtrim(k,2)+'.fits', s.(k)
writefits, out+'ps_positions.fits', positions & writefits, out+'ps_trace.fits', tracecoeffs
; ---------------- extended source (Luis's tutorial settings)
files = raw+['SPEC_SO2025D-006-0120.fits','SPEC_SO2025D-006-0121.fits']
kna_readts4fits, files, data, hdrinfo, var, KEYWORD=keywords, ROTATE=rotation, /PAIR, ITOT=itot
data = data/flat & var = var/flat^2
prof = mc_mkspatprof2d(data, omask, wavecal, spatcal, slith_arc/slith_pix, 1, awave, atrans, 0.4)
for i=0,norders-1 do writefits, out+'xs_prof_'+strtrim(orders[i],2)+'.fits', prof.(i)
pos = rebin([20.1],1,norders)
coeffs = fltarr(2,norders) & for i=0,norders-1 do coeffs[*,i] = [pos[0,i],0]
apsign = intarr(1)+1
bgr = [[2.,6.],[32.,36.]]
apradii = rebin([1.0],1,norders)
s = mc_extspec2d(data,var,omask,orders,wavecal,spatcal,ps,apsign,coeffs,apradii,XSBGINFO=bgr,BGDEG=bgorder)
for k=0,n_tags(s)-1 do writefits, out+'xs_spec_'+strtrim(k,2)+'.fits', s.(k)
print, 'DONE'
exit
