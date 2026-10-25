; Reference for the xmergeorders / xcleanspec / spex2text / xgetoffset options
; (validation/compare_mergeclean.py, tests/test_parity.py).
; xmergeorders: its own data procedures (XMERGEORDERS_CUTSCALESPEC,
; XMERGEORDERS_AUTOSCALE through tsr_autoscale, XMERGEORDERS_COMBINESPEC) on a
; state structure, the events (add order, 's' region, trims, scale field,
; merge, write) made by hand as the widget does them.
@setup
!path = getenv('TSR_ROOT')+'/validation/idl:'+!path
resolve_routine, 'xmergeorders', /COMPILE_FULL_FILE, /EITHER
nan = !values.f_nan
cs = readfits(out+'tc_NGC134.fits',/SILENT)
cs = mc_caaspecflags(cs)
orders = [3,4,5,6,7] & naps = 1 & ap = 0
state = {r:{trim:nan, trimdir:'Right', trimspec:'Green', scale:1., weighted:1, sum:0, cursormode:'None'}, d:{addspec:ptr_new(2), bluespec:ptr_new(2), greenspec:ptr_new(2), mergedspec:ptr_new(2), tempmergedspec:ptr_new(2)}, p:{reg:[[nan,nan],[nan,nan]]}, w:{scale_fld:[0L,0L], xmergeorders_base:0L}}
specscales = fltarr(5)+1.0
; anchor order 4 (XMERGEORDERS_MODANCHOR)
index = 1*naps+ap & nonan = where(finite(cs[*,1,index]) eq 1)
*state.d.bluespec = reform(cs[nonan,*,index]) & *state.d.mergedspec = reform(cs[nonan,*,index])
mergedorders = [4]
; ---- add 5: auto scale over 1.42-1.47, then a green cut 'Right' at 1.47 and a blue cut 'Left' at 1.44
addorder = 5 & z = where(orders eq addorder) & z = total(z)*naps+ap & nonan = where(finite(cs[*,1,z]) eq 1)
*state.d.addspec = reform(cs[nonan,*,z]) & *state.d.greenspec = reform(cs[nonan,*,z]) & *state.d.bluespec = *state.d.mergedspec & state.p.reg = nan
state.p.reg[*,0] = [1.42,0.] & state.p.reg[*,1] = [1.47,0.]
tsr_autoscale, state & xmergeorders_cutscalespec, state
print, 'SCALE5', state.r.scale, format='(a,e24.16)'
state.r.trimspec = 'Green' & state.r.trimdir = 'Right' & state.r.trim = 1.47 & xmergeorders_cutscalespec, state & state.r.trim = nan
state.r.trimspec = 'Blue' & state.r.trimdir = 'Left' & state.r.trim = 1.44 & xmergeorders_cutscalespec, state & state.r.trim = nan
xmergeorders_combinespec, state, CANCEL=cancel
mergedorders = [mergedorders,addorder] & *state.d.bluespec = reform(*state.d.tempmergedspec) & *state.d.mergedspec = reform(*state.d.tempmergedspec)
writefits, out+'mc_step1.fits', *state.d.mergedspec
specscales[where(orders eq addorder)] = state.r.scale & state.r.scale = 1.0 & state.p.reg[*] = nan & state.r.trim = nan
; ---- add 3: no overlap, the scale field set to 1.05
addorder = 3 & z = where(orders eq addorder) & z = total(z)*naps+ap & nonan = where(finite(cs[*,1,z]) eq 1)
*state.d.addspec = reform(cs[nonan,*,z]) & *state.d.greenspec = reform(cs[nonan,*,z]) & *state.d.bluespec = *state.d.mergedspec & state.p.reg = nan
state.r.scale = float('1.05') & xmergeorders_cutscalespec, state
xmergeorders_combinespec, state, CANCEL=cancel
mergedorders = [mergedorders,addorder] & *state.d.bluespec = reform(*state.d.tempmergedspec) & *state.d.mergedspec = reform(*state.d.tempmergedspec)
writefits, out+'mc_step2.fits', *state.d.mergedspec
specscales[where(orders eq addorder)] = state.r.scale & state.r.scale = 1.0 & state.p.reg[*] = nan & state.r.trim = nan
; ---- add 6: region selected right to left, 1.22 then 1.15 (it reaches into both cuts);
;      auto scale, then a green cut 'Right' at 1.20 and a blue cut 'Left' at 1.16
addorder = 6 & z = where(orders eq addorder) & z = total(z)*naps+ap & nonan = where(finite(cs[*,1,z]) eq 1)
*state.d.addspec = reform(cs[nonan,*,z]) & *state.d.greenspec = reform(cs[nonan,*,z]) & *state.d.bluespec = *state.d.mergedspec & state.p.reg = nan
state.p.reg[*,0] = [1.22,0.] & state.p.reg[*,1] = [1.15,0.]
tsr_autoscale, state & xmergeorders_cutscalespec, state
print, 'SCALE6', state.r.scale, format='(a,e24.16)'
state.r.trimspec = 'Green' & state.r.trimdir = 'Right' & state.r.trim = 1.20 & xmergeorders_cutscalespec, state & state.r.trim = nan
state.r.trimspec = 'Blue' & state.r.trimdir = 'Left' & state.r.trim = 1.16 & xmergeorders_cutscalespec, state & state.r.trim = nan
xmergeorders_combinespec, state, CANCEL=cancel
mergedorders = [mergedorders,addorder] & *state.d.bluespec = reform(*state.d.tempmergedspec) & *state.d.mergedspec = reform(*state.d.tempmergedspec)
writefits, out+'mc_step3.fits', *state.d.mergedspec
specscales[where(orders eq addorder)] = state.r.scale & state.r.scale = 1.0 & state.p.reg[*] = nan & state.r.trim = nan
; ---- XMERGEORDERS_WRITEFILE after orders 4, 5, 3, 6 (header of the combined file, which has the Spextool cards)
hdr = headfits(out+'comb_NGC134.fits') & file = 'mc_merged' & ifile = out+'tc_NGC134.fits'
stop_ = n_elements((*state.d.mergedspec)[*,1])-1
fxaddpar,hdr,'IRAFNAME',file & fxaddpar,hdr,'ORDERS','0' & fxaddpar,hdr,'NORDERS',1 & fxaddpar,hdr,'NAPS',1 & fxaddpar,hdr,'START',0 & fxaddpar,hdr,'STOP',stop_
history = 'This merged spectrum was derived from aperture '+string(ap+1,FORMAT='(i2.2)')+' of file '+strtrim(file_basename(ifile),2)+'.'
morders = mergedorders[sort(mergedorders)] & match,morders,orders,junk,idx
history = history+'  The scaled factors for orders '+strjoin(strcompress(morders,/RE),', ')+' are '+strjoin(strcompress(specscales[idx],/RE),', ')+' respectively.'
hdr = mc_addhistlabel(hdr,'Xmergeorders History') & history = mc_splittext(history,70) & sxaddhist,history,hdr
writefits, out+file+'.fits', *state.d.mergedspec, hdr
; ---- add 7: a green cut 'Left' at 0.97 puts the order wholly inside the merged range (MC_MERGESPEC 'Inside')
addorder = 7 & z = where(orders eq addorder) & z = total(z)*naps+ap & nonan = where(finite(cs[*,1,z]) eq 1)
*state.d.addspec = reform(cs[nonan,*,z]) & *state.d.greenspec = reform(cs[nonan,*,z]) & *state.d.bluespec = *state.d.mergedspec & state.p.reg = nan
state.r.trimspec = 'Green' & state.r.trimdir = 'Left' & state.r.trim = 0.97 & xmergeorders_cutscalespec, state & state.r.trim = nan
writefits, out+'mc_inside_blue.fits', *state.d.bluespec & writefits, out+'mc_inside_green.fits', *state.d.greenspec
xmergeorders_combinespec, state, CANCEL=cancel
mergedorders = [mergedorders,addorder] & *state.d.bluespec = reform(*state.d.tempmergedspec) & *state.d.mergedspec = reform(*state.d.tempmergedspec)
writefits, out+'mc_step4.fits', *state.d.mergedspec
writefits, out+'mc_inside.fits', *state.d.mergedspec
; ---- xcleanspec on that file: XCLEANSPEC_LOADSPEC, then edits (plot window events),
;      's', 'u', XCLEANSPEC_SMOOTHSPEC, XCLEANSPEC_WRITEFILE
mc_readspec,out+'mc_merged.fits',spc,chdr,obsmode,start,stop,norders,naps,corders,xunits,yunits,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,airmass,xtitle,ytitle
spc = mc_caaspecflags(spc)
espec = spc & pspec = spc & pidx = 0
sw = float(slitw_pix) & sc = (slitw_pix gt 2) ? 1.5:2.5
sgwin_txt = strtrim(string(slitw_pix*sc,FORMAT='(f4.1)'),2) & fwhm_txt = strtrim(sw,2)
print, 'DEFAULTS ', sgwin_txt, ' ', fwhm_txt
; Fix 1.27-1.29 (Flux), Fix 1.285-1.30 (Flux), no 's' between
xr = [1.27,1.29] & mc_tabinv,espec[*,0,pidx],xr,ii & ii = round(ii)
x = [espec[ii[0],0,pidx],espec[ii[1],0,pidx]] & y = [espec[ii[0],1,pidx],espec[ii[1],1,pidx]] & e = [espec[ii[0],2,pidx],espec[ii[1],2,pidx]]
coeff = poly_fit1d(x,y,1,/SILENT) & pspec[ii[0]:ii[1],1,pidx] = poly(espec[ii[0]:ii[1],0,pidx],coeff)
coeff = poly_fit1d(x,e,1,/SILENT) & pspec[ii[0]:ii[1],2,pidx] = poly(espec[ii[0]:ii[1],0,pidx],coeff) & pspec[ii[0]:ii[1],3,pidx] = pspec[ii[0]:ii[1],3,pidx] + 4
xr = [1.30,1.285] & xr = xr[sort(xr)] & mc_tabinv,espec[*,0,pidx],xr,ii & ii = round(ii)
x = [espec[ii[0],0,pidx],espec[ii[1],0,pidx]] & y = [espec[ii[0],1,pidx],espec[ii[1],1,pidx]] & e = [espec[ii[0],2,pidx],espec[ii[1],2,pidx]]
coeff = poly_fit1d(x,y,1,/SILENT) & pspec[ii[0]:ii[1],1,pidx] = poly(espec[ii[0]:ii[1],0,pidx],coeff)
coeff = poly_fit1d(x,e,1,/SILENT) & pspec[ii[0]:ii[1],2,pidx] = poly(espec[ii[0]:ii[1],0,pidx],coeff) & pspec[ii[0]:ii[1],3,pidx] = pspec[ii[0]:ii[1],3,pidx] + 4
espec = pspec
; Remove 1.82-1.80 (Error), 's'
reg = [1.82,1.80] & mc_tabinv,espec[*,0,pidx],reg,ii & ii = round(ii) & mn = ii[0] < ii[1] & mx = ii[0] > ii[1]
aa = findgen(n_elements(pspec[*,0,0])) & amin = min(aa,MAX=amax) & pspec[(mn > amin):(mx < amax),2,pidx] = !values.f_nan
espec = pspec
; Fix 1.55-1.56 (Flux), then 'u'
xr = [1.55,1.56] & mc_tabinv,espec[*,0,pidx],xr,ii & ii = round(ii)
x = [espec[ii[0],0,pidx],espec[ii[1],0,pidx]] & y = [espec[ii[0],1,pidx],espec[ii[1],1,pidx]] & coeff = poly_fit1d(x,y,1,/SILENT) & pspec[ii[0]:ii[1],1,pidx] = poly(espec[ii[0]:ii[1],0,pidx],coeff) & pspec[ii[0]:ii[1],3,pidx] = pspec[ii[0]:ii[1],3,pidx] + 4
pspec = espec
; Fix 2.10-2.12 (Error), not committed: the smoothing starts from espec
xr = [2.10,2.12] & mc_tabinv,espec[*,0,pidx],xr,ii & ii = round(ii)
x = [espec[ii[0],0,pidx],espec[ii[1],0,pidx]] & e = [espec[ii[0],2,pidx],espec[ii[1],2,pidx]] & coeff = poly_fit1d(x,e,1,/SILENT) & pspec[ii[0]:ii[1],2,pidx] = poly(espec[ii[0]:ii[1],0,pidx],coeff) & pspec[ii[0]:ii[1],3,pidx] = pspec[ii[0]:ii[1],3,pidx] + 4
writefits, out+'cl_edits.fits', espec
; XCLEANSPEC_SMOOTHSPEC, Savitzky-Golay, the default fields
spec = espec & sgwin = float(sgwin_txt) & sgdeg = float('2')
sgdeg = mc_crange(sgdeg,sgwin,'SG Degree',/KLT,CANCEL=cancel)
for i = 0, norders-1 do for j = 0, naps-1 do begin & z = where(finite(espec[*,1,i*naps+j]) eq 1) & cflux = espec[z,1,i*naps+j] & cerror = espec[z,2,i*naps+j] & cflux = savitzky_golay(cflux,sgwin,DEGREE=sgdeg,ERROR=cerror) & spec[z,1,i*naps+j] = cflux & spec[z,2,i*naps+j] = cerror & endfor
smthhistory='This spectrum has been convolved with a Savitzky-Golay filter of width '+strtrim(sgwin,2)+' pixels and degree '+strtrim(sgdeg,2)+'.'
print, "SGHIST|"+smthhistory+"|"
espec = spec & pspec = spec
whdr = chdr & whdr = mc_addhistlabel(whdr,'Xcleanspec History') & sxaddhist,mc_splittext(smthhistory,70),whdr
writefits, out+'cl_session.fits', espec, whdr
; the Gaussian on a fresh load, FWHM field as loaded
spec = spc & fwhm_kernel = float(fwhm_txt)
for i = 0, norders-1 do for j = 0, naps-1 do begin & z = where(finite(spc[*,1,i*naps+j]) eq 1) & mc_convolvespec,findgen(stop-start+1),spc[z,1,i*naps+j],fwhm_kernel,cflux,cerror,ERROR=spc[z,2,i*naps+j] & spec[z,1,i*naps+j] = cflux & spec[z,2,i*naps+j] = cerror & endfor
print, 'GHIST This spectrum has been convolved with a Gaussian of FWHM='+strcompress(fwhm_kernel,/RE)+' pixels.'
writefits, out+'cl_session_g.fits', spec
; MC_CRANGE /KLT: degree 4 with the width 3.6
r = mc_crange(4.,3.6,'SG Degree',/KLT,CANCEL=c4) & r = mc_crange(3.,3.6,'SG Degree',/KLT,CANCEL=c3)
print, 'CRANGE', c4, c3
; ---- SPEX2TEXT (writes into the current directory)
cd, out+'s2t'
spex2text, out+'comb_NGC134.fits' & spex2text, out+'spectra0116.fits' & spex2text, out+'mc_merged.fits'
spex2text, getenv('TSR_ROOT')+'/validation/py_out/HD225200_ci_spec.fits'
; ---- XGETOFFSET_ZOOM: the template and observed home order (order 4), three clicked lines
calinfo = mc_readlinecalfile(pkg+'data/OHXD_LineCal.fits')
whome = reform(calinfo.xcorspec[*,0]) & fhome = reform(calinfo.xcorspec[*,1])
k = mc_nantrim(whome,2) & whome = whome[k] & fhome = fhome[k]
s1 = readfits(out+'wc_spec1d_4.fits',/SILENT) & wobs = reform(s1[*,0]) & fobs = reform(s1[*,1])
d = {x1:whome, x2:wobs, y1:fhome, y2:fobs}
expr = 'p[0] + p[1]*x + gaussian(x,p[2:4])'
clicks = [[784.3,783.6],[1182.2,1181.4],[920.7,920.2]]
for c = 0, 2 do begin & g1 = 1. & g1 = clicks[0,c] & mc_tabinv,d.x1,g1,idx & botidx = 0 > (idx-15) & topidx = (idx+15) < (n_elements(d.x1)-1) & xwin = d.x1[botidx:topidx] & ywin = d.y1[botidx:topidx] & mc_tabinv,xwin,g1,idx & start = float([0.,0.,ywin[round(idx)],g1,2]) & fit1 = mpfitexpr(expr,xwin,ywin,dummy,start,/QUIET,/NOCOVAR) & pos1 = 1. & pos1 = fit1[3] & g2 = 1. & g2 = clicks[1,c] & mc_tabinv,d.x2,g2,idx & botidx = 0 > (idx-15) & topidx = (idx+15) < (n_elements(d.x2)-1) & xwin = d.x2[botidx:topidx] & ywin = d.y2[botidx:topidx] & mc_tabinv,xwin,g2,idx & start = float([0.,0.,ywin[round(idx)],g2,2]) & fit2 = mpfitexpr(expr,xwin,ywin,dummy,start,/QUIET,/NOCOVAR) & pos2 = 1. & pos2 = fit2[3] & offset = pos2-pos1 & print, 'GOFF', c, fit1, fit2, pos1, pos2, offset, format='(a,i2,13e24.16)' & endfor
print, 'DONE'
exit
