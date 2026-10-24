; Reference for the options beyond the default path: combined images (xspextool_combineimgs),
; A-mode extraction of a combined image, xtellcor units, edited H-line scales, combine apertures.
pro ref_extras
  sp = getenv('SPEXTOOL_DIR')
  raw  = getenv('TSR_ROOT')+'/data/raw/2025-11-02/'
  out  = getenv('TSR_ROOT')+'/validation/idl_out/'
  keywords = ['FILENAME','OBJECT','UT','DATE','EXPTIME','RA','DEC','FSAMPLE','COADDS','AIRMASS']
  ; ---- combine images: A-B, ABBA->ABAB, Robust Mean (Std Error), thresh 8, no scaling/bg
  files = raw+'SPEC_SO2025D-006-0'+['116','117','118','119']+'.fits'
  idx = intarr(4) & for i=0,1 do begin & if i mod 2 then begin & idx[i*2] = i*2+1 & idx[i*2+1] = i*2 & endif else begin & idx[i*2] = i*2 & idx[i*2+1] = i*2+1 & endelse & endfor
  files = files[idx]
  kna_readts4fits,files,idata,hdrinfo,ivar,KEYWORD=keywords,NIMAGES=nimages,/PAIR,ROTATE=0,BITINFO={lincormax:80000.0,lincormaxbit:0},BITMASK=bitmask,ITOT=itot
  mc_meancomb,idata,mean,mvar,ROBUST=8.0,OGOODBAD=ogoodbad
  writefits, out+'ci_mean.fits', reform(mean,1024,2048) & writefits, out+'ci_var.fits', reform(mvar,1024,2048)
  print, 'ITOT', total(itot)/2.
  ; ---- the same with the sky subtracted (mc_subspecsky, ybuffer 3), images read rotated
  mc_readflat, out+'flat.fits', flat,ncols,nrows,modename,ps,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms,rotation,OMASK=omask
  ; (its own variable names: the A-mode extraction below uses MEAN and MVAR)
  kna_readts4fits,files,bdata,bhdrinfo,bvar,KEYWORD=keywords,NIMAGES=bnimages,/PAIR,ROTATE=rotation,BITINFO={lincormax:80000.0,lincormaxbit:0},BITMASK=bbitmask,ITOT=bitot
  mc_subspecsky,bdata,bvar,edgecoeffs,norders,3,xranges
  mc_meancomb,bdata,bmean,bmvar,ROBUST=8.0,OGOODBAD=bogoodbad
  writefits, out+'ci_bg_mean.fits', bmean & writefits, out+'ci_bg_var.fits', bmvar
  ; ---- A-mode extraction of the combined image (as xspextool reads it from proc/)
  mc_readflat, out+'flat.fits', flat,ncols,nrows,modename,ps,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms,rotation,OMASK=omask
  mc_rdwavecal2d, out+'wavecal_NGC134.fits', wavecal, spatcal, wctype, xo2w, wdeg, odeg, homeorder, xy2w, xy2s, wxdeg, wydeg, sxdeg, sydeg, disp, ROTATE=rotation
  atr = readfits(sp+'/data/atran4000.fits',/SILENT) & aw = reform(atr[*,0]) & at = reform(atr[*,1])
  ; KNA_READTS4FITS with ROTATE=0 stores each frame into a [2048,1024] array in raw memory order:
  ; the values are right, the shape is not. Undo it to extract.
  data = rotate(reform(mean,1024,2048),rotation)/flat & var = rotate(reform(mvar,1024,2048),rotation)/flat^2
  prof = mc_mkspatprof2d(data, omask, wavecal, spatcal, slith_arc/slith_pix, 1, aw, at, 0.4)
  mc_findpeaks, prof, 2, intarr(norders)+1, positions, apsign, /AUTO
  tracecoeffs = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs, xranges, positions, 5, 5, 3.0, 0.1, 1)
  apsigndat = fltarr(2)
  for i=0,norders-1 do for j=0,1 do begin & profile = prof.(i) & med = median(profile[*,1]) & tabinv,profile[*,0],positions[j,i],ii & apsigndat[j] = (profile[ii,1] gt med) ? 1:-1 & endfor
  print, 'CI_POS', positions & print, 'CI_SIGN', fix(apsigndat)
  s = mc_extspec2d(data,var,omask,orders,wavecal,spatcal,ps,fix(apsigndat),tracecoeffs,replicate(1.0,2,norders),PSBGINFO=[2.5,5.0],BGDEG=1)
  for k=0,n_tags(s)-1 do writefits, out+'ci_spec_'+strtrim(k,2)+'.fits', s.(k)
  ; ---- units of the telluric spectrum
  t = readfits(out+'tc_tellspec.fits',/SILENT)
  jy = mc_chfunits(t[*,0,1],t[*,1,1],0,1,4,IERROR=t[*,2,1],OERROR=ejy)
  wm = mc_chfunits(t[*,0,1],t[*,1,1],0,1,0)
  writefits, out+'un_jy.fits', [[jy],[ejy]] & writefits, out+'un_wm.fits', wm
  ; B-V is (bmag-vmag) in FLOAT, as xtellcor passes it
  ; ---- edited H-line scales on order 4 (Br 11 at 1.68111 um pushed to 0.8 x scale)
  mc_readspec,out+'comb_HD225200.fits',std,stdhdr,stdobsmode,start,stop,stdnorders,stdnaps,stdorders,stdxunits,stdyunits,slith_pix,slith_arc,slitw_pix,slitw_arc,stdrp,stdairmass,xtitle,ytitle,/SILENT
  std = mc_caaspecflags(std)
  readcol,sp+'/data/HI.dat',hlines,hnames,FORMAT='D,A',COMMENT='#',DELIMITER='|',/SILENT
  restore, sp+'/data/lvega5.sav' & wvega = wvin/10000.
  scale = 0.74261242 & vshift = 33.111367
  i4 = 1
  mn = min(std[*,0,i4],MAX=mx,/NAN)
  z = where(hlines gt mn and hlines lt mx,count)
  arr = [mn,scale] & for j=0,count-1 do arr = [[arr],[hlines[z[j]],scale]] & arr = [[arr],[mx,scale]]
  srt = sort(arr[0,*]) & arr = arr[*,srt]
  k11 = where(abs(arr[0,*]-1.6811) lt 0.001) & arr[1,k11] = 0.8*scale
  writefits, out+'hl_cpoints.fits', arr
  scales = float(spline(arr[0,*],arr[1,*],std[*,0,i4],10.))   ; xmc_scalelines returns an FLTARR
  writefits, out+'hl_scales.fits', scales
  nkern = readfits(out+'tc_kern_4.fits',/SILENT)
  mc_mktellspec,std[*,0,i4],std[*,1,i4],std[*,2,i4],6.31,(6.33-6.31),nkern,scales,wvega,fvin,fcvin,fc2vin,vshift,tellcor,tellcor_error,scvega
  writefits, out+'hl_tell4.fits', [[tellcor],[tellcor_error]]
  ; ---- xmc_scalelines_estimatescale on order 4, region around Br 11, from the
  ;      untouched control points (as the panel starts)
  idx = mc_nantrim(std[*,0,i4],2)
  wave = reform(std[idx,0,i4]) & flux = reform(std[idx,1,i4]) & err = reform(std[idx,2,i4])
  arr[1,*] = scale
  scalespec = spline(arr[0,*],arr[1,*],wave,10.)
  mc_mktellspec,wave,flux,err,6.31,(6.33-6.31),nkern,scalespec,wvega,fvin,fcvin,fc2vin,vshift,tc
  ereg = [1.6770,1.6850]
  tel_flux = 1./tc
  z = where(wave gt ereg[0] and wave lt ereg[1])
  ii = where(hlines gt ereg[0] and hlines lt ereg[1],count)
  wguess = hlines[ii]
  tabinv,wave,wguess,ix & tabinv,wave,ereg[0],idlo & tabinv,wave,ereg[1],idhi
  fguess = tel_flux[ix] - 0.5*(tel_flux[idlo] + tel_flux[idhi])
  sguess = (ereg[1]-ereg[0])/5.
  cguess = poly_fit1d(wave[z],tel_flux[z],1,/SILENT)
  mguess = (tel_flux[idlo] - tel_flux[idhi])/(wave[idlo]-wave[idhi])
  result = mpfitpeak(wave[z],tel_flux[z],gcoeff,NTERMS=5,/LORENTZIAN,ESTIMATES=[fguess,wguess,sguess,cguess,mguess])
  wv = wvega*(1+vshift/2.99792458E5)
  wmin = min(wave,/NAN,MAX=wmax) & zv = where(wv gt wmin and wv lt wmax,count)
  nadd = round(n_elements(nkern)/2.) & vidx = findgen(count+2*nadd)+(zv[0]-nadd)
  magscale = 10.0^(-0.4*(6.31-0.03)) & ebmv = ((6.33-6.31) - 0.00) > 0.0 & avred = 3.10*ebmv & redfact = 10.0^(0.4*avred)
  magfact = replicate(redfact*magscale,n_elements(wave))
  nfvconv = convol((fvin[vidx]/fcvin[vidx]-1.0),nkern)
  linterp,wv[vidx],nfvconv,wave,rnfvconv
  linterp,wv,fcvin,wave,rfcvega
  mc_redden,wave,rfcvega,ebmv,rrfcvega
  s = (poly_smooth(flux,5)/(magfact*poly(wave,gcoeff[3:4])*rrfcvega))- 1.0
  s = temporary(s)/rnfvconv
  zc = where(arr[0,*] gt ereg[0] and arr[0,*] lt ereg[1])
  tabinv,wave,arr[0,zc],ic
  print, 'ESTSCALE', s[ic], gcoeff, format='(a,6e24.16)'
  writefits, out+'hl_est_gcoeff.fits', gcoeff
  writefits, out+'hl_est_tc.fits', tc
  writefits, out+'hl_est_scalespec.fits', scalespec
  writefits, out+'hl_est_est.fits', [fguess,wguess,sguess,cguess,mguess]
  ; ---- xmergeorders with an auto scale, a green trim, a gap and a blue trim:
  ;      anchor 4; add 5 (trim Green 'Right': kept below 1.47, auto scale over 1.42-1.47);
  ;      add 3 (no overlap); add 6 (trim Blue 'Right': merged kept below 1.30)
  cs = readfits(out+'tc_NGC134.fits',/SILENT)
  cs = mc_caaspecflags(cs)
  nn = where(finite(cs[*,1,1]) eq 1) & merged = reform(cs[nn,*,1])
  nn = where(finite(cs[*,1,2]) eq 1) & add = reform(cs[nn,*,2])
  trim = 1.47 & z = where(add[*,0] lt trim) & green = add[z,*] & blue = merged
  xr = [1.42,1.47]
  zb = where(blue[*,0] gt xr[0] and blue[*,0] lt xr[1]) & zy = where(green[*,0] gt xr[0] and green[*,0] lt xr[1])
  mc_interpspec,reform(green[zy,0]),reform(green[zy,1]),reform(blue[zb,0]),ngf,nge,IYERROR=reform(green[zy,2])
  bf = reform(blue[zb,1]) & be = reform(blue[zb,2])
  zz = where(finite(bf) eq 1 and finite(ngf) eq 1 and finite(be) eq 1 and finite(nge) eq 1)
  denom = be^2 + nge^2
  scale = 1. & scale[0] = total(ngf[zz]*bf[zz]/denom[zz])/total(ngf[zz]^2/denom[zz])
  print, 'MSCALE', scale, format='(a,e24.16)'
  green = add[z,*] & green[*,1] = green[*,1]*scale & green[*,2] = green[*,2]*abs(scale)
  mc_mergespec,blue[*,0],double(blue[*,1]),green[*,0],double(green[*,1]),owave,oflux,E1=double(blue[*,2]),E2=double(green[*,2]),BF1=blue[*,3],BF2=green[*,3],OERROR=oerror,OBITFLAG=obitflag
  merged = [[owave],[oflux],[oerror],[obitflag]]
  nn = where(finite(cs[*,1,0]) eq 1) & green = reform(cs[nn,*,0]) & blue = merged
  mc_mergespec,blue[*,0],double(blue[*,1]),green[*,0],double(green[*,1]),owave,oflux,E1=double(blue[*,2]),E2=double(green[*,2]),BF1=blue[*,3],BF2=green[*,3],OERROR=oerror,OBITFLAG=obitflag
  merged = [[owave],[oflux],[oerror],[obitflag]]
  nn = where(finite(cs[*,1,3]) eq 1) & green = reform(cs[nn,*,3])
  trim = 1.30 & z = where(merged[*,0] lt trim or finite(merged[*,0]) eq 0) & blue = merged[z,*]
  mc_mergespec,blue[*,0],double(blue[*,1]),green[*,0],double(green[*,1]),owave,oflux,E1=double(blue[*,2]),E2=double(green[*,2]),BF1=blue[*,3],BF2=green[*,3],OERROR=oerror,OBITFLAG=obitflag
  merged = [[owave],[oflux],[oerror],[obitflag]]
  writefits, out+'m_trims.fits', merged
  ; ---- xcleanspec: fix 1.27-1.29 (flux view), save, remove 1.82-1.80, save; cursor positions FLOAT
  spc = readfits(out+'m_NGC134.fits',/SILENT)
  espec = spc & pspec = spc
  xrange = [1.27,1.29]
  mc_tabinv,espec[*,0],xrange,idx & idx = round(idx)
  x = [espec[idx[0],0],espec[idx[1],0]] & y = [espec[idx[0],1],espec[idx[1],1]] & e = [espec[idx[0],2],espec[idx[1],2]]
  coeff = poly_fit1d(x,y,1,/SILENT) & pspec[idx[0]:idx[1],1] = poly(espec[idx[0]:idx[1],0],coeff)
  coeff = poly_fit1d(x,e,1,/SILENT) & pspec[idx[0]:idx[1],2] = poly(espec[idx[0]:idx[1],0],coeff)
  pspec[idx[0]:idx[1],3] = pspec[idx[0]:idx[1],3] + 4
  espec = pspec
  reg = [1.82,1.80]
  mc_tabinv,espec[*,0],reg,idx & idx = round(idx) & mn = idx[0] < idx[1] & mx = idx[0] > idx[1]
  n = n_elements(pspec[*,0]) & pspec[(mn > 0):(mx < (n-1)),1] = !values.f_nan & pspec[(mn > 0):(mx < (n-1)),2] = !values.f_nan
  espec = pspec
  writefits, out+'cl_edit.fits', espec
  ; ---- xcombspec statistics 1-8 on the four HD225200 spectra, order 5 (index 2), as loaded
  ;      by xcombspec (first spectrum's grid, FLTARR stack), unscaled
  sf = out+['spectra0116','spectra0117','spectra0118','spectra0119']+'.fits'
  mc_readspec,sf[0],first,hdr,obsmode,start,stop,norders,naps,orders,xunits,yunits,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,airmass,xtitle,ytitle,instr,/SILENT
  npix = n_elements(first[*,0,0]) & stk = fltarr(npix,4,4)*!values.f_nan
  jj = 2*naps
  for i=0,3 do begin & d = readfits(sf[i],/SILENT) & d = mc_caaspecflags(d) & z = where(finite(first[*,0,jj]) eq 1) & mc_interpspec,d[z,0,jj],d[z,1,jj],first[z,0,jj],nf,ne_,IYERROR=d[z,2,jj],LEAVENANS=1 & stk[z,0,i] = first[z,0,jj] & stk[z,1,i] = nf & stk[z,2,i] = ne_ & endfor
  spec = reform(stk[*,1,*]) & err = reform(stk[*,2,*]) & mask = intarr(npix,4)+1
  res = fltarr(npix,2,9)
  mc_meancomb,spec,m,v,/RMS,MASK=mask,ROBUST=8.0,OGOODBAD=og,/SILENT & res[*,0,1] = m & res[*,1,1] = sqrt(v)
  mc_meancomb,spec,m,v,MASK=mask,ROBUST=8.0,OGOODBAD=og,/SILENT & res[*,0,2] = m & res[*,1,2] = sqrt(v)
  mc_meancomb,spec,m,v,DATAVAR=err^2,MASK=mask,/SILENT & res[*,0,3] = m & res[*,1,3] = sqrt(v)
  mc_meancomb,spec,m,v,/RMS,MASK=mask,/SILENT & res[*,0,4] = m & res[*,1,4] = sqrt(v)
  mc_meancomb,spec,m,v,MASK=mask,/SILENT & res[*,0,5] = m & res[*,1,5] = sqrt(v)
  mc_medcomb,spec,m,v,/MAD,MASK=mask & res[*,0,6] = m & res[*,1,6] = sqrt(v)
  mc_medcomb,spec,m,v,MASK=mask & res[*,0,7] = m & res[*,1,7] = sqrt(v)
  m = total(spec*mask,2,/NAN) & v = total(mask*double(err)^2,2,/NAN) & res[*,0,8] = m & res[*,1,8] = sqrt(v)
  writefits, out+'cs_stats.fits', res
  ; ---- mc_findpeaks GUESS and FIXED on the point-source profiles
  ; (its own copies: blocks above reuse the flat's variable names)
  mc_readflat, out+'flat.fits', fpflat,fncols,fnrows,fmode,fps,fslith_pix,fslith_arc,fslitw_pix,fslitw_arc,frp,fnorders,forders,fedge,fxr,frms,frot,OMASK=fomask
  mc_rdwavecal2d, out+'wavecal_NGC134.fits', fwavecal, fspatcal, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12, ROTATE=frot
  prof = mc_mkspatprof2d(readfits(out+'ps_work.fits',/SILENT), fomask, fwavecal, fspatcal, fslith_arc/fslith_pix, 1, aw, at, 0.4)
  g = fltarr(2,fnorders) & g[0,*] = 9.5 & g[1,*] = 28.3
  mc_findpeaks, prof, 2, intarr(fnorders)+1, pg, sg, GUESS=g
  mc_findpeaks, prof, 2, intarr(fnorders)+1, pf, sf_, FIXED=g
  writefits, out+'fp_guess.fits', pg & writefits, out+'fp_fixed.fits', pf
  writefits, out+'fp_prof3.fits', prof.(0)
  print, 'FPSIGN', sg, sf_
  ; ---- A-Sky: raw 120 minus a combined image read as xspextool reads it (MRDFITS, rotated)
  skyp = getenv('TSR_ROOT')+'/validation/py_out/HD225200_combined.fits'
  kna_readts4fits,raw+'SPEC_SO2025D-006-0120.fits',adata,ahdr,avar,KEYWORD=keywords,ROTATE=frot,BITINFO={lincormax:80000.0,lincormaxbit:0},BITMASK=abm,ITOT=aitot
  sky = rotate(mrdfits(skyp,1,/SILENT),frot) & skyvar = rotate(mrdfits(skyp,2,/SILENT),frot)
  adata = temporary(adata)-sky & avar = temporary(avar)+skyvar
  adata = adata/fpflat & avar = avar/fpflat^2
  writefits, out+'as_work.fits', adata
  s = mc_extspec2d(adata,avar,fomask,forders,fwavecal,fspatcal,fps,[1],[[20.1,0],[20.1,0],[20.1,0],[20.1,0],[20.1,0]],[[1.0],[1.0],[1.0],[1.0],[1.0]])
  for k=0,n_tags(s)-1 do writefits, out+'as_spec_'+strtrim(k,2)+'.fits', s.(k)
  print, 'DONE'
end
