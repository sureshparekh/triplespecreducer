; The default reduction as one chain, every stage writing the file the IDL tools write,
; each stage reading only the files of the stage before it:
;   xspextool  (pairs 116-117, 118-119 HD225200; 120-121, 126-127 NGC134), MC_WRITESPEC
;   xcombspec  (order 5 scaled by XMC_SCALESPEC's default Median, Robust Weighted Mean 8),
;              XCOMBSPEC_WRITEFILE
;   xtellcor   (Deconvolution kernel on Pa delta, order 6; shifts on fixed windows),
;              XTELLCOR_WRITEFILE: corrected, telluric and modified-Vega files
;   xmergeorders (anchor 3, then 4, 5, 6, 7, scale 1, weighted), XMERGEORDERS_WRITEFILE
;   xcleanspec (Savitzky-Golay and Gaussian at the panel's defaults), XCLEANSPEC_WRITEFILE
; The case says how the point sources are extracted:
;   'gui'  as the TS4_SOAR panel does by default: PS BG Subtraction switches optimal
;          extraction on (PSFRADIUS 2.2, the SPATCOEFFS map, background 2.5/5/1)
;   'sum'  MC_EXTSPEC2D without PSFRADIUS: the aperture sum with the background, which
;          is what ref_dg.pro does (the TS4 panel cannot make it)
; The extended source is the same in both (A-B, 20.1 arcsec, radius 1; XSPEXTOOL_WRITESPEC
; passes no background regions for an A-B pair). Inputs typed as MC_CFLD types them;
; TRACEsIGTHRESH read with FIX() is 0. The flat and wavecal are validation/idl_out's.
; Batch IDL has no widgets: ref_combtell.pro's stand-ins are used (MC_CFLD reads its fields
; from a common block, XVSPEC does nothing) with its XCOMBSPEC state and scaling helpers.
; Outputs: validation/idl_out/chain_<case>/ and its ref_chain.log.
; Run: cd validation/idl; IDL_PATH="+$SPEXTOOL_DIR:<IDL_DEFAULT>" idl -e "ref_chain,'gui'"

pro ch_extract, files, fl, wc, aw, at, psmode, gui, opaths, appos_xs, flatname, wcname, keywords
  mc_readflat, fl, flat,ncols,nrows,modename,ps,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms,rotation,OMASK=omask
  mc_rdwavecal2d, wc, wavecal, spatcal, wctype, xo2w, wdeg, odeg, homeorder, xy2w, xy2s, wxdeg, wydeg, sxdeg, sydeg, disp, ROTATE=rotation
  kna_readts4fits, files, data, hdrinfo, var, KEYWORD=keywords, ROTATE=rotation, /PAIR, ITOT=itot
  data = data/flat & var = var/flat^2
  if psmode then begin
    prof = mc_mkspatprof2d(data, omask, wavecal, spatcal, slith_arc/slith_pix, 1L, aw, at, 0.4, $
                           SPATCOEFFS=spc, MEDSPATCOEFFS=mspc)
    mc_findpeaks, prof, 2, intarr(norders)+1, positions, apsign, /AUTO
    tracecoeffs = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs, xranges, positions, 5L, 5L, $
                                 3.0, float(fix('0.1')), 1L)
    apsigndat = fltarr(2)
    for i=0,norders-1 do for j=0,1 do begin & profile = prof.(i) & med = median(profile[*,1]) & tabinv,profile[*,0],positions[j,i],idx & apsigndat[j] = (profile[idx,1] gt med) ? 1:-1 & endfor
    apsign = fix(apsigndat) & naps = 2
    apradius = 1.0
    apradii = replicate(apradius,2,norders)
    psbginfo = [2.5,5.0]
    bgorder = 1.0
    if gui then begin
      psfradius = 2.2
      for i = 0, norders-1 do begin
        key = 'Order'+string(i,format='(i2.2)')
        smapinfo = (i eq 0) ? create_struct(key+'Y',spc.(i*2),key+'C',spc.(i*2+1)) : $
                   create_struct(smapinfo,key+'Y',spc.(i*2),key+'C',spc.(i*2+1))
      endfor
    endif
    s = mc_extspec2d(data,var,omask,orders,wavecal,spatcal,ps,apsign,tracecoeffs,apradii, $
                     PSBGINFO=psbginfo,BGDEG=bgorder,PSFRADIUS=psfradius,SMAPINFO=smapinfo)
    psbginfo = [2.5,5.0,bgorder]               ; XSPEXTOOL_WRITESPEC: [bgstart,bgwidth,bgorder]
  endif else begin
    positions = rebin([appos_xs],1,norders)
    tracecoeffs = fltarr(2,norders) & for i=0,norders-1 do tracecoeffs[*,i] = [positions[0,i],0]
    apsign = intarr(1)+1 & naps = 1
    apradii = rebin([1.0],1,norders)
    s = mc_extspec2d(data,var,omask,orders,wavecal,spatcal,ps,apsign,tracecoeffs,apradii,XSBGINFO=[[2.,6.],[32.,36.]],BGDEG=bgorder_xs)
    apradius = [1.0]
  endelse
  ntags = n_tags(s) & sizes = lonarr(2,ntags)
  for i = 0,ntags-1 do sizes[*,i] = size(s.(i),/DIMEN)
  spec = dblarr(max(sizes), max(sizes[1,*]), ntags)*!values.f_nan
  for i = 0,ntags-1 do spec[0:(sizes[0,i]-1),*,i] = s.(i)
  junk = tag_exist(hdrinfo[0].vals,'FILENAME',INDEX=idx)
  for side = 0, 1 do begin
    want = (side eq 0) ? 1:-1
    pos = where(apsign eq want,pnaps)
    if pnaps eq 0 then continue
    hdr0 = hdrinfo[side] & aimage = hdr0.vals.(idx) & sky = hdrinfo[1-side].vals.(idx)
    hdr0.vals.(idx) = file_basename(opaths[side])+'.fits'
    a = replicas(apsign,norders) & z = where(a eq want)
    ; XSPEXTOOL_WRITESPEC, A-B: PSBGINFO for a point source only, no XSBGINFO
    mc_writespec,spec[*,*,z],opaths[side]+'.fits',aimage,itot[side],sky,flatname,pnaps,orders,hdr0,positions[pos,*],apradius,modename,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,'um','DN / s','!7k!5 (!7l!5m)','f (!5DN s!u-1!N)','TS4_SOAR','v4.1',PSFRADIUS=psfradius,PSBGINFO=psbginfo,WAVECAL=wcname,WAVETYPE='Vacuum',DISP=disp,LINEARITY=0,LINCORMAX=80000L
  endfor
end

pro ref_chain, case_
  forward_function ct_xcstate
  common ct_fld, fv
  sp = getenv('SPEXTOOL_DIR')
  raw  = getenv('TSR_ROOT')+'/data/raw/2025-11-02/'
  idl_out = getenv('TSR_ROOT')+'/validation/idl_out/'
  out = idl_out+'chain_'+case_+'/'
  file_mkdir, out
  !quiet = 1
  set_plot, 'Z'
  gui = case_ eq 'gui'
  resolve_routine, 'ref_combtell', /COMPILE_FULL_FILE    ; MC_CFLD, XVSPEC, CT_XCSTATE ...
  resolve_routine, 'xcombspec', /COMPILE_FULL_FILE
  resolve_routine, 'xmc_scalespec', /COMPILE_FULL_FILE
  resolve_routine, 'xmc_maskspec', /COMPILE_FULL_FILE
  resolve_routine, 'xtellcor', /COMPILE_FULL_FILE
  resolve_routine, 'xmergeorders', /COMPILE_FULL_FILE
  resolve_routine, 'xcleanspec', /COMPILE_FULL_FILE
  keywords = ['FILENAME','OBJECT','UT','DATE','EXPTIME','RA','DEC','FSAMPLE','COADDS','AIRMASS']
  flatname = 'flat.fits' & wcname = 'wavecal_NGC134.fits'
  fl = idl_out+flatname & wc = idl_out+wcname
  atr = readfits(sp+'/data/atran4000.fits',/SILENT) & aw = reform(atr[*,0]) & at = reform(atr[*,1])
  print, 'CASE ', case_
  ; ---- xspextool
  ch_extract, raw+['SPEC_SO2025D-006-0116.fits','SPEC_SO2025D-006-0117.fits'], fl, wc, aw, at, 1, gui, out+['spectra0116','spectra0117'], 0, flatname, wcname, keywords
  ch_extract, raw+['SPEC_SO2025D-006-0118.fits','SPEC_SO2025D-006-0119.fits'], fl, wc, aw, at, 1, gui, out+['spectra0118','spectra0119'], 0, flatname, wcname, keywords
  ch_extract, raw+['SPEC_SO2025D-006-0120.fits','SPEC_SO2025D-006-0121.fits'], fl, wc, aw, at, 0, gui, out+['spectra0120','spectra0121'], 20.1, flatname, wcname, keywords
  ch_extract, raw+['SPEC_SO2025D-006-0126.fits','SPEC_SO2025D-006-0127.fits'], fl, wc, aw, at, 0, gui, out+['spectra0126','spectra0127'], 20.1, flatname, wcname, keywords
  print, 'EXTRACTED'
  ; ---- xcombspec
  fv = strarr(16) & fv[0] = '8.0'                         ; Robust threshold
  st = ct_xcstate(out+['spectra0116','spectra0117','spectra0118','spectra0119']+'.fits','No',0)
  st.r.combinestat = 0 & st.r.scaleorderidx = 2           ; Robust Weighted Mean; order 5
  ct_scalespec, st, 'Median', [!values.f_nan,!values.f_nan], 0
  ct_write, st, out, 'comb_HD225200'
  st = ct_xcstate(out+['spectra0120','spectra0126']+'.fits','No',0)
  st.r.combinestat = 0 & st.r.scaleorderidx = 2
  ct_scalespec, st, 'Median', [!values.f_nan,!values.f_nan], 0
  ct_write, st, out, 'comb_NGC134'
  ; ---- xtellcor: the kernel, the telluric spectra and the shifts as ref_dg.pro makes them
  mc_readspec,out+'comb_HD225200.fits',std,stdhdr,stdobsmode,start,stop,stdnorders,stdnaps,stdorders,stdxunits,stdyunits,slith_pix,slith_arc,slitw_pix,slitw_arc,stdrp,stdairmass,xtitle,ytitle,/SILENT
  std = mc_caaspecflags(std)
  mc_readspec,out+'comb_NGC134.fits',obj,objhdr,objobsmode,start,stop,objnorders,objnaps,objorders,objxunits,objyunits,slith_pix,slith_arc,slitw_pix,slitw_arc,objrp,objairmass,/SILENT
  obj = mc_caaspecflags(obj)
  stddisp = mc_fxpar(stdhdr,'DISP*')
  restore, sp+'/data/lvega5.sav'
  wvega = wvin/10000. & cf2vega = fc2vin & cfvega = fcvin & fvega = fvin
  bmag = 6.33 & vmag = 6.31
  idx6 = where(stdorders eq 6)
  i0 = mc_nantrim(std[*,0,idx6],2)
  swave = (reform(std[*,0,idx6]))[i0] & sflux = (reform(std[*,1,idx6]))[i0]
  normreg = [[0.990,0.997],[1.013,1.020]]
  goodbad = intarr(n_elements(swave))
  for i=0,1 do begin & tabinv,swave,normreg[*,i],ii & ii = ii[sort(ii)] & goodbad[round(ii[0]):round(ii[1])] = 1 & endfor
  ncoeffs = robustpoly1d(swave,sflux,3,3,0.01,/GAUSSJ,IGOODBAD=goodbad,/SILENT)
  nflux = float(sflux/poly(swave,ncoeffs))
  linereg = [0.997,1.013]
  z = where(swave gt linereg[0] and swave lt linereg[1])
  mc_vegacorr,swave[z],nflux[z],wvega,fvega,cfvega,vshift,redshift
  print, 'VSHIFT', vshift, format='(a,e24.16)'
  wv2 = wvega*(1.0+redshift) & fnvega = fvega/cfvega
  mc_vegaconv,swave,nflux,wv2,cfvega,fnvega,linereg[0],linereg[1],wline,kernel,scale,maxdev,rmsdev,TWF=10.
  print, 'KERNEL', n_elements(kernel), scale, maxdev, rmsdev, wline, format='(a,i6,4e24.16)'
  writefits, out+'tc_kernel.fits', kernel
  vdisp = 1.49774e-5
  ndat = n_elements(kernel) & kx_vegapix = findgen(ndat)
  result = gaussfit(kx_vegapix,kernel,a,NTERMS=3)
  kx_spexpix = (kx_vegapix - a[1])*vdisp/total(stddisp[idx6])
  for i=0,stdnorders-1 do begin & kx_spexwave = kx_spexpix*total(stddisp[i]) & del = max(kx_spexwave,MIN=min)-min & npix = del/vdisp & if not npix mod 2 then npix = npix + 1 & kx_vegawave = (findgen(npix)-npix/2.)*vdisp & linterp,kx_spexwave,kernel,kx_vegawave,nkern & nkern = nkern/total(nkern) & writefits, out+'tc_kern_'+strtrim(stdorders[i],2)+'.fits', nkern & endfor
  npix = n_elements(std[*,0,0])
  match, stdorders, objorders, stdidx
  tellspec = std[*,*,stdidx] & vegaspec = std[*,*,stdidx] & vegaspec[*,2,*] = 1.0
  for i = 0, objnorders-1 do begin
     zz = where(stdorders eq objorders[i])
     nkern = readfits(out+'tc_kern_'+strtrim(stdorders[zz],2)+'.fits',/SILENT)
     mc_mktellspec,std[*,0,zz],std[*,1,zz],std[*,2,zz],vmag,(bmag-vmag),nkern,replicate(scale,npix),wvega,fvega,cfvega,cf2vega,vshift,tellcor,tellcor_error,scvega
     tellspec[*,1,i] = tellcor & tellspec[*,2,i] = tellcor_error & vegaspec[*,1,i] = scvega
  endfor
  sr = [[2.00,2.06],[1.43,1.48],[1.11,1.15],[0.93,0.97],[0.94,1.00]]
  shifts = fltarr(objnorders)
  for i=0,objnorders-1 do begin
    ow = obj[*,0,i] & oflx = obj[*,1,i] & t = mc_nantrim(ow,2) & ow = ow[t] & oflx = oflx[t]
    tw = tellspec[*,0,i] & tf = tellspec[*,1,i] & t = mc_nantrim(tw,2) & tw = tw[t] & tf = tf[t]
    mc_interpspec,tw,tf,ow,nflx
    smooth = median(oflx,4) & oflx = oflx/max(smooth,/NAN)
    z = where(ow gt sr[0,i] and ow lt sr[1,i])
    sh = findgen(151)/50.+(0.-1.5) & rmsv = fltarr(151) & x2 = findgen(n_elements(nflx))
    for k=0,150 do begin & mc_interpspec,x2+sh[k],nflx,x2,y,e & d = oflx*y & d = d/median(d,/EVEN) & mc_moments,d[z],m,v,sdv,/SILENT & rmsv[k] = sdv & endfor
    mn = min(rmsv,mi)
    if mi lt 5 or mi gt 145 then shifts[i] = 0. else begin
      c = poly_fit1d(sh[(mi-5):(mi+5)],rmsv[(mi-5):(mi+5)],2,/SILENT)
      shifts[i] = -c[1]/2./c[2]
    endelse
  endfor
  print, 'SHIFTS', shifts, format='(a,5e24.16)'
  ; XTELLCOR_WRITEFILE (state as ref_combtell.pro builds it; vshift, the deviations and
  ; the airmass are FLOAT fields of XTELLCOR's state)
  vshift = float(vshift) & maxdev = float(maxdev) & rmsdev = float(rmsdev)
  stdairmass = float(stdairmass)
  w = {objoname_fld:[0L,9L], stdspectra_fld:[0L,10L], xtellcor_base:0L}
  rr = {bmag:bmag, continue:4L, maxdev:maxdev, method:'Deconvolution', $
        nytitle:'!5f!D!7k!N!5 (ergs s!E-1!N cm!E-2!N '+'!5!sA!r!u!9 %!5!n'+'!E-1!N', $
        objpath:out, propstdunc:1L, rmsdev:rmsdev, shift:ptr_new(rebin(shifts,objnorders,objnaps)), $
        stdairmass:stdairmass, telluricoutput:1, units:'ergs s-1 cm-2 A-1', vegaoutput:1, $
        vmag:vmag, vshift:vshift}
  d = {dairmass:0.0, objhdr:ptr_new(objhdr), objnaps:objnaps, objorders:ptr_new(objorders), $
       objspec:ptr_new(obj), stdhdr:ptr_new(stdhdr), stdorders:ptr_new(stdorders), $
       tellspec:ptr_new(tellspec), vegaspec:ptr_new(vegaspec)}
  d.dairmass = stdairmass-objairmass
  xs = {w:w, r:rr, d:d}
  fv[9] = 'tc_NGC134' & fv[10] = 'comb_HD225200.fits'
  xtellcor_writefile, xs
  print, 'TELLCOR'
  ; ---- xmergeorders: XMERGEORDERS_LOADSPEC reads the corrected file
  mc_readspec,out+'tc_NGC134.fits',cs,chdr,obsmode,start,stop,norders,naps,corders,xunits,yunits,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,airmass,xtitle,ytitle,/SILENT
  cs = mc_caaspecflags(cs)
  i3 = where(finite(cs[*,1,0]) eq 1)
  merged = reform(cs[i3,*,0])
  for o=4,7 do begin & k = where(corders eq o) & nn = where(finite(cs[*,1,k]) eq 1) & g = reform(cs[nn,*,k]) & mc_mergespec,merged[*,0],double(merged[*,1]),g[*,0],double(g[*,1]),owave,oflux,E1=double(merged[*,2]),E2=double(g[*,2]),BF1=merged[*,3],BF2=g[*,3],OERROR=oerror,OBITFLAG=obitflag & merged = [[owave],[oflux],[oerror],[obitflag]] & endfor
  fv[11] = 'tc_NGC134.fits' & fv[12] = 'merged_NGC134'
  ms = {w:{ispectrum_fld:[0L,11L], oname_fld:[0L,12L]}, $
        r:{weighted:1, ap:0, mergedorders:ptr_new([3,4,5,6,7]), orders:ptr_new(corders), $
           specscales:ptr_new(fltarr(norders)+1.0), path:out}, $
        d:{mergedspec:ptr_new(merged), hdr:ptr_new(chdr)}}
  xmergeorders_writefile, ms
  print, 'MERGED'
  ; ---- xcleanspec: XCLEANSPEC_LOADSPEC, the default fields, XCLEANSPEC_SMOOTHSPEC's loop
  mc_readspec,out+'merged_NGC134.fits',spc,mhdr,obsmode,start,stop,norders,naps,morders,xunits,yunits,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,airmass,xtitle,ytitle,/SILENT
  spc = mc_caaspecflags(spc)
  sw = float(slitw_pix) & sc = (slitw_pix gt 2) ? 1.5:2.5
  sgwin_txt = strtrim(string(slitw_pix*sc,FORMAT='(f4.1)'),2) & fwhm_txt = strtrim(sw,2)
  spec = spc & sgwin = float(sgwin_txt) & sgdeg = float('2')
  for i = 0, norders-1 do for j = 0, naps-1 do begin & z = where(finite(spc[*,1,i*naps+j]) eq 1) & cflux = spc[z,1,i*naps+j] & cerror = spc[z,2,i*naps+j] & cflux = savitzky_golay(cflux,sgwin,DEGREE=sgdeg,ERROR=cerror) & spec[z,1,i*naps+j] = cflux & spec[z,2,i*naps+j] = cerror & endfor
  hist = 'This spectrum has been convolved with a Savitzky-Golay filter of width '+strtrim(sgwin,2)+' pixels and degree '+strtrim(sgdeg,2)+'.'
  fv[13] = 'clean_sg_NGC134'
  cl = {w:{oname_fld:[0L,13L]}, r:{smthhistory:hist, path:out}, d:{hdr:ptr_new(mhdr), espec:ptr_new(spec)}}
  xcleanspec_writefile, cl
  spec = spc & fwhm_kernel = float(fwhm_txt)
  for i = 0, norders-1 do for j = 0, naps-1 do begin & z = where(finite(spc[*,1,i*naps+j]) eq 1) & mc_convolvespec,findgen(stop-start+1),spc[z,1,i*naps+j],fwhm_kernel,cflux,cerror,ERROR=spc[z,2,i*naps+j] & spec[z,1,i*naps+j] = cflux & spec[z,2,i*naps+j] = cerror & endfor
  hist = 'This spectrum has been convolved with a Gaussian of FWHM='+strcompress(fwhm_kernel,/RE)+' pixels.'
  fv[14] = 'clean_g_NGC134'
  cl = {w:{oname_fld:[0L,14L]}, r:{smthhistory:hist, path:out}, d:{hdr:ptr_new(mhdr), espec:ptr_new(spec)}}
  xcleanspec_writefile, cl
  print, 'DONE'
end
