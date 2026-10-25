; Point-source extraction as xspextool's Do All Steps makes it with the TS4_SOAR
; defaults: PSBGSUB=1 switches optimal extraction on (the 'PS BG Subtraction'
; event sets state.r.optext[0]), so XSPEXTOOL_PSEXTRACTSPEC passes PSFRADIUS
; (PSPSFRAD 2.2) and SMAPINFO (the SPATCOEFFS of MC_MKSPATPROF2D, or
; MEDSPATCOEFFS with 'Median Profile') to MC_EXTSPEC2D. Every user input is
; typed as MC_CFLD types the widget text: radii and background as FLOAT,
; trace step/sum/degree as LONG, TRACEsIGTHRESH read by FIX() as 0.
; Cases (prefix of the output files):
;   op  Auto, traced, optimal (the GUI default)    om  the same, Median Profile
;   nb  PS BG Subtraction off: sum, no background
;   gs  Guess 9.5,28.3, traced, optimal            fx  Fix 9.5,28.3 (Use Ap Positions)
;   it  Import Trace of op's trace file            os  PS Orders 3,5,6 only (Auto)
;   am  A mode, frame 116, Fix 9.8,20.0         (fs: Fix with PS Orders 3,5,6 -- IDL hangs)
; Run: cd validation/idl; IDL_PATH="+$SPEXTOOL_DIR:<IDL_DEFAULT>" idl -e ref_optimal
pro op_smap, coeffs, doorders, smapinfo
  ; XSPEXTOOL_PSEXTRACTSPEC: the selected orders' coefficient pairs
  idx = 0
  for i = 0, n_elements(doorders)-1 do begin
     if doorders[i] ne 1 then continue
     key = 'Order'+string(i,format='(i2.2)')
     if idx eq 0 then smapinfo = create_struct(key+'Y',coeffs.(i*2),key+'C',coeffs.(i*2+1)) $
     else smapinfo = create_struct(smapinfo,key+'Y',coeffs.(i*2),key+'C',coeffs.(i*2+1))
     idx = idx + 1
  endfor
end

function op_signs, prof, appos, doorders, psnaps
  ; XSPEXTOOL_PSDEFINEAP
  apsigndat = fltarr(psnaps)
  for i = 0, n_elements(doorders)-1 do begin
     if doorders[i] ne 1 then continue
     profile = prof.(i)
     for j = 0, psnaps-1 do begin
        med = median(profile[*,1])
        tabinv, profile[*,0], appos[j,i], idx
        apsigndat[j] = (profile[idx,1] gt med) ? 1:-1
     endfor
  endfor
  return, fix(apsigndat)
end

pro op_extract, tag, data, var, omask, orders, wavecal, spatcal, ps, apsign, tracecoeffs, $
                doorders, psnaps, coeffs, out, OPT=opt, BG=bg
  z = where(doorders eq 1, norders)
  apradius = 1.0
  apradii = replicate(apradius, psnaps, norders)
  if keyword_set(BG) then begin
     psbginfo = [2.5,5.0]
     bgorder = 1.0
  endif
  if keyword_set(OPT) then begin
     psfradius = 2.2
     op_smap, coeffs, doorders, smapinfo
  endif
  t0 = systime(1)
  s = mc_extspec2d(data,var,omask,orders[z],wavecal,spatcal,ps,apsign,tracecoeffs,apradii, $
                   PSBGINFO=psbginfo,BGDEG=bgorder,PSFRADIUS=psfradius,SMAPINFO=smapinfo, $
                   CANCEL=cancel)
  print, 'EXTRACT ', tag, ' cancel', cancel, ' seconds', systime(1)-t0
  for k = 0, n_tags(s)-1 do writefits, out+tag+'_spec_'+strtrim(k,2)+'.fits', s.(k)
  print, 'SIGN_'+tag, apsign
end

pro op_writetrace, path, orders, doorders, naps, tracecoeffs
  ; XSPEXTOOL_PSTRACESPEC's trace file
  z = where(doorders eq 1, gnorders)
  goodorders = orders[z]
  openw, lun, path, /get_lun
  printf, lun, 'OBSMODE=XD'
  printf, lun, 'ORDERS='+strjoin(strtrim(goodorders,2),',')
  printf, lun, 'NAPS='+strtrim(naps,2)
  for i = 0, gnorders-1 do begin
     idx = i*naps
     for j = 0, naps-1 do begin
        printf, lun, 'ORDER'+string(goodorders[i],format='(i2.2)')+'_AP'+string(j+1,format='(i2.2)')+ $
                '= '+strjoin(strtrim(tracecoeffs[*,idx],2),'  ')
        idx = idx + 1
     endfor
  endfor
  close, lun
  free_lun, lun
end

function op_fixtrace, appos, doorders, psnaps
  ; XSPEXTOOL_PSTRACESPEC with Use Ap Positions (appos[j,i], i counting the selected orders)
  tracecoeffs = fltarr(2,total(doorders)*psnaps)
  l = 0
  for i = 0, total(doorders)-1 do begin
     for j = 0, psnaps-1 do begin
        tracecoeffs[*,l] = [appos[j,i],0]
        l = l + 1
     endfor
  endfor
  return, tracecoeffs
end

pro ref_optimal
  sp = getenv('SPEXTOOL_DIR')
  raw = getenv('TSR_ROOT')+'/data/raw/2025-11-02/'
  out = getenv('TSR_ROOT')+'/validation/idl_out/'
  !quiet = 1
  set_plot, 'Z'
  keywords = ['FILENAME','OBJECT','UT','DATE','EXPTIME','RA','DEC','FSAMPLE','COADDS','AIRMASS']
  mc_readflat, out+'flat.fits', flat,ncols,nrows,modename,ps,slith_pix,slith_arc,slitw_pix, $
               slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms,rotation,OMASK=omask
  mc_rdwavecal2d, out+'wavecal_NGC134.fits', wavecal, spatcal, wctype, xo2w, wdeg, odeg, $
                  homeorder, xy2w, xy2s, wxdeg, wydeg, sxdeg, sydeg, disp, ROTATE=rotation
  atr = readfits(sp+'/data/atran4000.fits',/SILENT)
  awave = reform(atr[*,0])
  atrans = reform(atr[*,1])
  ybuffer = 1L
  atmosthresh = 0.4
  psnaps = 2
  step = 5L & sumap = 5L & fitdeg = 1L & winthresh = 3.0 & sigthresh = float(fix('0.1'))
  all = intarr(norders)+1
  ; the triangulations XSPEXTOOL_LOADIMAGES stores for MC_TRACETOAP2D
  for i = 0, norders-1 do begin
     z = where(omask eq orders[i])
     triangulate, double(wavecal[z]), double(spatcal[z]), tri
     name = 'ORD'+string(orders[i],FORMAT='(I2.2)')
     trist = (i eq 0) ? create_struct(name,tri):create_struct(trist,name,tri)
  endfor

  ; ---------------- A-B pair 116-117
  files = raw+['SPEC_SO2025D-006-0116.fits','SPEC_SO2025D-006-0117.fits']
  kna_readts4fits, files, data, hdrinfo, var, KEYWORD=keywords, ROTATE=rotation, /PAIR, ITOT=itot
  data = data/flat & var = var/flat^2
  t0 = systime(1)
  prof = mc_mkspatprof2d(data, omask, wavecal, spatcal, slith_arc/slith_pix, ybuffer, awave, atrans, $
                         atmosthresh, SPATCOEFFS=spc, MEDSPATCOEFFS=mspc)
  print, 'PROFILES seconds', systime(1)-t0
  for i = 0, norders-1 do begin
     writefits, out+'op_spc_y_'+strtrim(i,2)+'.fits', spc.(i*2)
     writefits, out+'op_spc_c_'+strtrim(i,2)+'.fits', spc.(i*2+1)
     writefits, out+'op_mspc_c_'+strtrim(i,2)+'.fits', mspc.(i*2+1)
  endfor

  ; op / om / nb: Auto, traced
  mc_findpeaks, prof, psnaps, all, positions, apsign, /AUTO
  tc = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs, xranges, positions, step, sumap, $
                      winthresh, sigthresh, fitdeg)
  writefits, out+'op_trace.fits', tc
  sign = op_signs(prof, positions, all, psnaps)
  op_extract, 'op', data, var, omask, orders, wavecal, spatcal, ps, sign, tc, all, psnaps, spc, out, /OPT, /BG
  op_extract, 'om', data, var, omask, orders, wavecal, spatcal, ps, sign, tc, all, psnaps, mspc, out, /OPT, /BG
  op_extract, 'nb', data, var, omask, orders, wavecal, spatcal, ps, sign, tc, all, psnaps, spc, out

  ; it: the trace file op's trace would be stored in, imported
  op_writetrace, out+'op_trace.dat', orders, all, psnaps, tc
  mc_readtrace, out+'op_trace.dat', tmode, torders, tnaps, tfitorder, itc
  tdo = intarr(norders)
  match, orders, torders, z
  tdo[z] = 1
  mc_tracetoap2d, omask, wavecal, spatcal, trist, itc, psnaps, orders, tdo, iappos, struc
  writefits, out+'it_appos.fits', iappos
  writefits, out+'it_trace.fits', itc
  mc_findpeaks, prof, psnaps, tdo, junk, isign0, FIXED=iappos
  isign = op_signs(prof, iappos, tdo, psnaps)
  print, 'ITSIGN_FINDPEAKS', isign0
  op_extract, 'it', data, var, omask, orders, wavecal, spatcal, ps, isign, itc, tdo, psnaps, spc, out, /OPT, /BG

  ; gs: Guess positions typed '9.5,28.3' (one list, rebinned to every order)
  g = float(strsplit('9.5,28.3',',',/EXTRACT))
  g = rebin(g, n_elements(g), norders)
  mc_findpeaks, prof, n_elements(g[*,0]), all, gpos, gsign, GUESS=g
  writefits, out+'gs_positions.fits', gpos
  gtc = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs, xranges, gpos, step, sumap, $
                       winthresh, sigthresh, fitdeg)
  writefits, out+'gs_trace.fits', gtc
  gs = op_signs(prof, gpos, all, psnaps)
  op_extract, 'gs', data, var, omask, orders, wavecal, spatcal, ps, gs, gtc, all, psnaps, spc, out, /OPT, /BG

  ; fx: Fix, which ticks Use Ap Positions
  mc_findpeaks, prof, n_elements(g[*,0]), all, fpos, fsign, FIXED=g
  ftc = op_fixtrace(fpos, all, psnaps)
  writefits, out+'fx_trace.fits', ftc
  op_writetrace, out+'fx_trace.dat', orders, all, psnaps, ftc
  fs = op_signs(prof, fpos, all, psnaps)
  op_extract, 'fx', data, var, omask, orders, wavecal, spatcal, ps, fs, ftc, all, psnaps, spc, out, /OPT, /BG

  ; os: Auto on PS Orders 3,5,6 only
  sel = intarr(norders) & sel[[0,2,3]] = 1
  mc_findpeaks, prof, psnaps, sel, spos, ssign, /AUTO
  z = where(sel eq 1)
  stc = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs[*,*,z], xranges[*,z], spos[*,z], step, $
                       sumap, winthresh, sigthresh, fitdeg)
  writefits, out+'os_positions.fits', spos
  writefits, out+'os_trace.fits', stc
  op_writetrace, out+'os_trace.dat', orders, sel, psnaps, stc
  ss = op_signs(prof, spos, sel, psnaps)
  op_extract, 'os', data, var, omask, orders, wavecal, spatcal, ps, ss, stc, sel, psnaps, spc, out, /OPT, /BG

  ; fs (Fix on PS Orders 3,5,6) is in its own run, REF_OPTIMAL_FS: IDL does not return from it
  ; ---------------- am: A mode, frame 116 alone, Fix '9.8,20.0' (Use Ap Positions), optimal.
  ; (Auto with PSNAPS 2 finds a second peak at 8.3 arcsec, next to the star at 9.9: the
  ; apertures overlap and MC_MKAPMASK1D stops with its dialog, so Auto is not a case.)
  kna_readts4fits, files[0], adata, ahdr, avar, KEYWORD=keywords, ROTATE=rotation, ITOT=aitot
  adata = adata/flat & avar = avar/flat^2
  aprof = mc_mkspatprof2d(adata, omask, wavecal, spatcal, slith_arc/slith_pix, ybuffer, awave, atrans, $
                          atmosthresh, SPATCOEFFS=aspc, MEDSPATCOEFFS=amspc)
  ag = float(strsplit('9.8,20.0',',',/EXTRACT))
  ag = rebin(ag, n_elements(ag), norders)
  mc_findpeaks, aprof, psnaps, all, apos, asign0, FIXED=ag
  atc = op_fixtrace(apos, all, psnaps)
  asg = op_signs(aprof, apos, all, psnaps)
  op_extract, 'am', adata, avar, omask, orders, wavecal, spatcal, ps, asg, atc, all, psnaps, aspc, out, /OPT, /BG
  print, 'DONE'
end
