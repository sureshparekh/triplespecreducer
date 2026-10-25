; The 1D (SpeX-type) extraction xspextool makes when the Full Wavecal Name field is blank:
; XSPEXTOOL_LOADIMAGES loads MC_SIMWAVECAL2D's images (wctype '1D'), the profiles and
; traces are the 2D routines on them, and the extraction is MC_EXTPSSPEC / MC_EXTXSSPEC
; with the bad pixel mask (fakebmask.fits, rotated), the linearity bitmask and Bad Pix
; Thresh 7. Inputs typed as MC_CFLD types them (see ref_optimal.pro). No ATRAN: the
; atmosphere mask needs a wavecal file.
;   p1o  point source 116-117, PS BG Subtraction on (optimal, PSF 2.2, background)
;   p1s  the same with PS BG Subtraction off (sum, no background)
;   x1f  extended 120-121 at 20.1, radius 1, background 2-6,32-36 deg 0; the PS
;        Optimal flag is on, so XSPEXTOOL_XSEXTRACTSPEC hands it the spatial map
;        (bad pixels fixed)
;   x1n  the same with PS BG Subtraction off (no spatial map) and XS BG Subtraction off
; Run: cd validation/idl; IDL_PATH="+$SPEXTOOL_DIR:<IDL_DEFAULT>" idl -e ref_1d
pro d1_smap, coeffs, doorders, smapinfo
  idx = 0
  for i = 0, n_elements(doorders)-1 do begin
     if doorders[i] ne 1 then continue
     key = 'Order'+string(i,format='(i2.2)')
     if idx eq 0 then smapinfo = create_struct(key+'Y',coeffs.(i*2),key+'C',coeffs.(i*2+1)) $
     else smapinfo = create_struct(smapinfo,key+'Y',coeffs.(i*2),key+'C',coeffs.(i*2+1))
     idx = idx + 1
  endfor
end

pro d1_write, tag, s, out
  for k = 0, n_tags(s)-1 do writefits, out+tag+'_spec_'+strtrim(k,2)+'.fits', s.(k)
  print, 'WROTE ', tag, n_tags(s)
end

pro ref_1d
  sp = getenv('SPEXTOOL_DIR')
  raw = getenv('TSR_ROOT')+'/data/raw/2025-11-02/'
  out = getenv('TSR_ROOT')+'/validation/idl_out/'
  pkg = sp+'/instruments/TS4_SOAR/'
  !quiet = 1
  set_plot, 'Z'
  keywords = ['FILENAME','OBJECT','UT','DATE','EXPTIME','RA','DEC','FSAMPLE','COADDS','AIRMASS']
  mc_readflat, out+'flat.fits', flat,ncols,nrows,modename,ps,slith_pix,slith_arc,slitw_pix, $
               slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms,rotation,OMASK=omask
  wavecal = mc_simwavecal2d(ncols,nrows,edgecoeffs,xranges,slith_arc)
  spatcal = reform(wavecal[*,*,1])
  wavecal = reform(wavecal[*,*,0])
  bdpxmask = rotate(byte(readfits(pkg+'data/fakebmask.fits',/SILENT)),rotation)
  bitinfo = {lincormax:80000L,lincormaxbit:0}
  bdpxthresh = 7.0
  ybuffer = 1L
  psnaps = 2
  all = intarr(norders)+1
  step = 5L & sumap = 5L & fitdeg = 1L & winthresh = 3.0 & sigthresh = float(fix('0.1'))

  ; ---------------- point source 116-117
  files = raw+['SPEC_SO2025D-006-0116.fits','SPEC_SO2025D-006-0117.fits']
  kna_readts4fits, files, data, hdrinfo, var, KEYWORD=keywords, ROTATE=rotation, /PAIR, $
                   BITINFO=bitinfo, BITMASK=bitmask, ITOT=itot
  data = data/flat & var = var/flat^2
  writefits, out+'d1_bitmask.fits', bitmask
  prof = mc_mkspatprof2d(data, omask, wavecal, spatcal, slith_arc/slith_pix, ybuffer, awave, $
                         atrans, atmosthresh, SPATCOEFFS=spc, MEDSPATCOEFFS=mspc)
  for i = 0, norders-1 do begin
     writefits, out+'d1_prof_'+strtrim(i,2)+'.fits', prof.(i)
     writefits, out+'d1_spc_c_'+strtrim(i,2)+'.fits', spc.(i*2+1)
  endfor
  mc_findpeaks, prof, psnaps, all, positions, apsign, /AUTO
  writefits, out+'d1_positions.fits', positions
  tc = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs, xranges, positions, step, sumap, $
                      winthresh, sigthresh, fitdeg)
  writefits, out+'d1_trace.fits', tc
  apsigndat = fltarr(psnaps)
  for i = 0, norders-1 do begin
     profile = prof.(i)
     for j = 0, psnaps-1 do begin
        med = median(profile[*,1])
        tabinv, profile[*,0], positions[j,i], idx
        apsigndat[j] = (profile[idx,1] gt med) ? 1:-1
     endfor
  endfor
  sign = fix(apsigndat)
  print, 'SIGN', sign
  apradius = 1.0
  d1_smap, spc, all, smap
  t0 = systime(1)
  s = mc_extpsspec(data,var,wavecal,edgecoeffs,tc,norders,psnaps,xranges,slith_arc, $
                   replicate(apradius,psnaps),sign,SPATCOEFFS=smap,PSBGINFO=[2.5,5.0], $
                   BGORDER=1.0,PSFWIDTH=2.2,BDPXMASK=bdpxmask,BITMASK=bitmask, $
                   BDPIXTHRESH=bdpxthresh,CANCEL=cancel)
  print, 'P1O cancel', cancel, systime(1)-t0
  d1_write, 'p1o', s, out
  s = mc_extpsspec(data,var,wavecal,edgecoeffs,tc,norders,psnaps,xranges,slith_arc, $
                   replicate(apradius,psnaps),sign,SPATCOEFFS=nosmap,PSBGINFO=nobg, $
                   BGORDER=nobgorder,PSFWIDTH=nopsf,BDPXMASK=bdpxmask,BITMASK=bitmask, $
                   BDPIXTHRESH=bdpxthresh,CANCEL=cancel)
  print, 'P1S cancel', cancel
  d1_write, 'p1s', s, out

  ; ---------------- extended source 120-121
  files = raw+['SPEC_SO2025D-006-0120.fits','SPEC_SO2025D-006-0121.fits']
  kna_readts4fits, files, xdata, xhdr, xvar, KEYWORD=keywords, ROTATE=rotation, /PAIR, $
                   BITINFO=bitinfo, BITMASK=xbitmask, ITOT=xitot
  xdata = xdata/flat & xvar = xvar/flat^2
  xprof = mc_mkspatprof2d(xdata, omask, wavecal, spatcal, slith_arc/slith_pix, ybuffer, awave, $
                          atrans, atmosthresh, SPATCOEFFS=xspc, MEDSPATCOEFFS=xmspc)
  xtc = fltarr(2,norders)
  for i = 0, norders-1 do xtc[*,i] = [20.1,0]
  bg = strsplit('2-6,32-36',',',/EXTRACT)
  for i = 0, n_elements(bg)-1 do begin
     junk = float(strsplit(bg[i],'-',/EXTRACT))
     bgr = (i eq 0) ? junk:[[bgr],[junk]]
  endfor
  bgfit = 0L
  apr = float(strsplit('1.0',',',/EXTRACT))
  d1_smap, xspc, all, xsmap
  s = mc_extxsspec(xdata,xvar,wavecal,edgecoeffs,xtc,norders,1,xranges,slith_arc,apr, $
                   SPATCOEFFS=xsmap,BGR=bgr,BGORDER=bgfit,BDPXMASK=bdpxmask,BITMASK=xbitmask, $
                   BDPIXTHRESH=bdpxthresh,CANCEL=cancel)
  print, 'X1F cancel', cancel
  d1_write, 'x1f', s, out
  s = mc_extxsspec(xdata,xvar,wavecal,edgecoeffs,xtc,norders,1,xranges,slith_arc,apr, $
                   SPATCOEFFS=nosmap,BGR=nobgr,BGORDER=nobgfit,BDPXMASK=bdpxmask, $
                   BITMASK=xbitmask,BDPIXTHRESH=bdpxthresh,CANCEL=cancel)
  print, 'X1N cancel', cancel
  d1_write, 'x1n', s, out
  print, 'DONE'
end
