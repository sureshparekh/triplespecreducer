; Reference for the xcombspec and xtellcor options beyond the default path: Prune
; Spectra (removed spectra, masked pixel ranges), the scale modes of XMC_SCALESPEC
; (Median, Spectrum, None, a chosen range, one order per aperture), Correct Spectral
; Shape after pruning and with two apertures, xcombspec,/BASIC, the history text;
; xtellcor_basic (a standard of any type, blackbody restore, units), xtellcor's three
; output files (corrected, telluric, modified Vega) and xtellcor_finish.
;
; The real Spextool routines run here: XCOMBSPEC_COMBINESPEC, XCOMBSPEC_CORRECTSPEC,
; XCOMBSPEC_WRITEFILE, XTELLCOR_BASIC_TELLCOR, XTELLCOR_WRITEFILE and
; XTELLCOR_FINISH_TELLCOR, called with a state structure built as the GUIs build it.
; Batch IDL has no widgets, so the widget-free parts of the routines that need them
; (XCOMBSPEC_LOADSPEC, _SCALESPEC, _PRUNESPEC, XTELLCOR_BASIC_LOADSPEC, the GETSHIFT
; routines, XMC_FINDSHIFT_AUTOFIND) are run below without their widget calls, and the choices
; a user makes in XMC_SCALESPEC, XMC_MASKSPEC and XMC_FINDSHIFT are arguments. MC_CFLD
; reads the "fields" from a common block and XVSPEC (the viewer) does nothing.
;
; Run: cd validation/idl; SPEXTOOL_DIR=... TSR_ROOT=... \
;      IDL_PATH="+$SPEXTOOL_DIR:<IDL_DEFAULT>" idl -e "set_plot,'Z' & ref_combtell"
; Outputs: validation/idl_out/ct_*.fits and the printed log (validation/idl_out/ref_combtell.log).

function mc_cfld, field_id, type, EMPTY=empty, CANCEL=cancel
  ; the widget text of field field_id[1], typed as MC_CFLD types it
  common ct_fld, fv
  cancel = 0
  val = fv[field_id[1]]
  case type of
     7: value = val
     1: value = byte(val)
     2: value = fix(val)
     3: value = long(val)
     4: value = float(val)
     5: value = double(val)
  endcase
  return, value
end

pro xvspec, file, _EXTRA=ex
end

pro mc_getfonts, buttonfont, textfont, CANCEL=cancel
  ; the fonts of the widgets (it opens a window to find them)
  cancel = 0 & buttonfont = '' & textfont = ''
end

function ct_xcstate, fullpaths, combineaps, basic
  ; XCOMBSPEC's state (the fields its routines use; TS4_SOAR's instrument values)
  ; filled by XCOMBSPEC_LOADSPEC's statements, its widget calls left out
  w = {basic:basic, modspeccontinue:ptr_new(2), rthresh_fld:[0L,0L], path_fld:[0L,1L], $
       outfile_fld:[0L,2L], xcombspec_base:0L}
  r = {airmass:'AIRMASS', ap:0, combineaps:combineaps, combinestat:0, combspec:ptr_new(2), $
       corspec:ptr_new(2), files:ptr_new(fltarr(2)), fnaps:0, ha:'HA', irafname:'FILENAME', $
       itot:0.0, medrms:ptr_new(2), naps:0, nfiles:0, norders:0, npix:0, nspec:0, $
       orders:ptr_new(2), pixmask:ptr_new(2), posangle:'ROTANG', pruneorderidx:0L, $
       scales:ptr_new(2), scaleorder:ptr_new(0), scaleorderidx:0L, spcmask:ptr_new(2), time:'UT'}
  p = {xtitle:'', ytitle:['','','']}
  d = {ospec:ptr_new(fltarr(2)), wspec:ptr_new(fltarr(2)), hdrinfo:ptr_new(strarr(2))}
  state = {w:w,r:r,d:d,p:p}
  mc_readspec,fullpaths[0],first,hdr,obsmode,start,stop,norders,naps,orders,$
              xunits,yunits,slith_pix,slith_arc,slitw_pix,slitw_arc,$
              rp,airmass,xtitle,ytitle,instr,/SILENT,CANCEL=cancel
  state.r.npix     = n_elements(first[*,0,0])
  state.r.norders  = norders
  state.r.naps     = naps
  state.p.xtitle   = xtitle
  *state.r.orders  = orders
  state.r.nfiles   = n_elements(fullpaths)
  *state.w.modspeccontinue = intarr(naps)
  state.r.ap = 0
  if state.r.combineaps eq 'No' then begin
     state.r.fnaps   = state.r.naps
     *state.r.medrms = fltarr(state.r.norders,state.r.fnaps)+!values.f_nan
     *state.r.scales = fltarr(state.r.nfiles,state.r.fnaps)+1
     state.r.nspec   = state.r.nfiles
     array   = fltarr(state.r.npix,4,state.r.nfiles)*!values.f_nan
     spcmask = intarr(state.r.nfiles,state.r.norders,state.r.fnaps)+1
     pmask   = intarr(state.r.npix,state.r.nfiles)+1
  endif else begin
     state.r.fnaps   = 1
     *state.r.medrms = fltarr(state.r.norders)+!values.f_nan
     *state.r.scales = fltarr(state.r.nfiles*state.r.naps)
     state.r.nspec   = state.r.nfiles*state.r.naps
     array   = fltarr(state.r.npix,4,state.r.nfiles*state.r.naps)+!values.f_nan
     spcmask = intarr(state.r.naps*state.r.nfiles,state.r.norders)+1
     pmask = intarr(state.r.npix,state.r.nfiles*state.r.naps)+1
  endelse
  *state.r.scaleorder   = fltarr(state.r.naps)+!values.f_nan
  state.r.scaleorderidx = 0
  state.r.pruneorderidx = 0
  *state.r.corspec      = intarr(state.r.fnaps)
  key     = 'Order'+string(00,FORMAT='(i2.2)')
  spc     = create_struct(key,array)
  pixmask = create_struct(key,pmask)
  for i = 1, state.r.norders-1 do begin
     key     = 'Order'+string(i,FORMAT='(i2.2)')
     spc     = create_struct(spc,key,array)
     pixmask = create_struct(pixmask,key,pmask)
  endfor
  hdrinfo = replicate(mc_gethdrinfo(hdr),state.r.nfiles)
  spc     = replicate(spc,state.r.fnaps)
  pixmask = replicate(pixmask,state.r.fnaps)
  itot = 0
  for i = 0, state.r.nfiles-1 do begin
     data = readfits(fullpaths[i],hdr,/SILENT)
     data = mc_caaspecflags(data,CANCEL=cancel)
     tmp = fxpar(hdr,'ITOT')
     if state.r.combineaps ne 'Yes' then itot = itot+tmp else itot = itot+tmp*2
     for j = 0, state.r.norders-1 do begin
        for k = 0, state.r.naps-1 do begin
           z = where(finite(first[*,0,j*state.r.naps+k]) eq 1)
           mc_interpspec,data[z,0,j*state.r.naps+k],data[z,1,j*state.r.naps+k], $
                         first[z,0,j*state.r.naps+k],newflux,newerror, $
                         IYERROR=data[z,2,j*state.r.naps+k],LEAVENANS=1
           mc_interpflagspec,data[z,0,j*state.r.naps+k],byte(data[z,3,j*state.r.naps+k]),$
                             first[z,0,j*state.r.naps+k],newflag
           if state.r.combineaps ne 'Yes' then begin
              spc[k].(j)[z,0,i] = first[z,0,j*state.r.naps+k]
              spc[k].(j)[z,1,i] = newflux
              spc[k].(j)[z,2,i] = newerror
              spc[k].(j)[z,3,i] = float(newflag)
           endif else begin
              spc[0].(j)[z,0,i*state.r.naps+k] = first[z,0,j*state.r.naps+k]
              spc[0].(j)[z,1,i*state.r.naps+k] = newflux
              spc[0].(j)[z,2,i*state.r.naps+k] = newerror
              spc[0].(j)[z,3,i*state.r.naps+k] = float(newflag)
           endelse
        endfor
     endfor
  endfor
  *state.d.ospec    = spc
  *state.d.wspec    = spc
  *state.r.spcmask  = spcmask
  *state.r.pixmask  = pixmask
  *state.d.hdrinfo  = hdrinfo
  state.r.itot = itot
  for i = 0,state.r.nfiles-1 do begin
     file = strmid(fullpaths[i],strpos(fullpaths[i],'/',/REVERSE_S)+1)
     sfile = (i eq 0) ? file:[sfile,file]
  endfor
  *state.r.files = sfile
  xcombspec_combinespec,state,CANCEL=cancel
  return, state
end

pro ct_xmc_scalespec, wave, stack, files, mask, scales, wrange, mode, srange, maxspec
  ; XMC_SCALESPEC as a user drives it: it opens on 'Median' over its default range;
  ; 's' sets a range; 'Scale to' None or Spectrum (then the spectrum from its list)
  common xmc_scalespec_state, state
  xmc_scalespec_initcommon,wave,stack,files,mask
  xmc_scalespec_scalespec
  if finite(srange[0]) then begin
     state.r.srange = srange
     xmc_scalespec_scalespec
  endif
  if mode eq 'None' then begin
     state.r.scaleto = 'None'
     xmc_scalespec_scalespec
  endif
  if mode eq 'Spectrum' then begin
     state.r.scaleto = 'None'
     xmc_scalespec_scalespec
     state.r.scaleto = 'Spectrum'
     state.r.maxspec = maxspec
     xmc_scalespec_scalespec
  endif
  scales = state.r.scales
  wrange = state.r.srange
end

pro ct_scalespec, state, mode, srange, maxspec
  ; XCOMBSPEC_SCALESPEC, its widget calls left out
  if (*state.w.modspeccontinue)[state.r.ap] gt 0 then begin
     print, 'SCALE REFUSED'
     return
  endif
  index = state.r.scaleorderidx
  spc = (*state.d.ospec)[state.r.ap]
  good = where((*state.r.spcmask)[*,index,state.r.ap] eq 1)
  zorder = total(where(*state.r.orders eq (*state.r.orders)[index]))
  (*state.r.scaleorder)[state.r.ap] = (*state.r.orders)[index]
  if state.r.combineaps eq 'Yes' then begin
     files = *state.r.files+' Ap '+string(findgen(state.r.naps)+1,FORMAT='(I2.2)')
  endif else files = *state.r.files
  ct_xmc_scalespec,reform(spc.(zorder)[*,0,0]),reform(spc.(zorder)[*,1,*]), $
                   files,reform((*state.r.spcmask)[*,index,state.r.ap]),scales,wrange, $
                   mode,srange,maxspec
  for i = 0, state.r.norders-1 do begin
     s = reform(spc.(i)[*,1,*])
     e = reform(spc.(i)[*,2,*])
     for k = 0, n_elements((*state.r.spcmask)[*,i,state.r.ap])-1 do begin
        s[*,k] = scales[k] * s[*,k]
        e[*,k] = scales[k] * e[*,k]
     endfor
     spc.(i)[*,1,*] = s
     spc.(i)[*,2,*] = e
     mc_meancomb,reform(s[*,good]),mean,var,/RMS
     (*state.r.medrms)[i,state.r.ap] = median(sqrt(var),/EVEN)
     (*state.r.scales)[*,state.r.ap] = scales
  endfor
  (*state.d.wspec)[state.r.ap] = spc
  xcombspec_combinespec,state,CANCEL=cancel
  (*state.w.modspeccontinue)[state.r.ap]=1
  print, 'SCALES', state.r.ap, scales
end

pro ct_xmc_maskspec, wave, stack, specmask, pixmask, remove, regions, regspec
  ; XMC_MASKSPEC as a user drives it: 'Remove' ticks, then for each range the 'Mask'
  ; ticks (regspec[*,r]) and an 's' selection (XMC_MASKSPEC_WINEVENT 'Select'); Accept
  common xmc_maskspec_state, state
  xmc_maskspec_initcommon,wave,stack
  if remove[0] ge 0 then for i = 0, n_elements(remove)-1 do begin
     state.r.plotspec[remove[i]] = 0
     state.r.maskspec[remove[i]] = 0
     state.r.removespec[remove[i]] = 1
  endfor
  nreg = (regions[0] lt 0) ? 0:n_elements(regions)/2
  for r = 0, nreg-1 do begin
     state.r.maskspec = regspec[*,r]
     reg = fltarr(2) & reg[*] = regions[*,r]           ; state.p.reg is FLOAT
     tmp = reg[sort(reg)]
     z = where(state.r.maskspec eq 1,cnt)
     if cnt ne 0 then begin
        good = mc_nantrim(state.d.owave,2,CANCEL=cancel)
        tabinv,state.d.owave[good],tmp,idx
        idx = idx+good[0]
        state.r.pixmask[idx[0]:idx[1],z] = !values.f_nan
        print, 'MASKIDX', idx, format='(a,2f16.8)'
     endif
  endfor
  specmask = state.r.removespec le 0
  pixmask = state.r.pixmask
  z = where(finite(pixmask) eq 0,cnt)
  if cnt ne 0 then pixmask[z] = 0.0
  pixmask = fix(pixmask)
end

pro ct_prunespec, state, index, remove, regions, regspec
  ; XCOMBSPEC_PRUNESPEC, its widget calls left out
  state.r.pruneorderidx = index
  spc    = (*state.d.ospec)[state.r.ap]
  zorder = total(where(*state.r.orders eq (*state.r.orders)[index]))
  y = reform(spc.(zorder)[*,1,*])
  s = size(y,/DIMEN)
  scales = rebin(reform((*state.r.scales)[*,state.r.ap],1,s[1]),s[0],s[1])
  ct_xmc_maskspec,reform(spc.(zorder)[*,0,0]),reform(spc.(zorder)[*,1,*])*scales,$
                  specmask,pixmask,remove,regions,regspec
  (*state.r.spcmask)[*,index,state.r.ap] = specmask
  (*state.r.pixmask)[state.r.ap].(zorder) = pixmask
  xcombspec_combinespec,state,CANCEL=cancel
end

pro ct_correct, state
  ; the real XCOMBSPEC_CORRECTSPEC; says whether MC_SPECCOR ran
  ; (its refusal is a DIALOG_MESSAGE, which needs a display: tested here instead)
  if (*state.w.modspeccontinue)[state.r.ap] ge 3 then begin
     print, 'CORRECT REFUSED'
     return
  endif
  before = (*state.r.corspec)[state.r.ap]
  xcombspec_correctspec,state
  print, 'CORRECTED', state.r.ap, (*state.r.corspec)[state.r.ap], before
end

pro ct_savestack, state, name
  ; the working stacks of the current aperture (the corrected spectra, made with IDL's FFT)
  spc = (*state.d.wspec)[state.r.ap]
  for j = 0, state.r.norders-1 do writefits, name+'_'+strtrim(j,2)+'.fits', spc.(j)
end

pro ct_write, state, out, oname
  common ct_fld, fv
  fv[1] = out & fv[2] = oname
  xcombspec_combinespec,state,CANCEL=cancel
  xcombspec_writefile,state
end

pro ct_mkinputs, out
  ; inputs: IDL's spectra 116-119 with flags (bits 0 and 2 on a pattern of each file's
  ; own), and three two-aperture files made of pairs of them
  for i = 0, 3 do begin
     d = readfits(out+'spectra0'+strtrim(116+i,2)+'.fits',h,/SILENT)
     s = size(d,/DIMEN)
     n = make_array(s[0],4,s[2],TYPE=size(d,/TYPE))
     n[*,0:2,*] = d
     x = lindgen(s[0])
     for k = 0, s[2]-1 do n[*,3,k] = ((x mod 5) eq i)*1 + ((x mod 7) eq (i+k) mod 7)*4
     writefits, out+'ct_flag'+strtrim(i,2)+'.fits', n, h
  endfor
  pairs = [[116,117],[118,119],[117,118]]
  for p = 0, 2 do begin
     a = readfits(out+'spectra0'+strtrim(pairs[0,p],2)+'.fits',h,/SILENT)
     b = readfits(out+'spectra0'+strtrim(pairs[1,p],2)+'.fits',/SILENT)
     s = size(a,/DIMEN)
     n = make_array(s[0],s[1],2*s[2],TYPE=size(a,/TYPE))
     for j = 0, s[2]-1 do begin & n[*,*,2*j] = a[*,*,j] & n[*,*,2*j+1] = b[*,*,j] & endfor
     fxaddpar, h, 'NAPS', 2
     writefits, out+'ct_2ap'+strtrim(p,2)+'.fits', n, h
  endfor
end

pro ref_combtell_xc, out
  common ct_fld, fv
  ; ---- XC1: four flagged spectra; prune three orders; every statistic; scale; correct
  st = ct_xcstate(out+'ct_flag'+['0','1','2','3']+'.fits','No',0)
  ns = 4
  rs = intarr(ns,2) & rs[[0,2],0] = 1 & rs[[0,2,3],1] = 1
  ct_prunespec, st, 2, [1], [[1.20,1.205],[1.25,1.252]], rs        ; order 5
  ct_prunespec, st, 0, [3], [-1], 0                                ; order 3: the last removed
  r1 = intarr(ns,1) & r1[[1,2,3],0] = 1
  ct_prunespec, st, 1, [0], [1.60,1.62], r1                        ; order 4: the first removed
  for s = 0, 8 do begin
     st.r.combinestat = s
     xcombspec_combinespec,st,CANCEL=cancel
     writefits, out+'ct_xc1_s'+strtrim(s,2)+'.fits', *st.r.combspec
  endfor
  st.r.combinestat = 0
  st.r.scaleorderidx = 2
  ct_scalespec, st, 'Median', [!values.f_nan,!values.f_nan], 0
  writefits, out+'ct_xc1_scaled.fits', *st.r.combspec
  ct_scalespec, st, 'Median', [!values.f_nan,!values.f_nan], 0    ; refused
  ; MC_SPECCOR's robust mean as XCOMBSPEC_CORRECTSPEC calls it, for each order
  for j = 0, st.r.norders-1 do begin
     spc = (*st.d.wspec)[0]
     sm = float((*st.r.spcmask)[*,j,0]) & z = where(sm eq 1)
     stack = reform(spc.(j)[*,1,*])
     mc_meancomb,stack[*,z],meanspec,MASK=float((*st.r.pixmask)[0].(j)),OGOODBAD=npixmask,ROBUST=3,/SILENT
     writefits, out+'ct_xc1_sc_'+strtrim(j,2)+'.fits', [[meanspec],[npixmask]]
  endfor
  ct_correct, st
  ct_correct, st                                                   ; refused
  ct_savestack, st, out+'ct_xc1_cor'
  ct_write, st, out, 'ct_xc1'
  ; ---- XC2: Spectrum mode (spectrum 3), a chosen range in order 5; Mean (Std Error); /BASIC
  st = ct_xcstate(out+'ct_flag'+['0','1','2','3']+'.fits','No',1)
  st.r.combinestat = 5
  st.r.scaleorderidx = 2
  ct_scalespec, st, 'Spectrum', [1.21,1.31], 2
  ct_write, st, out, 'ct_xc2'
  ; ---- XC3: 'None' in XMC_SCALESPEC on order 4; Median (MAD); a correction first, then
  ; scaling refused
  st = ct_xcstate(out+'ct_flag'+['0','1','2','3']+'.fits','No',0)
  st.r.combinestat = 6
  ct_correct, st
  st.r.scaleorderidx = 1
  ct_scalespec, st, 'None', [!values.f_nan,!values.f_nan], 0
  st.r.scaleorderidx = 1
  st2 = ct_xcstate(out+'ct_flag'+['0','1','2','3']+'.fits','No',0)
  st2.r.combinestat = 6
  st2.r.scaleorderidx = 1
  ct_scalespec, st2, 'None', [!values.f_nan,!values.f_nan], 0
  ct_write, st2, out, 'ct_xc3'
  ; ---- XC4: two apertures, not combined: each aperture its own scale order; the
  ; correction (MC_SPECCOR refuses the 3-D pixel mask); a removed spectrum in aperture 2
  st = ct_xcstate(out+'ct_2ap'+['0','1','2']+'.fits','No',0)
  st.r.ap = 0 & st.r.scaleorderidx = 2
  ct_scalespec, st, 'Median', [!values.f_nan,!values.f_nan], 0
  st.r.ap = 1 & st.r.scaleorderidx = 1
  ct_scalespec, st, 'Median', [!values.f_nan,!values.f_nan], 0
  st.r.ap = 0
  ct_correct, st
  st.r.ap = 1
  r2 = intarr(3,1) & r2[[1,2],0] = 1
  ct_prunespec, st, 3, [0], [1.05,1.06], r2                         ; order 6
  ct_write, st, out, 'ct_xc4'
  ; ---- XC5: the two apertures combined (six spectra); the sixth removed from order 5
  st = ct_xcstate(out+'ct_2ap'+['0','1','2']+'.fits','Yes',0)
  ct_prunespec, st, 2, [5], [-1], 0
  st.r.combinestat = 4
  st.r.scaleorderidx = 2
  ct_scalespec, st, 'Median', [!values.f_nan,!values.f_nan], 0
  ct_correct, st
  ct_savestack, st, out+'ct_xc5_cor'
  ct_write, st, out, 'ct_xc5'
end


function ct_findshift, objwave, objflux, telflux, awave, atrans, initshift, srange
  ; XMC_FINDSHIFT as a user drives it: XMC_FINDSHIFT_INITCOMMON, the computing part of
  ; XMC_FINDSHIFT_SETUP, a range chosen with 's' (when given), then 'Auto Find'
  ; (XMC_FINDSHIFT_AUTOFIND, its widget calls left out), and Accept
  common xmc_findshift_state, state
  xmc_findshift_initcommon,objwave,objflux,telflux,awave,atrans,INITSHIFT=initshift
  smooth = median(state.d.objflux,4)
  state.d.objflux = state.d.objflux/max(smooth,/NAN)
  smooth = median(1./state.d.telflux,4)
  state.d.dtelflux = 1./state.d.telflux/max(smooth,/NAN)
  state.p.plot1xrange = [min(state.d.objwave,/NAN,MAX=max),max]
  state.r.srange         = state.p.plot1xrange
  if finite(srange[0]) then state.r.srange = srange
  z = where(state.d.objwave gt state.r.srange[0] and state.d.objwave lt $
            state.r.srange[1],count)
  shifts = findgen(151)/50.+(state.r.shift-1.5)
  rms    = fltarr(151)
  for i = 0, 150 do begin
     state.r.shift = shifts[i]
     xmc_findshift_shiftspec
     mc_moments,state.d.divspec[z],mean,var,stddev,/SILENT
     rms[i] = stddev
  endfor
  min = min(rms,minidx)
  del = 5
  if minidx lt del or minidx gt 150-del then begin
     ; IDL stops here with a subscript error: the user would type a shift; 0 is typed
     print, 'EDGE', minidx
     return, 0.0
  endif
  coeff = poly_fit1d(shifts[(minidx-del):(minidx+del)],$
                     rms[(minidx-del):(minidx+del)],2,/SILENT)
  state.r.shift = -coeff[1]/2./coeff[2]
  return, state.r.shift
end

function ct_atran, sp, objrp
  ; the atmosphere as the LOADSPEC routines choose it
  if objrp ge 10000.0 then trp = round(objrp/10000.)*10000
  if objrp lt 10000.0 and objrp ge 1000.0 then trp = round(objrp/1000.)*1000
  if objrp lt 1000.0 and objrp ge 100.0 then trp = round(objrp/100.)*100
  if objrp lt 100.0 and objrp ge 10.0 then trp = 100
  if objrp eq 0 then trp=2000
  return, readfits(filepath('atran'+strtrim(trp,2)+'.fits',ROOT_DIR=sp,SUBDIR='data'),/SILENT)
end

function ct_basicstate, sp, stdfile, objfile, restore, units
  ; XTELLCOR_BASIC's state, filled by XTELLCOR_BASIC_LOADSPEC's statements
  w = {stdspectra_fld:[0L,3L], objoname_fld:[0L,4L], temp_fld:[0L,5L], vmag_fld:[0L,6L], $
       xtellcor_basic_base:0L}
  r = {nytitle:'', objpath:'', restore:restore, shiftobjorder:0, shiftobjap:0, $
       shifts:ptr_new(2), units:units, xtitle:''}
  d = {awave:ptr_new(2), atrans:ptr_new(2), std:ptr_new(fltarr(2)), stdorders:ptr_new(fltarr(2)), $
       objorders:ptr_new(fltarr(2)), obj:ptr_new(fltarr(2)), objnaps:0, objhdr:ptr_new(fltarr(2)), $
       stdhdr:ptr_new(fltarr(2))}
  state = {w:w,r:r,d:d}
  mc_readspec,stdfile,std,stdhdr,stdobsmode,start,stop,stdnorders,stdnaps,$
              stdorders,stdxunits,stdyunits,slith_pix,slith_arc,slitw_pix, $
              slitw_arc,stdrp,stdairmass,stdxtitle,CANCEL=cancel
  std = mc_caaspecflags(std,CANCEL=cancel)
  mc_readspec,objfile,obj,objhdr,objobsmode,start,stop,objnorders,objnaps,$
              objorders,objxunits,objyunits,slith_pix,slith_arc,slitw_pix, $
              slitw_arc,objrp,objairmass,CANCEL=cancel
  obj = mc_caaspecflags(obj,CANCEL=cancel)
  if objairmass eq 0.0 then objairmass = fxpar(objhdr,'TCS_AM')
  *state.d.stdorders = stdorders
  *state.d.objorders = objorders
  state.d.objnaps    = objnaps
  *state.d.std       = std
  *state.d.obj       = obj
  *state.d.objhdr    = objhdr
  *state.d.stdhdr    = stdhdr
  state.r.xtitle     = stdxtitle
  state.r.shiftobjap = 0
  state.r.shiftobjorder = objorders[0]
  *state.r.shifts = fltarr(objnorders,objnaps)
  print, 'BASIC AIRMASS', stdairmass, objairmass, abs(stdairmass-objairmass) gt 0.1
  spec = ct_atran(sp, objrp)
  *state.d.awave = reform(spec[*,0])
  *state.d.atrans = reform(spec[*,1])
  return, state
end

pro ct_basic_getshift, state, order, ap, srange
  ; XTELLCOR_BASIC_GETSHIFT, its widget call left out
  state.r.shiftobjorder = order & state.r.shiftobjap = ap
  zobj = where((*state.d.objorders) eq state.r.shiftobjorder)
  obj_wave = (*state.d.obj)[*,0,zobj*state.d.objnaps+state.r.shiftobjap]
  obj_flux = (*state.d.obj)[*,1,zobj*state.d.objnaps+state.r.shiftobjap]
  trim = mc_nantrim(obj_wave,2)
  obj_wave = obj_wave[trim]
  obj_flux = obj_flux[trim]
  zstd = where((*state.d.stdorders) eq state.r.shiftobjorder)
  tel_wave = (*state.d.std)[*,0,zstd]
  tel_flux = 1./(*state.d.std)[*,1,zstd]
  trim = mc_nantrim(tel_wave,2)
  tel_wave = tel_wave[trim]
  tel_flux = tel_flux[trim]
  mc_interpspec,tel_wave,tel_flux,obj_wave,new_flux,CANCEL=cancel
  shift = ct_findshift(obj_wave,obj_flux,new_flux,*state.d.awave,*state.d.atrans,$
                       total((*state.r.shifts)[zobj,state.r.shiftobjap]),srange)
  (*state.r.shifts)[zobj,state.r.shiftobjap] = shift
end

function ct_finishstate, sp, telfile, objfile
  ; XTELLCOR_FINISH's state, filled by XTELLCOR_FINISH_LOADSPEC's statements
  w = {objoname_fld:[0L,7L], telspectra_fld:[0L,8L], xtellcor_finish_base:0L}
  r = {dairmass:0.0, plotobjorder:0, plotobjap:0, shift:ptr_new(2), xtitle:'', ytitle:''}
  d = {telspec:ptr_new(2), objspec:ptr_new(2), objhdr:ptr_new(2), telhdr:ptr_new(2), $
       telorders:ptr_new(2), objorders:ptr_new(2), objnaps:0, telhistory:ptr_new(2), $
       awave:ptr_new(2), atrans:ptr_new(2)}
  state = {w:w,r:r,d:d}
  mc_readspec,telfile,tel,telhdr,telobsmode,start,stop,telnorders,telnaps,$
              telorders,telxunits,telyunits,slith_pix,slith_arc,slitw_pix, $
              slitw_arc,rp,telairmass,xtitle,ytitle,HISTORY=history,CANCEL=cancel
  mc_readspec,objfile,obj,objhdr,objobsmode,start,stop,objnorders,objnaps,$
              objorders,objxunits,objyunits,slith_pix,slith_arc,slitw_pix, $
              slitw_arc,objrp,objairmass,CANCEL=cancel
  *state.d.telspec = tel
  *state.d.objspec = obj
  *state.d.objhdr  = objhdr
  *state.d.telhdr  = telhdr
  *state.d.telorders = telorders
  *state.d.objorders = objorders
  state.d.objnaps    = objnaps
  state.r.plotobjorder = objorders[0]
  *state.d.telhistory = strjoin(history)
  *state.r.shift = fltarr(objnorders,objnaps)
  state.r.dairmass = telairmass-objairmass
  state.r.xtitle = xtitle
  state.r.ytitle = ytitle
  spec = ct_atran(sp, objrp)
  *state.d.awave = reform(spec[*,0])
  *state.d.atrans = reform(spec[*,1])
  print, 'FINISH AIRMASS', telairmass, objairmass, state.r.dairmass, abs(telairmass-objairmass) gt 0.1
  return, state
end

pro ct_finish_getshift, state, order, ap, srange
  ; XTELLCOR_FINISH_GETSHIFT, its widget call left out
  state.r.plotobjorder = order & state.r.plotobjap = ap
  zz = where((*state.d.objorders) eq state.r.plotobjorder)
  obj_wave = (*state.d.objspec)[*,0,zz*state.d.objnaps+state.r.plotobjap]
  obj_flux = (*state.d.objspec)[*,1,zz*state.d.objnaps+state.r.plotobjap]
  trim = mc_nantrim(obj_wave,2)
  obj_wave = obj_wave[trim]
  obj_flux = obj_flux[trim]
  z        = where((*state.d.telorders) eq state.r.plotobjorder)
  tel_wave = (*state.d.telspec)[*,0,z]
  tel_flux = (*state.d.telspec)[*,1,z]
  trim = mc_nantrim(tel_wave,2)
  tel_wave = tel_wave[trim]
  tel_flux = tel_flux[trim]
  mc_interpspec,tel_wave,tel_flux,obj_wave,new_flux,CANCEL=cancel
  shift = ct_findshift(obj_wave,obj_flux,new_flux,*state.d.awave,*state.d.atrans, $
                       total((*state.r.shift)[zz,state.r.plotobjap]),srange)
  (*state.r.shift)[zz,state.r.plotobjap] = shift
end

pro ref_combtell_tc, sp, out
  common ct_fld, fv
  nan2 = [!values.f_nan,!values.f_nan]
  sr = [[2.00,2.06],[1.43,1.48],[1.11,1.15],[0.93,0.97],[0.94,1.00]]
  ; ---- xtellcor's three files (XTELLCOR_WRITEFILE) from ref_dg.pro's kernel, scale, shifts
  mc_readspec,out+'comb_HD225200.fits',std,stdhdr,stdobsmode,start,stop,stdnorders,stdnaps,stdorders,stdxunits,stdyunits,slith_pix,slith_arc,slitw_pix,slitw_arc,stdrp,stdairmass,xtitle,ytitle,/SILENT
  std = mc_caaspecflags(std)
  mc_readspec,out+'comb_NGC134.fits',obj,objhdr,objobsmode,start,stop,objnorders,objnaps,objorders,objxunits,objyunits,slith_pix,slith_arc,slitw_pix,slitw_arc,objrp,objairmass,/SILENT
  obj = mc_caaspecflags(obj)
  restore, sp+'/data/lvega5.sav'
  wvega = wvin/10000. & cf2vega = fc2vin & cfvega = fcvin & fvega = fvin
  r = {bmag:6.33, vmag:6.31, vshift:0., scale:0., maxdev:0.0, rmsdev:0.0, stdairmass:0.}
  openr, lun, out+'ref_dg.log', /GET_LUN & line = '' & shifts = fltarr(objnorders)
  while not eof(lun) do begin
     readf, lun, line & t = strsplit(line,/EXTRACT)
     if t[0] eq 'VSHIFT' then r.vshift = double(t[1])
     if t[0] eq 'KERNEL' then begin & r.scale = double(t[2]) & r.maxdev = double(t[3]) & r.rmsdev = double(t[4]) & endif
     if t[0] eq 'SHIFTS' then shifts[*] = double(t[1:*])
  endwhile
  free_lun, lun
  npix = n_elements(std[*,0,0])
  match, stdorders, objorders, stdidx
  tellspec = std[*,*,stdidx] & vegaspec = std[*,*,stdidx] & vegaspec[*,2,*] = 1.0
  for i = 0, objnorders-1 do begin
     z = where(stdorders eq objorders[i])
     nkern = readfits(out+'tc_kern_'+strtrim(stdorders[z],2)+'.fits',/SILENT)
     mc_mktellspec,std[*,0,z],std[*,1,z],std[*,2,z],r.vmag,(r.bmag-r.vmag),nkern,replicate(r.scale,npix),wvega,fvega,cfvega,cf2vega,r.vshift,tellcor,tellcor_error,scvega
     tellspec[*,1,i] = tellcor & tellspec[*,2,i] = tellcor_error & vegaspec[*,1,i] = scvega
  endfor
  r.stdairmass = stdairmass
  w = {objoname_fld:[0L,9L], stdspectra_fld:[0L,10L], xtellcor_base:0L}
  rr = {bmag:r.bmag, continue:4L, maxdev:r.maxdev, method:'Deconvolution', $
        nytitle:'!5f!D!7k!N!5 (ergs s!E-1!N cm!E-2!N '+'!5!sA!r!u!9 %!5!n'+'!E-1!N', $
        objpath:out, propstdunc:1L, rmsdev:r.rmsdev, shift:ptr_new(rebin(shifts,objnorders,objnaps)), $
        stdairmass:r.stdairmass, telluricoutput:1, units:'ergs s-1 cm-2 A-1', vegaoutput:1, $
        vmag:r.vmag, vshift:r.vshift}
  d = {dairmass:0.0, objhdr:ptr_new(objhdr), objnaps:objnaps, objorders:ptr_new(objorders), $
       objspec:ptr_new(obj), stdhdr:ptr_new(stdhdr), stdorders:ptr_new(stdorders), $
       tellspec:ptr_new(tellspec), vegaspec:ptr_new(vegaspec)}
  d.dairmass = stdairmass-objairmass
  xs = {w:w, r:rr, d:d}
  fv[9] = 'ct_tc_NGC134' & fv[10] = 'comb_HD225200.fits'
  xtellcor_writefile, xs
  ; ---- XTELLCOR_BASIC, HD 225200 taken as a standard of no particular type
  ; B1: blackbody restore at 9480 K, V = 6.31, ergs s-1 cm-2 A-1; every shift found
  ; over the whole order
  bs = ct_basicstate(sp, out+'comb_HD225200.fits', out+'comb_NGC134.fits', 'Yes', 'ergs s-1 cm-2 A-1')
  for i = 0, 4 do ct_basic_getshift, bs, (*bs.d.objorders)[i], 0, nan2
  print, 'BSHIFTS1', *bs.r.shifts, format='(a,5e24.16)'
  bs.r.objpath = out & fv[3] = 'comb_HD225200.fits' & fv[4] = 'ct_bas1' & fv[5] = '9480' & fv[6] = '6.31'
  xtellcor_basic_tellcor, bs
  ; B2: no restore (a ratio), Jy (converted all the same); a standard without order 5,
  ; so the orders after it stay uncorrected; one shift, in a chosen range
  s5 = readfits(out+'comb_HD225200.fits',h5,/SILENT)
  s5 = s5[*,*,[0,1,3,4]] & fxaddpar, h5, 'ORDERS', '3,4,6,7' & fxaddpar, h5, 'NORDERS', 4
  writefits, out+'ct_std_no5.fits', s5, h5
  bs = ct_basicstate(sp, out+'ct_std_no5.fits', out+'comb_NGC134.fits', 'No', 'Jy')
  ct_basic_getshift, bs, 4, 0, sr[*,1]
  print, 'BSHIFTS2', *bs.r.shifts, format='(a,5e24.16)'
  bs.r.objpath = out & fv[3] = 'ct_std_no5.fits' & fv[4] = 'ct_bas2'
  xtellcor_basic_tellcor, bs
  ; B3: restore at 5800 K, V = 7.5, W m-2 Hz-1; shifts in the ranges of ref_dg.pro
  bs = ct_basicstate(sp, out+'comb_HD225200.fits', out+'spectra0122.fits', 'Yes', 'W m-2 Hz-1')
  for i = 0, 4 do ct_basic_getshift, bs, (*bs.d.objorders)[i], 0, sr[*,i]
  print, 'BSHIFTS3', *bs.r.shifts, format='(a,5e24.16)'
  bs.r.objpath = out & fv[3] = 'comb_HD225200.fits' & fv[4] = 'ct_bas3' & fv[5] = '5800.' & fv[6] = '7.5'
  xtellcor_basic_tellcor, bs
  ; ---- XTELLCOR_FINISH: the telluric file written above applied to spectrum 122
  fs = ct_finishstate(sp, out+'ct_tc_NGC134_tellspec.fits', out+'spectra0122.fits')
  for i = 0, 4 do ct_finish_getshift, fs, (*fs.d.objorders)[i], 0, (i eq 1) ? sr[*,1]:nan2
  print, 'FSHIFTS', *fs.r.shift, format='(a,5e24.16)'
  fv[7] = 'ct_fin_0122' & fv[8] = 'ct_tc_NGC134_tellspec.fits'
  cd, out, CURRENT=here
  xtellcor_finish_tellcor, fs
  cd, here
end

pro ref_combtell
  common ct_fld, fv
  sp = getenv('SPEXTOOL_DIR')
  out = getenv('TSR_ROOT')+'/validation/idl_out/'
  fv = ['8.0', '', '', '', '', '', '', '', '', '', '', '']
  resolve_routine, 'xcombspec', /COMPILE_FULL_FILE
  resolve_routine, 'xmc_scalespec', /COMPILE_FULL_FILE
  resolve_routine, 'xmc_maskspec', /COMPILE_FULL_FILE
  resolve_routine, 'xtellcor', /COMPILE_FULL_FILE
  resolve_routine, 'xtellcor_basic', /COMPILE_FULL_FILE
  resolve_routine, 'xtellcor_finish', /COMPILE_FULL_FILE
  resolve_routine, 'xmc_findshift', /COMPILE_FULL_FILE, /IS_FUNCTION
  resolve_routine, 'ct_after', /COMPILE_FULL_FILE                 ; no plotting
  ct_mkinputs, out
  ref_combtell_xc, out
  ref_combtell_tc, sp, out
  print, 'DONE'
end
