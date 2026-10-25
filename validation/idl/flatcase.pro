; One run of ts41cals_mkflat (kna_ts41cals2d.pro) without the display, for
; the branch tests (ref_flat_branches.pro): lidfiles may be '' (no darks).
pro flatcase, lowfiles, lidfiles, flatoname, out, pkg, keywords
  modeinfo = mc_readmodefile(pkg+'data/XD.dat')
  ncols = 2048 & nrows = 1024
  kna_readts4fits,lowfiles,data,hdrinfo,var,KEYWORDS=keywords,ROTATE=modeinfo.rotation,NIMAGES=nimages
  if strtrim(lidfiles[0],2) ne '' then begin
     kna_readts4fits,lidfiles,darkdata,darkhdrinfo,darkvar,KEYWORDS=keywords,ROTATE=modeinfo.rotation,NIMAGES=darknimages
     if darknimages ne 1 then begin
        mc_combimgs,darkdata,2,dark,darkvar,/UPDATE,THRESH=thresh,CANCEL=cancel
     endif else dark = darkdata
     dark = rebin(reform(dark,ncols,nrows,1),ncols,nrows,nimages)
     darkvar = rebin(reform(darkvar,ncols,nrows,1),ncols,nrows,nimages)
     data = temporary(data)-dark
     var  = temporary(var)+darkvar
  endif
  if nimages ne 1 then begin
     mc_combimgs,data,2,med,mvar,/UPDATE,THRESH=thresh,IVCUBE=var,CANCEL=cancel
  endif else med = data
  writefits, out+flatoname+'_med.fits', med
  ts4_findorders,med,modeinfo.guesspos,modeinfo.findrange,modeinfo.step,$
                 modeinfo.slith_pix,modeinfo.slith_pix_range,modeinfo.edgedeg,$
                 modeinfo.norm_ybuffer,modeinfo.flatfrac,modeinfo.comwin,$
                 edgecoeffs,xranges,CANCEL=cancel
  if nimages ne 1 then begin
     ka_scaleflats, data, edgecoeffs, xranges, var, CANCEL=cancel
     mc_combimgs,data,2,mean,var,/UPDATE,THRESH=thresh,IVCUBE=var,CANCEL=cancel
  endif else mean = data
  mean = mc_normspecflat(mean,edgecoeffs,xranges,modeinfo.slith_arc,$
                         modeinfo.norm_nxg,modeinfo.norm_nyg,$
                         modeinfo.oversamp,modeinfo.norm_ybuffer, $
                         RMS=rms, IVAR=var, OVAR=ovar, MODEL=model,CANCEL=cancel)
  mask = byte(readfits(pkg+'data/fakebmask.fits',/SILENT))
  z = where(rotate(mask,modeinfo.rotation) eq 0,count)
  if count ne 0 then mean[z] = 1.0
  junk = tag_exist(hdrinfo[*].vals,'FILENAME',INDEX=idx)
  history='This flat was created by 1. scaling '+strjoin(hdrinfo[*].vals.(idx), ', ')+ $
          ' using the median and then 2. Robust Mean combining the resulting images.'
  newhdr = hdrinfo[0]
  newhdr.vals.(idx) = flatoname+'.fits'
  bitmask=bytarr(ncols,nrows)
  ; state.w.version is data/Version.dat ('4.1'); state.r.lincormax the LONG of LINCORMAX
  mc_writeflat,mc_unrotate(mean,modeinfo.rotation),mc_unrotate(ovar,modeinfo.rotation),$
               mc_unrotate(bitmask,modeinfo.rotation),newhdr,modeinfo.rotation,$
               fix(modeinfo.orders),edgecoeffs,xranges,modeinfo.ps,$
               modeinfo.slith_pix,modeinfo.slith_arc,2.4,1.1,$
               modeinfo.modename,rms,3500,'4.1',out+flatoname+'.fits',LINCORMAX=80000L,HISTORY=history
  print, 'WROTE ', flatoname
end
