; Reference for xcombspec, xtellcor, xmergeorders and xcleanspec, run on spectra
; extracted exactly as xspextool writes them. All interactive choices are fixed here
; (see validation/compare_dg.py, which makes the same ones).
pro ref_extract_pair, files, fl, wc, aw, at, psmode, opaths, appos_xs
  common refcom, raw, out, keywords, flatname, wcname
  mc_readflat, fl, flat,ncols,nrows,modename,ps,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,norders,orders,edgecoeffs,xranges,rms,rotation,OMASK=omask
  mc_rdwavecal2d, wc, wavecal, spatcal, wctype, xo2w, wdeg, odeg, homeorder, xy2w, xy2s, wxdeg, wydeg, sxdeg, sydeg, disp, ROTATE=rotation
  kna_readts4fits, files, data, hdrinfo, var, KEYWORD=keywords, ROTATE=rotation, /PAIR, ITOT=itot
  data = data/flat & var = var/flat^2
  if psmode then begin
    prof = mc_mkspatprof2d(data, omask, wavecal, spatcal, slith_arc/slith_pix, 1, aw, at, 0.4)
    mc_findpeaks, prof, 2, intarr(norders)+1, positions, apsign, /AUTO
    tracecoeffs = mc_tracespec2d(data, wavecal, spatcal, edgecoeffs, xranges, positions, 5, 5, 3.0, 0.1, 1)
    apsigndat = fltarr(2)
    for i=0,norders-1 do for j=0,1 do begin & profile = prof.(i) & med = median(profile[*,1]) & tabinv,profile[*,0],positions[j,i],idx & apsigndat[j] = (profile[idx,1] gt med) ? 1:-1 & endfor
    apsign = fix(apsigndat) & naps = 2
    apradii = replicate(1.0,2,norders)
    s = mc_extspec2d(data,var,omask,orders,wavecal,spatcal,ps,apsign,tracecoeffs,apradii,PSBGINFO=[2.5,5.0],BGDEG=1)
    psbginfo = [2.5,5.0,1] & apradius = 1.0
  endif else begin
    positions = rebin([appos_xs],1,norders)
    tracecoeffs = fltarr(2,norders) & for i=0,norders-1 do tracecoeffs[*,i] = [positions[0,i],0]
    apsign = intarr(1)+1 & naps = 1
    apradii = rebin([1.0],1,norders)
    s = mc_extspec2d(data,var,omask,orders,wavecal,spatcal,ps,apsign,tracecoeffs,apradii,XSBGINFO=[[2.,6.],[32.,36.]],BGDEG=bgorder)
    xsbginfo = ['2-6,32-36','0'] & apradius = [1.0]
  endelse
  ntags = n_tags(s) & sizes = lonarr(2,ntags)
  for i = 0,ntags-1 do sizes[*,i] = size(s.(i),/DIMEN)
  spec = dblarr(max(sizes), max(sizes[1,*]), ntags)*!values.f_nan
  for i = 0,ntags-1 do spec[0:(sizes[0,i]-1),*,i] = s.(i)
  junk = tag_exist(hdrinfo[0].vals,'FILENAME',INDEX=idx)
  pos = where(apsign eq 1,pos_naps)
  if pos_naps ne 0 then begin
    hdr0 = hdrinfo[0] & aimage = hdr0.vals.(idx) & sky = hdrinfo[1].vals.(idx)
    hdr0.vals.(idx) = file_basename(opaths[0])+'.fits'
    a = replicas(apsign,norders) & z = where(a eq 1)
    if psmode then mc_writespec,spec[*,*,z],opaths[0]+'.fits',aimage,itot[0],sky,flatname,pos_naps,orders,hdr0,positions[pos,*],apradius,modename,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,'um','DN / s','!7k!5 (!7l!5m)','f (!5DN s!u-1!N)','TS4_SOAR','v4.1',PSBGINFO=psbginfo,WAVECAL=wcname,WAVETYPE='Vacuum',DISP=disp,LINEARITY=0,LINCORMAX=80000L $
    else mc_writespec,spec[*,*,z],opaths[0]+'.fits',aimage,itot[0],sky,flatname,pos_naps,orders,hdr0,positions[pos,*],apradius,modename,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,'um','DN / s','!7k!5 (!7l!5m)','f (!5DN s!u-1!N)','TS4_SOAR','v4.1',XSBGINFO=xsbginfo,WAVECAL=wcname,WAVETYPE='Vacuum',DISP=disp,LINEARITY=0,LINCORMAX=80000L
  endif
  neg = where(apsign eq -1,neg_naps)
  if neg_naps ne 0 then begin
    hdr1 = hdrinfo[1] & aimage = hdr1.vals.(idx) & sky = hdrinfo[0].vals.(idx)
    hdr1.vals.(idx) = file_basename(opaths[1])+'.fits'
    a = replicas(apsign,norders) & z = where(a eq -1)
    mc_writespec,spec[*,*,z],opaths[1]+'.fits',aimage,itot[1],sky,flatname,neg_naps,orders,hdr1,positions[neg,*],apradius,modename,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,'um','DN / s','!7k!5 (!7l!5m)','f (!5DN s!u-1!N)','TS4_SOAR','v4.1',PSBGINFO=psbginfo,WAVECAL=wcname,WAVETYPE='Vacuum',DISP=disp,LINEARITY=0,LINCORMAX=80000L
  endif
end

pro ref_combine, files, scaleorder, oname
  ; xcombspec: load + interp to first, scale (Median, middle 80% of order), Robust Weighted Mean thr 8
  mc_readspec,files[0],first,hdr,obsmode,start,stop,norders,naps,orders,xunits,yunits,slith_pix,slith_arc,slitw_pix,slitw_arc,rp,airmass,xtitle,ytitle,instr,/SILENT
  ; first is not passed through mc_caaspecflags: xcombspec_loadspec keeps it as read (double)
  nfiles = n_elements(files) & npix = n_elements(first[*,0,0])
  spc = fltarr(npix,4,nfiles,norders,naps)*!values.f_nan
  ; xcombspec_loadspec: REPLICATE of the first header, never updated per file
  hdrs = replicate(mc_gethdrinfo(hdr),nfiles)
  itot = 0
  for i=0,nfiles-1 do begin
    data = readfits(files[i],h,/SILENT)
    data = mc_caaspecflags(data)
    tmp = fxpar(h,'ITOT') & itot = itot+tmp
    for j=0,norders-1 do begin
      for k=0,naps-1 do begin
        z = where(finite(first[*,0,j*naps+k]) eq 1)
        mc_interpspec,data[z,0,j*naps+k],data[z,1,j*naps+k],first[z,0,j*naps+k],nflx,nerr,IYERROR=data[z,2,j*naps+k],LEAVENANS=1
        mc_interpflagspec,data[z,0,j*naps+k],byte(data[z,3,j*naps+k]),first[z,0,j*naps+k],nfl
        spc[z,0,i,j,k] = first[z,0,j*naps+k]
        spc[z,1,i,j,k] = nflx
        spc[z,2,i,j,k] = nerr
        spc[z,3,i,j,k] = float(nfl)
      endfor
    endfor
  endfor
  zo = where(orders eq scaleorder)
  sspec = fltarr(npix,4,norders*naps)
  for k=0,naps-1 do begin & wave = reform(spc[*,0,0,zo,k]) & wr = max(wave,MIN=wmin,/NAN)-wmin & srange = [wmin+0.1*wr, max(wave,/NAN)-0.1*wr] & z = where(wave gt srange[0] and wave lt srange[1]) & scales = mc_getspecscale(reform(spc[z,1,*,zo,k])) & print, 'SCALES', scales & for j=0,norders-1 do for i=0,nfiles-1 do begin & spc[*,1,i,j,k] = spc[*,1,i,j,k]*scales[i] & spc[*,2,i,j,k] = spc[*,2,i,j,k]*scales[i] & endfor & for j=0,norders-1 do begin & spec = reform(spc[*,1,*,j,k]) & err = reform(spc[*,2,*,j,k]) & mask = intarr(npix,nfiles)+1 & mc_meancomb,spec,mean,mvar,DATAVAR=err^2,MASK=mask,ROBUST=8.0,OGOODBAD=ogoodbad,/SILENT & sspec[*,0,j*naps+k] = spc[*,0,0,j,k] & sspec[*,1,j*naps+k] = mean & sspec[*,2,j*naps+k] = sqrt(mvar) & endfor & endfor
  avehdr = mc_avehdrs(hdrs,TIME_OBS='UT',POSANGLE='ROTANG',HA='HA',AIRMASS='AIRMASS')
  fxhmake,newhdr,sspec
  names = tag_names(avehdr.vals)
  for i=0,n_tags(avehdr.vals)-1 do if names[i] ne 'HISTORY' then fxaddpar,newhdr,names[i],avehdr.vals.(i),avehdr.coms.(i)
  r = {itot:0.0} & r.itot = itot          ; state.r.itot is a FLOAT field
  fxaddpar,newhdr,'ITOT',r.itot,' Total integration time (sec)'
  fxaddpar,newhdr,'FILENAME',file_basename(oname)
  writefits, oname, sspec, newhdr
end

pro ref_dg
  common refcom, raw, out, keywords, flatname, wcname
  sp = getenv('SPEXTOOL_DIR')
  raw  = getenv('TSR_ROOT')+'/data/raw/2025-11-02/'
  out  = getenv('TSR_ROOT')+'/validation/idl_out/'
  ; the instrument file's KEYWORD list, as xspextool reads it (mc_readinstrfile)
  keywords = ['FILENAME','OBJECT','UT','DATE','EXPTIME','RA','DEC','FSAMPLE','COADDS','AIRMASS']
  flatname = 'flat.fits' & wcname = 'wavecal_NGC134.fits'
  fl = out+flatname & wc = out+wcname
  atr = readfits(sp+'/data/atran4000.fits',/SILENT) & aw = reform(atr[*,0]) & at = reform(atr[*,1])
  ; ---- extraction (xspextool), written as Spextool spectra
  ref_extract_pair, raw+['SPEC_SO2025D-006-0116.fits','SPEC_SO2025D-006-0117.fits'], fl, wc, aw, at, 1, out+['spectra0116','spectra0117'], 0
  ref_extract_pair, raw+['SPEC_SO2025D-006-0118.fits','SPEC_SO2025D-006-0119.fits'], fl, wc, aw, at, 1, out+['spectra0118','spectra0119'], 0
  ref_extract_pair, raw+['SPEC_SO2025D-006-0120.fits','SPEC_SO2025D-006-0121.fits'], fl, wc, aw, at, 0, out+['spectra0120','spectra0121'], 20.1
  ref_extract_pair, raw+['SPEC_SO2025D-006-0126.fits','SPEC_SO2025D-006-0127.fits'], fl, wc, aw, at, 0, out+['spectra0126','spectra0127'], 20.1
  ; ---- xcombspec
  ref_combine, out+['spectra0116','spectra0117','spectra0118','spectra0119']+'.fits', 5, out+'comb_HD225200.fits'
  ref_combine, out+['spectra0120','spectra0126']+'.fits', 5, out+'comb_NGC134.fits'
  ; ---- xtellcor (Deconvolution on Pa delta, order 6)
  mc_readspec,out+'comb_HD225200.fits',std,stdhdr,stdobsmode,start,stop,stdnorders,stdnaps,stdorders,stdxunits,stdyunits,slith_pix,slith_arc,slitw_pix,slitw_arc,stdrp,stdairmass,xtitle,ytitle,/SILENT
  std = mc_caaspecflags(std)
  mc_readspec,out+'comb_NGC134.fits',obj,objhdr,objobsmode,start,stop,objnorders,objnaps,objorders,objxunits,objyunits,slith_pix,slith_arc,slitw_pix,slitw_arc,objrp,objairmass,/SILENT
  obj = mc_caaspecflags(obj)
  print, 'AIRMASS', stdairmass, objairmass
  stddisp = mc_fxpar(stdhdr,'DISP*') & print, 'STDDISP', stddisp
  restore, sp+'/data/lvega5.sav'
  wvega = wvin/10000. & cf2vega = fc2vin & cfvega = fcvin & fvega = fvin
  bmag = 6.33 & vmag = 6.31
  ; xmc_conkern on order 6
  idx6 = where(stdorders eq 6)
  i0 = mc_nantrim(std[*,0,idx6],2)
  swave = (reform(std[*,0,idx6]))[i0] & sflux = (reform(std[*,1,idx6]))[i0]
  normreg = [[0.990,0.997],[1.013,1.020]]
  goodbad = intarr(n_elements(swave))
  for i=0,1 do begin & tabinv,swave,normreg[*,i],ii & ii = ii[sort(ii)] & goodbad[round(ii[0]):round(ii[1])] = 1 & endfor
  ncoeffs = robustpoly1d(swave,sflux,3,3,0.01,/GAUSSJ,IGOODBAD=goodbad,/SILENT)
  nflux = float(sflux/poly(swave,ncoeffs))   ; xmc_conkern stores it in a FLOAT structure field
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
  ; scales: the control points untouched (all = scale)
  tellspec = std & npix = n_elements(std[*,0,0])
  for i=0,stdnorders-1 do begin & nkern = readfits(out+'tc_kern_'+strtrim(stdorders[i],2)+'.fits',/SILENT) & scales = replicate(scale, npix) & mc_mktellspec,std[*,0,i],std[*,1,i],std[*,2,i],vmag,(bmag-vmag),nkern,scales,wvega,fvega,cfvega,cf2vega,vshift,tellcor,tellcor_error,scvega & tellspec[*,1,i] = tellcor & tellspec[*,2,i] = tellcor_error & endfor
  writefits, out+'tc_tellspec.fits', tellspec
  ; shifts (autofind) on fixed windows
  sr = [[2.00,2.06],[1.43,1.48],[1.11,1.15],[0.93,0.97],[0.94,1.00]]
  shifts = fltarr(objnorders)
  for i=0,objnorders-1 do begin
    ow = obj[*,0,i] & oflx = obj[*,1,i] & t = mc_nantrim(ow,2) & ow = ow[t] & oflx = oflx[t]
    tw = tellspec[*,0,i] & tf = tellspec[*,1,i] & t = mc_nantrim(tw,2) & tw = tw[t] & tf = tf[t]
    mc_interpspec,tw,tf,ow,nflx
    ; XMC_FINDSHIFT_SETUP normalizes the object by the maximum of its 4-pixel running median
    smooth = median(oflx,4) & oflx = oflx/max(smooth,/NAN)
    z = where(ow gt sr[0,i] and ow lt sr[1,i])
    sh = findgen(151)/50.+(0.-1.5) & rmsv = fltarr(151) & x2 = findgen(n_elements(nflx))
    for k=0,150 do begin
      mc_interpspec,x2+sh[k],nflx,x2,y,e
      d = oflx*y & d = d/median(d,/EVEN)
      mc_moments,d[z],m,v,sdv,/SILENT
      rmsv[k] = sdv
    endfor
    mn = min(rmsv,mi)
    if mi lt 5 or mi gt 145 then begin
      print, 'EDGE', i, mi
      shifts[i] = 0.
    endif else begin
      c = poly_fit1d(sh[(mi-5):(mi+5)],rmsv[(mi-5):(mi+5)],2,/SILENT)
      shifts[i] = -c[1]/2./c[2]
    endelse
  endfor
  print, 'SHIFTS', shifts, format='(a,5e24.16)'
  corspec = obj
  for i=0,objnorders-1 do begin & mc_interpspec,tellspec[*,0,i],tellspec[*,1,i],obj[*,0,i],nflux,nerror,IYERROR=tellspec[*,2,i] & x = findgen(n_elements(obj[*,0,i])) & mc_interpspec,x+shifts[i],nflux,x,nnflux,nnerror,IYERROR=nerror & corspec[*,1,i] = nnflux*obj[*,1,i] & corspec[*,2,i] = sqrt( nnflux^2 * obj[*,2,i]^2 + obj[*,1,i]^2 * nnerror^2 ) & endfor
  writefits, out+'tc_NGC134.fits', corspec
  ; ---- xmergeorders: anchor 3, add 4,5,6,7, scale 1, weighted
  i3 = where(finite(corspec[*,1,0]) eq 1)
  merged = reform(corspec[i3,*,0])
  for o=4,7 do begin & k = where(objorders eq o) & nn = where(finite(corspec[*,1,k]) eq 1) & g = reform(corspec[nn,*,k]) & mc_mergespec,merged[*,0],double(merged[*,1]),g[*,0],double(g[*,1]),owave,oflux,E1=double(merged[*,2]),E2=double(g[*,2]),BF1=merged[*,3],BF2=g[*,3],OERROR=oerror,OBITFLAG=obitflag & merged = [[owave],[oflux],[oerror],[obitflag]] & endfor
  writefits, out+'m_NGC134.fits', merged
  ; ---- xcleanspec: Gaussian FWHM 2.4 and Savitzky-Golay 7/2
  z = where(finite(merged[*,1]) eq 1)
  mc_convolvespec,findgen(n_elements(merged[*,0])),merged[z,1],2.4,cf,ce,ERROR=merged[z,2]
  writefits, out+'cl_gauss.fits', [[cf],[ce]]
  cf = savitzky_golay(merged[z,1],7,DEGREE=2,ERROR=ce2)
  ce2 = merged[z,2] & cf = savitzky_golay(merged[z,1],7,DEGREE=2,ERROR=ce2)
  writefits, out+'cl_sg.fits', [[cf],[ce2]]
  print, 'DONE'
end
