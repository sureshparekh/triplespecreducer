; Reference flat, reproducing ts41cals_mkflat (kna_ts41cals2d.pro) without the display.
@setup
modeinfo = mc_readmodefile(pkg+'data/XD.dat')
lowfiles = raw+'SPEC_DFLAT_'+string(indgen(50)+2,FORMAT='(I4.4)')+'.fits'
lidfiles = raw+'SPEC_DFLAT_'+string(indgen(52)+52,FORMAT='(I4.4)')+'.fits'
kna_readts4fits,lowfiles,data,hdrinfo,var,KEYWORDS=keywords,ROTATE=modeinfo.rotation,NIMAGES=nimages
kna_readts4fits,lidfiles,darkdata,darkhdrinfo,darkvar,KEYWORDS=keywords,ROTATE=modeinfo.rotation,NIMAGES=darknimages
mc_combimgs,darkdata,2,dark,darkvar,THRESH=thresh,CANCEL=cancel
ncols = 2048 & nrows = 1024
dark = rebin(reform(dark,ncols,nrows,1),ncols,nrows,nimages)
darkvar = rebin(reform(darkvar,ncols,nrows,1),ncols,nrows,nimages)
data = temporary(data)-dark
var  = temporary(var)+darkvar
mc_combimgs,data,2,med,mvar,THRESH=thresh,IVCUBE=var,CANCEL=cancel
writefits, out+'flat_med.fits', med
ts4_findorders,med,modeinfo.guesspos,modeinfo.findrange,modeinfo.step,$
               modeinfo.slith_pix,modeinfo.slith_pix_range,modeinfo.edgedeg,$
               modeinfo.norm_ybuffer,modeinfo.flatfrac,modeinfo.comwin,$
               edgecoeffs,xranges,CANCEL=cancel
print, 'EDGES', edgecoeffs
print, 'XRANGES', xranges
ka_scaleflats, data, edgecoeffs, xranges, var, CANCEL=cancel
mc_combimgs,data,2,mean,var,THRESH=thresh,IVCUBE=var,CANCEL=cancel
writefits, out+'flat_mean.fits', mean
mean = mc_normspecflat(mean,edgecoeffs,xranges,modeinfo.slith_arc,$
                       modeinfo.norm_nxg,modeinfo.norm_nyg,$
                       modeinfo.oversamp,modeinfo.norm_ybuffer, $
                       RMS=rms, IVAR=var, OVAR=ovar, MODEL=model,CANCEL=cancel)
bmask = readfits(pkg+'data/fakebmask.fits',/SILENT)
z = where(rotate(bmask,modeinfo.rotation) eq 0,count)
print, 'BADPIX', count
if count ne 0 then mean[z] = 1.0
print, 'RMS', rms
; kna_ts41cals2d: the history names the lamps-on frames; LINCORMAX is xspextool's LONG
history='This flat was created by 1. scaling '+strjoin(hdrinfo[*].vals.(0), ', ')+ $
        ' using the median and then 2. Robust Mean combining the resulting images.'
newhdr = hdrinfo[0]
newhdr.vals.(0) = 'flat.fits'
bitmask=bytarr(ncols,nrows)
mc_writeflat,mc_unrotate(mean,modeinfo.rotation),mc_unrotate(ovar,modeinfo.rotation),$
             mc_unrotate(bitmask,modeinfo.rotation),newhdr,modeinfo.rotation,$
             fix(modeinfo.orders),edgecoeffs,xranges,modeinfo.ps,$
             modeinfo.slith_pix,modeinfo.slith_arc,2.4,1.1,$
             modeinfo.modename,rms,3500,'v4.1',out+'flat.fits',LINCORMAX=80000L,HISTORY=history
writefits, out+'flat_model.fits', model
print, 'DONE'
exit
