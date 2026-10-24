for v=0,5 do begin & x = [24.3, 25.0, 25.7, 26.0, 0.0, 3.5][v] & print, x, x mod 2, not (x mod 2), (not (x mod 2)) ? 'TRUE':'FALSE' & endfor
sig = findgen(20)^1.5 & k = [0.2,0.5,0.3]
print, 'ET', convol(sig,k,/CENTER,/EDGE_TRUNCATE)
print, 'ETe', convol(sig,[0.1,0.2,0.4,0.3],/CENTER,/EDGE_TRUNCATE)
print, 'PSF', psf_gaussian(NDIMEN=1,NPIXEL=11,FWHM=2.4,/NORMALIZE)
exit
