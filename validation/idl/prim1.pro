!quiet=1
print, 'MED_EVEN_NO', median([1.,2.,3.,4.])
print, 'MED_EVEN_YES', median([1.,2.,3.,4.],/EVEN)
print, 'MED_NAN', median([1.,!values.f_nan,3.,4.,10.])
print, 'MED_NAN_EVEN', median([1.,!values.f_nan,3.,4.],/EVEN)
a = findgen(7,6)^1.5
print, 'MEDFILT'
print, median(a,3)
print, 'MEDDIM', median(findgen(4,3)^2,DIMENSION=1)
print, 'SUB', (findgen(10))[2.7:5.9]
print, 'ROUND', round([0.5,1.5,2.5,-0.5,-1.5])
b = indgen(3,2)
print, 'ROT3'
print, rotate(b,3)
print, 'FIX', fix(2.9), fix(-2.9)
exit
