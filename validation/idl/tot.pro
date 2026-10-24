x=[1e8,1.,-1e8] & print, total(x), total(double(x))
y=replicate(0.1,10000000L) & print, total(y), total(y,/DOUBLE)
print, machar() 
exit
