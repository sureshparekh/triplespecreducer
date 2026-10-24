img = findgen(10,8)^1.3
print, 'NANCOORD', interpolate(img,[2.5,!values.f_nan,3.2],[1.5,2.0,!values.f_nan],/CUBIC)
print, 'CUB1', interpolate(img,[2.3],[3.7],/CUBIC), interpolate(img,[2.3],[3.7],CUBIC=-1.0), interpolate(img,[2.3],[3.7],CUBIC=-0.5)
x = dblarr(30) & y = dblarr(30) & for i=0,29 do x[i] = 1.0d + (i mod 6)*0.123456789d & for i=0,29 do y[i] = (i/6)*0.37
v = x+y
triangulate, x, y, tri
dw = (max(x)-min(x))/5.0d
g = trigrid(x,y,v,tri,[dw,0.37],[min(x),min(y),max(x),max(y)],XGRID=xg,YGRID=yg)
print, 'NXINT', n_elements(xg), n_elements(yg), (max(x)-min(x))/dw
exit
