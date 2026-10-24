x = [0d,1d,0d,1d,0.5d] & y = [0d,0d,1d,1d,0.5d] & v = x+y
triangulate, x, y, tri
foreach_d = [0.2d, 0.2d*(1+1d-12), 0.2d*(1-1d-12), 0.3d, 0.25d, 1d/3d, 0.1d, 0.7d]
for i=0,n_elements(foreach_d)-1 do begin & g = trigrid(x,y,v,tri,[foreach_d[i],0.5d],[0d,0d,1d,1d],XGRID=xg,YGRID=yg) & print, foreach_d[i], 1d/foreach_d[i], n_elements(xg), xg[n_elements(xg)-1], FORMAT='(F22.17,F22.17,I5,F22.17)' & endfor
exit
