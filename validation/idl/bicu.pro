!quiet=1
sp = getenv('SPEXTOOL_DIR')
!path = expand_path('+'+sp) + ':' + !path
z=[1.,2.,4.,3.] & dz1=[0.5,0.2,-0.3,0.1] & dz2=[0.3,-0.1,0.2,0.4] & dz12=[0.01,0.02,-0.03,0.05]
x = findgen(7)#replicate(1.,5) & y = replicate(1.,7)#findgen(5)
r = mc_bicuval(z,dz1,dz2,dz12,0.0,7,0.0,5,x,y)
print, r, FORMAT='(7F12.6)'
exit
