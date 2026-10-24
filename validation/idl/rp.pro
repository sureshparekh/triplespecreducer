!path=expand_path('+'+getenv('SPEXTOOL_DIR'))+':'+!path
ty = 427.132572279766d + 3d*dindgen(6) & tx = [402.499145300d, 402.475191d, 402.461438d, 402.410198d, 402.364340d, 401.899965554d]
c = robustpoly1d(ty,tx,1,4,0.1,/SILENT,OGOODBAD=g)
print, 'IDLC', c, g, FORMAT='(A,2F20.12,6I3)'
print, 'IDLG', poly(ty[5],c), FORMAT='(A,F20.12)'
exit
