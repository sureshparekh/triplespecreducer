!path=expand_path('+'+getenv('SPEXTOOL_DIR'))+':'+!path
f = getenv('TSR_ROOT')+'/data/raw/2025-11-02/SPEC_SO2025D-006-0116.fits'
kna_readts4fits,[f],d0,h,v,KEYWORDS=['FILENAME'],ROTATE=0
r = readfits(f,/SILENT)
print, 'SIZE0', size(d0,/DIMEN), '  RAW', size(r,/DIMEN)
x = r/5.0
print, 'EQUAL_AS_RAW_LAYOUT', array_equal(reform(d0,1024,2048), x)
print, 'SAMPLE', d0[0:3,0], x[0:3,0]
writefits, getenv('TSR_ROOT')+'/validation/idl_out/rot0_116.fits', d0
exit
