; Shared setup for the IDL reference runs: Spextool 4.1 + TS4_SOAR on the path.
; Set SPEXTOOL_DIR (the Spextool 4.1 folder) and TSR_ROOT (this repository) first.
sp = getenv('SPEXTOOL_DIR')
!path = expand_path('+'+sp) + ':' + !path
!quiet = 1
set_plot, 'Z'
raw  = getenv('TSR_ROOT')+'/data/raw/2025-11-02/'
out  = getenv('TSR_ROOT')+'/validation/idl_out/'
pkg  = sp+'/instruments/TS4_SOAR/'
keywords = ['FILENAME','OBJECT','UT','DATE','EXPTIME','RA','DEC','FSAMPLE','COADDS','AIRMASS']
