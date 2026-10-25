; Runs the real XMERGEORDERS_AUTOSCALE on a state structure without its
; widgets: the scale is stored in state.r.scale before the procedure
; reaches its first widget call (the scale field), which fails in batch
; mode and is caught here. The caller then makes the procedure's last
; data step itself (XMERGEORDERS_CUTSCALESPEC).
pro tsr_autoscale, state
  catch, err
  if err ne 0 then begin
     catch, /cancel
     return
  endif
  xmergeorders_autoscale, state
end
