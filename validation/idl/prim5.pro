xs=[24.3, 25.0, 25.7, 26.0, 0.0, 3.5]
for v=0,5 do print, xs[v], xs[v] mod 2, not (xs[v] mod 2), (not (xs[v] mod 2)) ? 'TRUE':'FALSE'
np=25.0 & if not np mod 2 then np = np + 1 & print, 'NP25', np
np=24.0 & if not np mod 2 then np = np + 1 & print, 'NP24', np
np=24.6 & if not np mod 2 then np = np + 1 & print, 'NP246', np
exit
