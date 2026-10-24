y=replicate(0.1,50000L) & print, total(y), FORMAT='(F20.6)'
z=randomu(5,37) & writefits,'/tmp/claude-1000/-home-suresh-PhD-soartriplespec/fa57ebd4-0322-469c-a343-592a6933fbd0/scratchpad/z37.fits',z & print, total(z), total(z^2), FORMAT='(2F20.10)'
print, !cpu.tpool_min_elts
exit
