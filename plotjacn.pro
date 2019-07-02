PRO plotjacn, fname, nparts=nparts, sgn=sgn
IF N_ELEMENTS(fname) EQ 0 THEN fname = '37065Y022_jacobian.cdf'
IF N_ELEMENTS(sgn) EQ 0 THEN sgn=0

WINDOW,0,XSIZE=1347,YSIZE=1010

ncid = NCDF_OPEN(fname)
ptcid = NCDF_DIMID(ncid, 'nptcl')
IF ptcid EQ -1 THEN BEGIN
    PRINT,'Dimension not found'
    RETURN
ENDIF
NCDF_DIMINQ, ncid, ptcid, dummy, nparts
PRINT,nparts,' particles in file'
pdata = DBLARR(3,nparts)

;Read toroidal angular momentum
vid = NCDF_VARID(ncid, 'pphi')
IF vid EQ -1 THEN BEGIN
    PRINT,'Variable pphi not found'
    RETURN
ENDIF
NCDF_VARGET, ncid, vid, buf, COUNT=[1,nparts]
pdata[0,*] = buf

;Read magnetic moment
vid = NCDF_VARID(ncid, 'mu')
IF vid EQ -1 THEN BEGIN
    PRINT,'Variable mu not found'
    RETURN
ENDIF
NCDF_VARGET, ncid, vid, buf, COUNT=[1,nparts]
pdata[1,*] = buf

;Read kinetic energy
vid = NCDF_VARID(ncid, 'E')
IF vid EQ -1 THEN BEGIN
    PRINT,'Variable E not found'
    RETURN
ENDIF
NCDF_VARGET, ncid, vid, buf, COUNT=[1,nparts]
pdata[2,*] = buf

IF sgn NE 0 THEN BEGIN
 ;Read sign of v||
 vid = NCDF_VARID(ncid, 'sgn_v')
 IF vid EQ -1 THEN BEGIN
    PRINT,'Variable sgn_v not found'
    RETURN
 ENDIF
 sgnv = BYTARR(nparts)
 NCDF_VARGET, ncid, vid, sgnv, COUNT=[1,nparts]
 IF sgn GT 0 THEN $
   svidx = WHERE(sgnv EQ sgn, count) $
 ELSE $
   svidx = WHERE(sgnv EQ sgn+256, count)
 PRINT,count,' particles have sgn(v) =',sgn
 buf = pdata[0,svidx]  &  pdata[0,0:count-1] = buf
 buf = pdata[1,svidx]  &  pdata[1,0:count-1] = buf
 buf = pdata[2,svidx]  &  pdata[2,0:count-1] = buf
 nparts = count
ENDIF
NCDF_CLOSE, ncid

PRINT, 'particle data read.'

;Plots
!P.MULTI = [0,3,2]

;1D Histograms
maxbins = MIN([FLOOR(nparts/100.0), 128])
!P.TITLE = '1D Histogram, '+STRCOMPRESS(STRING(nparts/maxbins))+' ppb'
!P.CHARSIZE = 3.2

!Y.TITLE = '!9II!3f'
!Y.RANGE = [0,0]

;Toroidal angular momentum
pmax = MAX(pdata[0,*])  &  pmin = MIN(pdata[0,*])
dxip = maxbins/(pmax - pmin)
prescale = FLOOR(ALOG10(pmax - pmin))
!X.TITLE = 'Angular momentum P!D!4u!3!N (kg m!E2!N s!E-1!N) x 10!E'+STRCOMPRESS(STRING(-prescale))+'!N'
ptitle = !X.TITLE
!X.STYLE = 1
PRINT, 'toroidal angular momentum in',maxbins,' bins of width ',1.0/dxip
phist = DBLARR(maxbins+1)
FOR ipart=0L,nparts-1L DO BEGIN
    ibin = FLOOR(dxip*(pdata(0,ipart) - pmin))
    phist[ibin] = phist[ibin] + 1.0
ENDFOR
phist[maxbins - 1] = phist[maxbins - 1] + phist[maxbins]
ptot = MEAN(phist[0:maxbins-1])*maxbins
phist[0:maxbins-1] = phist[0:maxbins-1]/ptot
px = (FINDGEN(maxbins)/dxip + pmin)*(1d-1^prescale)
PLOT, px, phist[0:maxbins-1], PSYM=10

;Magnetic moment
phist = 0.0*phist
mmax = MAX(pdata(1,*))  &  mmin = MIN(pdata(1,*))
dxim = maxbins/(mmax - mmin)
mrescale = FLOOR(ALOG10(mmax - mmin))
!X.TITLE = 'Magnetic moment !4l!3 (A m!E2!N) x 10!E'+STRCOMPRESS(STRING(-mrescale))+'!N'
mtitle = !X.TITLE
PRINT, 'toroidal angular momentum in',maxbins,' bins of width ',1.0/dxim
FOR ipart=0L,nparts-1L DO BEGIN
    ibin = FLOOR(dxim*(pdata(1,ipart) - mmin))
    phist[ibin] = phist[ibin] + 1.0
ENDFOR
phist[maxbins - 1] = phist[maxbins - 1] + phist[maxbins]
ptot = MEAN(phist[0:maxbins-1])*maxbins
phist[0:maxbins-1] = phist[0:maxbins-1]/ptot
PLOT, (FINDGEN(maxbins)/dxim + mmin)*(1D-1^mrescale), phist[0:maxbins-1], PSYM=10

;Kinetic Energy
phist = 0.0*phist
emax = MAX(pdata(2,*))  &  emin = MIN(pdata(2,*))
dxie = maxbins/(emax - emin)
erescale = FLOOR(ALOG10(emax - emin))
!X.TITLE = 'Kinetic energy E (J) x 10!E'+STRCOMPRESS(STRING(-erescale))+'!N'
etitle = !X.TITLE
PRINT, 'kinetic energy in',maxbins,' bins of width ',1.0/dxie
FOR ipart=0L,nparts-1L DO BEGIN
    ibin = FLOOR(dxie*(pdata(2,ipart) - emin))
    phist[ibin] = phist[ibin] + 1.0
ENDFOR
phist[maxbins - 1] = phist[maxbins - 1] + phist[maxbins]
ptot = MEAN(phist[0:maxbins-1])*maxbins
phist[0:maxbins-1] = phist[0:maxbins-1]/ptot
PLOT, (FINDGEN(maxbins)/dxie + emin)*(1D-1^erescale), phist[0:maxbins-1], PSYM=10

;2D Histograms
LOADCT, 9

maxbins = MIN([96, FLOOR(0.05*SQRT(nparts))])
PRINT, maxbins,' x ',maxbins,' bins'
dxip = maxbins/(pmax - pmin)
dxim = maxbins/(mmax - mmin)
thist = DBLARR(maxbins+1,maxbins+1)
!P.TITLE = '2D Histogram,'+STRCOMPRESS(STRING(nparts/(maxbins*maxbins)))+' ppb'
FOR ipart=0L,nparts-1L DO BEGIN
    ibinm = FLOOR(dxim*(pdata(1,ipart) - mmin))
    ibinp = FLOOR(dxip*(pdata(0,ipart) - pmin))
    thist[ibinp,ibinm] = thist[ibinp,ibinm] + 1.0
ENDFOR
FOR ibin=0,maxbins DO BEGIN
    thist[ibin,maxbins-1] = thist[ibin,maxbins-1] + thist[ibin,maxbins]
ENDFOR
FOR ibin=0,maxbins-1 DO BEGIN
    thist[maxbins-1,ibin] = thist[maxbins-1,ibin] + thist[maxbins,ibin]
ENDFOR
PRINT, '2d tot:',MEAN(thist[0:maxbins-1,0:maxbins-1])*maxbins^2
PRINT, '2d peak:',MAX(thist[0:maxbins-1,0:maxbins-1])
!X.TITLE = ptitle  &  !Y.TITLE = mtitle
!y.range = [0,0]
PRINT,'maxbins=',maxbins,', dxim=',dxim
PRINT, 'mmax=',mmax,mmax-mmin
PRINT, 'offset =',(maxbins - 1.0)/dxim
PRINT, 'ymax=',(1d-1^mrescale)*((maxbins-1.0)/dxim + mmin)
CONTOUR,ALOG(thist[0:maxbins-1,0:maxbins-1]+1),(1d-1^prescale)*(FINDGEN(maxbins)/dxip + pmin),$
  (1d-1^mrescale)*(FINDGEN(maxbins)/dxim + mmin), /FILL, NLEVELS=255
;PRINT, 'Pearson correlation coefficient for p/m: ',CORRELATE(pdata[0,*], pdata[1,*])

dxie = maxbins/(emax - emin)
thist = 0.0*thist
FOR ipart=0L,nparts-1L DO BEGIN
  ibine = FLOOR(dxie*(pdata(2,ipart) - emin))
  ibinp = FLOOR(dxip*(pdata(0,ipart) - pmin))
  thist[ibinp,ibine] = thist[ibinp,ibine] + 1.0
ENDFOR
FOR ibin=0,maxbins DO BEGIN
    thist[ibin,maxbins-1] = thist[ibin,maxbins-1] + thist[ibin,maxbins]
ENDFOR
FOR ibin=0,maxbins-1 DO BEGIN
    thist[maxbins-1,ibin] = thist[maxbins-1,ibin] + thist[maxbins,ibin]
ENDFOR
PRINT, '2d tot:',MEAN(thist[0:maxbins-1,0:maxbins-1])*maxbins^2
PRINT, '2d peak:',MAX(thist[0:maxbins-1,0:maxbins-1])
!Y.TITLE = etitle
CONTOUR,ALOG(thist[0:maxbins-1,0:maxbins-1]+1),(1d-1^prescale)*(FINDGEN(maxbins)/dxip + pmin),$
  (1d-1^erescale)*(FINDGEN(maxbins)/dxie + emin), /FILL, NLEVELS=255
;PRINT, 'Pearson correlation coefficient for p/e: ',CORRELATE(pdata[0,*], pdata[2,*])

thist = 0.0*thist
FOR ipart=0L,nparts-1L DO BEGIN
    ibine = FLOOR(dxie*(pdata(2,ipart) - emin))
    ibinm = FLOOR(dxim*(pdata(1,ipart) - mmin))
    thist[ibine,ibinm] = thist[ibine,ibinm] + 1.0
ENDFOR
FOR ibin=0,maxbins DO BEGIN
    thist[ibin,maxbins-1] = thist[ibin,maxbins-1] + thist[ibin,maxbins]
ENDFOR
FOR ibin=0,maxbins-1 DO BEGIN
    thist[maxbins-1,ibin] = thist[maxbins-1,ibin] + thist[maxbins,ibin]
ENDFOR
PRINT, '2d tot:',MEAN(thist[0:maxbins-1,0:maxbins-1])*maxbins^2
PRINT, '2d peak:',MAX(thist[0:maxbins-1,0:maxbins-1])
!X.TITLE = etitle  &  !Y.TITLE = mtitle
CONTOUR,ALOG(thist[0:maxbins-1,0:maxbins-1]+1),(1d-1^erescale)*(FINDGEN(maxbins)/dxie + emin),$
  (1d-1^mrescale)*(FINDGEN(maxbins)/dxim + mmin), /FILL, NLEVELS=255
 ;PRINT, 'Pearson correlation coefficient for e/m: ',CORRELATE(pdata[2,*], pdata[1,*])

END
