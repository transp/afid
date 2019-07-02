PRO contours, cname, fname, der=der

 IF N_ELEMENTS(der) EQ 0 THEN der = 0

 OPENR, 1, cname
 READF, 1, mumin, mumax
 READF, 1, pmin, pmax
 READF, 1, cmin, emax
 READF, 1, nkbins, npbins
 A = FLTARR(npbins,nkbins)
 READF, 1, A
 CLOSE, 1
 morig = MAX(A)
 PRINT,'A ranges from ',MIN(A),' to ',MAX(A)

 mnzk = nkbins-1
 FOR ik=nkbins-1,0,-1 DO BEGIN
  m = MAX(A[*,ik])
  IF m EQ 0.0 THEN mnzk = ik
 ENDFOR
 PRINT, 'mnzk=',mnzk,' / ',nkbins

 ;WINDOW, 0, XSIZE=1024, YSIZE=570
 SET_PLOT,'PS'
 DEVICE,/COLOR,BITS=8,/INCHES,XSIZE=12.25,YSIZE=5.8
 LOADCT, 33
 !P.MULTI = [0,3,1]
 !P.REGION = [0.0,0.0,0.46,1.0]
 !P.CHARSIZE = 2.0
 xex = FLOOR(ALOG10(emax))
 PRINT, 'A  xex=',xex
 !X.TITLE = 'KE x 10!E'+STRCOMPRESS(STRING(-xex))+'!N'
 yex = FLOOR(ALOG10(pmax - pmin))
 !Y.TITLE = 'P!D!4u!3!N x 10!E'+STRCOMPRESS(STRING(-yex))+'!N'

 PRINT, 'E from ',cmin*mumin,' to ',emax,': range = ',emax - cmin*mumin
 PRINT, 'P from ',pmin,' to ',pmax,': range =',pmax-pmin
 x = (0.1^xex)*((emax - cmin*mumin)*FINDGEN(nkbins)/(nkbins - 1.0) + cmin*mumin)
 !X.RANGE = [MIN(x), x[mnzk]]
 !X.STYLE = 1
 PRINT, 'x range = ',!X.RANGE
 y = (0.1^yex)*(pmin + (pmax - pmin)*FINDGEN(npbins)/(npbins - 1.0))
 !Y.RANGE = [MIN(y), MAX(y)]
 PRINT, 'y range = ',!Y.RANGE
 !Y.STYLE = 1
 nlf=253
 nl = 9

;Take derivative?
 IF der GT 0 THEN BEGIN
   sf = 2.339/MAX(A)
   A = A * sf
   PRINT,'A ranges from ',MIN(A),' to ',MAX(A)
   IF der EQ 1 THEN BEGIN
     !P.TITLE = '!9d!3H/!9d!3P: '+STRCOMPRESS(STRING(mumin))+' < !4l!3 < '+STRCOMPRESS(STRING(mumax))
     C = (A(1:npbins-1,*) - A(0:npbins-2,*)) / (y[1] - y[0])
     xp = x
     yp = 0.5*(y[0:npbins-2] + y[1:npbins-1])
   ENDIF ELSE BEGIN
     !P.TITLE = '!9d!3H/!9d!3E: '+STRCOMPRESS(STRING(mumin))+' < !4l!3 < '+STRCOMPRESS(STRING(mumax))
     C = (A(*,1:nkbins-1) - A(*,0:nkbins-2)) / (x[1] - x[0])
     xp = 0.5*(x[0:nkbins-2] + x[1:nkbins-1])
     yp = y
   ENDELSE
   PRINT, 'Bin derivative ranges from ',MIN(C),' to ',MAX(C)
   scal = MAX(ABS(C))
   scal = 34.3
   flvs = DBLARR(nlf)
   flvs[nlf/2] = 0.0
   flvs[nlf/2 + 1:nlf-1] = 0.1 * scal * 10.0 ^ (FINDGEN(nlf/2)/((nlf/2) - 1.0))
   flvs[0:nlf/2 - 1] = -REVERSE(flvs[nlf/2 + 1:nlf-1])
   CONTOUR, TRANSPOSE(C), xp, yp, LEVELS=flvs, /CELL_FILL
   lvs = DBLARR(nl)
   lvs[nl/2] = 0.0
   lvs[nl/2 + 1:nl-1] = 0.1 * scal * 10.0 ^ (FINDGEN(nl/2)/((nl/2) - 1.0))
   lvs[0:nl/2 - 1] = -REVERSE(lvs[nl/2 + 1:nl-1])
   clab = INTARR(nl)  &  clab[nl/2] = 1
   CONTOUR, TRANSPOSE(C), xp, yp, LEVELS=lvs, /FOLLOW, /OVERPLOT, C_LABELS=clab
   PRINT, 'lvs = ',lvs
   maxa = MAX(C)
 ENDIF ELSE BEGIN
   sf = 2048.0/MAX(A)
   A = A * sf
   PRINT,'A ranges from ',MIN(A),' to ',MAX(A)
   mna = MEAN(A)
   PRINT, 'Mean = ',mna
   !P.TITLE = STRCOMPRESS(STRING(FLOOR(nkbins)))+' x'+STRCOMPRESS(STRING(FLOOR(npbins)))+$
     ' Histogram: '+STRCOMPRESS(STRING(mumin))+' < !4l!3 < '+STRCOMPRESS(STRING(mumax))
   flvs = MAX(A) ^ (FINDGEN(nlf)/(nlf-1.0))
   CONTOUR, TRANSPOSE(A), x, y, LEVELS=flvs, /CELL_FILL
   lvs = MAX(A) ^ (FINDGEN(nl)/(nl-1.0))
   PRINT, 'A contours at ',lvs
   CONTOUR, TRANSPOSE(A), x, y, LEVELS=lvs, /FOLLOW, /OVERPLOT, C_LABELS=INTARR(nl)
   maxa = MAX(A)
 ENDELSE

;Draw legend
 !P.REGION = [0.48,0.0,0.54,1.0]
 !X.TITLE = ''
 xold = !X.RANGE
 !X.RANGE = [0,1]
 xcsold = !X.CHARSIZE
 !X.CHARSIZE = 0.01
 xtkold = !X.ticks
 !X.ticks = 1
 xmrgold = !X.MARGIN
 !X.MARGIN = 3.0
 !y.title = ''
 yold = !Y.RANGE
 ytkold = !Y.ticks
 llvs = (morig/2048.0)*flvs
 IF der EQ 0 THEN BEGIN
   !Y.range = [llvs[0], morig]
   !P.TITLE = 'Level'
   !Y.TICKS= 10; CEIL(ALOG(MAX(llvs))/ALOG(2.0))
   ;PRINT, morig*(2.0^FINDGEN(!Y.TICKS+1))/(2.0^!Y.TICKS)
   CONTOUR, [1,1] # llvs, [0,1], llvs, /FILL, LEVELS=llvs, CHARSIZE=1.5, /YLOG, $
          YTICKV=morig*(2.0^FINDGEN(!Y.TICKS+1))/(2.0^!Y.TICKS); [2.0^FINDGEN(!Y.TICKS), morig]
 ENDIF ELSE BEGIN
   !Y.RANGE = [MIN(flvs), MAX(flvs)]
   !P.TITLE = 'Derivative'
   CONTOUR, [1,1] # flvs, [0,1], flvs, /FILL, LEVELS=flvs, CHARSIZE=1.5
 ENDELSE
 !Y.RANGE = yold
 !X.CHARSIZE = xcsold
 !X.MARGIN =  xmrgold
 !X.ticks = xtkold
 !X.RANGE = xold
 !Y.TICKS = ytkold

 OPENR, 1, fname
 READF, 1, mu
 READF, 1, pmin, pmax
 READF, 1, cmin, emax
 READF, 1, nvals
 B = DBLARR(nvals,nvals)
 READF, 1, B
 CLOSE, 1
 PRINT,'B ranges from ',MIN(B),' to ',MAX(B)
 PRINT, 'E from ',cmin*mu,' to ',emax,': range = ',emax - cmin*mu
 ;PRINT, 'P from ',pmin,' to ',pmax,': range =',pmax-pmin

 ;xex = FLOOR(ALOG10(emax))
 PRINT, 'B  xex=',xex
 !X.TITLE = 'KE x 10!E'+STRCOMPRESS(STRING(-xex))+'!N'
 ;yex = FLOOR(ALOG10(pmax - pmin))
 !Y.TITLE = 'P!D!4u!3!N x 10!E'+STRCOMPRESS(STRING(-yex))+'!N'
 x = (0.1^xex)*((emax - cmin*mu)*FINDGEN(nvals)/(nvals - 1.0) + cmin*mu)
 y = (0.1^yex)*(pmin + (pmax - pmin)*FINDGEN(nvals)/(nvals - 1.0))
 ;!Y.RANGE = [MIN(y), MAX(y)]
 ;!Y.STYLE = 1
 !P.REGION = [0.54,0.0,1.0,1.0]
 ;!X.RANGE = [x[0], x[nvals-1]]
 npcoefs = 5*npbins/8
 IF npcoefs LT 6 THEN npcoefs=6
 nkcoefs = 5*nkbins/8
 IF nkcoefs LT 6 THEN nkcoefs=6

;Plot derivative?
 IF der GT 0 THEN BEGIN
   sf = maxa/MAX(B)
   B = B * sf
   PRINT,'B ranges from ',MIN(B),' to ',MAX(B)
   !P.TITLE = STRCOMPRESS(STRING(FLOOR(nkcoefs)))+' x'+STRCOMPRESS(STRING(FLOOR(npcoefs)))+$
     ' B-spline`: !4l!3 = '+STRCOMPRESS(STRING(mu))
   PRINT, 'Spline derivative ranges from ',MIN(B),' to ',MAX(B)
   scal = MAX(ABS(C))
   CONTOUR, TRANSPOSE(B), x, y, LEVELS=flvs, /CELL_FILL
   PRINT, 'lvs =',lvs
   CONTOUR, TRANSPOSE(B), x, y, LEVELS=lvs, /FOLLOW, /OVERPLOT, C_LABELS=clab
 ENDIF ELSE BEGIN
   sf = mna/MEAN(B)
   B = B * sf
   PRINT,'B ranges from ',MIN(B),' to ',MAX(B)
   PRINT, 'Mean = ',MEAN(B)
   !P.TITLE = STRCOMPRESS(STRING(FLOOR(nkcoefs)))+' x'+STRCOMPRESS(STRING(FLOOR(npcoefs)))+$
     ' B-spline: !4l!3 = '+STRCOMPRESS(STRING(mu))
   PRINT, 'B  x range = ',!X.RANGE
   PRINT, 'B  y range = ',!Y.RANGE
   CONTOUR, TRANSPOSE(B), x, y, LEVELS=flvs, /CELL_FILL
   PRINT, 'B contours at ',lvs
   CONTOUR, TRANSPOSE(B), x, y, LEVELS=[0.0, lvs], /FOLLOW, /OVERPLOT, C_LABELS=[1, INTARR(nl)]
 ENDELSE

 DEVICE,/CLOSE
 SET_PLOT,'x'
END
