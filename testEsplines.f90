PROGRAM testEsplines
  IMPLICIT NONE
  REAL*8 pmin, pmax, mumin, mumax, cmin, emin, emax, mu
  REAL*8 val, dke, dp, ddp, ddk
  INTEGER, PARAMETER :: npts=400
  INTEGER nmubinsn, nmubinsp, ike, ipp, der, sgnv, iorder
  NAMELIST /params/ mu, der, sgnv, iorder

  !Read input options:
  ! mu = magnetic moment
  ! der = 0: distribution function f(p,mu,E)
  ! der = 1: df/dp
  ! der = 2: df/dE
  mu=0.0; der=0; sgnv=1; iorder=0  ! Default values
  OPEN(8, FILE="testEsplines.in", STATUS="OLD", ACTION="READ")
  READ(8, NML=params)
  CLOSE(8)
  PRINT *,'mu=',mu,',  der=',der,', sign(v)=',sgnv,', iorder=',iorder

  !Read spline data
  PRINT *,'Reading data...'
  CALL pspline_init(nmubinsn, nmubinsp)
  PRINT *,nmubinsn,' negative v bins and ',nmubinsp,' positive v bins read, normalized.'
  CALL pspline_set_order(iorder)

  !Find global bounds
  CALL getpspline3bounds(pmin, pmax, mumin, mumax, cmin, emax)
  PRINT *,'Global p_phi from ',pmin,' to ',pmax
  PRINT *,'       mu from ',mumin,' to ',mumax
  PRINT *,'       KE from ',cmin*mumin,' to ',emax
  IF ((mu.LT.mumin).OR.(mu.GT.mumax)) STOP

  !Find local bounds
  CALL getpspline2bounds(mu, sgnv, pmin, pmax, cmin, emax)
  PRINT *,'Local p_phi from ',pmin,' to ',pmax
  PRINT *,'      KE from ',cmin*mu,' to ',emax

  !Set increment for interpolation grid
  dp = (pmax - pmin)/(npts - 1)
  emin = cmin*mu
  dke = (emax - emin)/(npts - 1)

  IF (der.eq.0) THEN
     CALL getpdf(pmin + (npts/2)*dp, mu, emin + (npts/2)*dke, sgnv, val)
  ELSE
     CALL getpdfd(pmin + (npts/2)*dp, mu, emin + (npts/2)*dke, sgnv, val, ddp, ddk)
     IF (der.eq.1) THEN
        PRINT *,'Differentiating with respect to Pphi.'
        val = ddp
     ELSE
        PRINT *,'Differentiating with respect to E.'
        val = ddk
     END IF
  END IF
  PRINT *,'middle val=',val

  !Create output file, write header
  OPEN(9, FILE="mufine", STATUS="REPLACE", ACTION="WRITE")
  WRITE(9,'(E16.8)') mu
  WRITE(9,'(2E16.8)') pmin, pmax
  WRITE(9,'(2E16.8)') cmin, emax
  WRITE(9,'(I12)') npts

  !Loop over pphi, ke to interpolate and write values to file
  DO ike=1,npts
     DO ipp=1,npts
        IF (der.eq.0) THEN
           CALL getpdf(pmin + ipp*dp, mu, emin + ike*dke, sgnv, val)
        ELSE IF (der.eq.1) THEN
           CALL getdfdp(pmin + ipp*dp, mu, emin + ike*dke, sgnv, val, ddp)
           val = ddp
        ELSE
           CALL getdfde(pmin + ipp*dp, mu, emin + ike*dke, sgnv, val, ddk)
           val = ddk
        END IF
        WRITE(9,'(E16.8)') val
     END DO !ipp
  END DO !ike

  CLOSE(9)
  CALL pspline_free
END PROGRAM testEsplines
