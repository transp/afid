PROGRAM mcgen
  USE plasma_state_mod
  IMPLICIT NONE
  INTRINSIC TRIM
  CHARACTER*96 froot
  REAL(KIND=rspec) e_C, m_kg
  REAL(KIND=rspec) vmax,vlel_min,vlel_max,vps_max
  INTEGER :: nparts, ierr
  INTEGER :: nprocs, nthreads, OMP_GET_NUM_PROCS, OMP_GET_NUM_THREADS
  NAMELIST /strings/ froot
  NAMELIST /params/ nparts

!$OMP PARALLEL DEFAULT(SHARED)
  nthreads = OMP_GET_NUM_THREADS()
  nprocs = OMP_GET_NUM_PROCS()
!$OMP END PARALLEL
  WRITE(0,'(A,I3,A,I4,A)') &
       'mcgen running on',nthreads,' thread(s),',nprocs,' proc(s).'

  !Get file root from namelist
  OPEN(8, FILE="mcgen.in", STATUS="OLD", ACTION="READ")
  READ(8, NML=strings)
  READ(8, NML=params)
  CLOSE(8)
  PRINT *,'froot = ',TRIM(froot)
  PRINT *,'nparts = ',nparts

  !Open plasma state, read magnetic axis
  CALL ps_get_plasma_state(ierr,filename=TRIM(froot)//'_ps_ts1_state.cdf')
  PRINT *,'ierr = ',ierr
  CALL ckerr(ierr,'ps_get_plasma_state')
  WRITE(0,'(i5,A,i5,A)') ps%nR,' x ',ps%nZ,' grid'
  WRITE(0,'(A,f6.3,A,f6.3,A)') ' Axis at (',ps%R_axis,', ',ps%Z_axis,')'
  WRITE(0,'(A,f6.3,A)')' |B| = ',ps%B_axis,' on axis.'
  PRINT *,'lcfs R bnds: ',ps%R_min_lcfs,', ',ps%R_max_lcfs
  PRINT *,'lcfs Z bnds: ',ps%Z_min_lcfs,', ',ps%Z_max_lcfs
  PRINT *,'B at lcfs from ',ps%B_min_lcfs,' to ',ps%B_max_lcfs

  !Get particle info from condensed NetCDF file
  CALL my_nc_read(froot, e_C, m_kg, vmax)
  vlel_min = -vmax;  vlel_max = vmax
  vps_max = vmax**2

  !Perform sanity check
  CALL sanity(froot, e_C, m_kg)

  !Generate Jacobian
  CALL writejac_omp(froot, nparts, nthreads, &
       ps%R_min_lcfs,ps%R_max_lcfs,ps%Z_min_lcfs,ps%Z_max_lcfs,&
       vlel_min,vlel_max,vps_max,&
       e_C,m_kg)

END PROGRAM mcgen

!--------------------------------------------------------------------
!Read some attributes from NetCDF file
SUBROUTINE my_nc_read(froot, e_C, m_kg, vmax)
  USE netcdf
  IMPLICIT NONE
  INTRINSIC TRIM
  CHARACTER*96, INTENT(IN) :: froot
  INTEGER, PARAMETER :: rspec = SELECTED_REAL_KIND(9,99)
  REAL(KIND=rspec), INTENT(OUT) :: e_C, m_kg, vmax
  INTEGER ncid, mid, ierr

  ierr = nf90_open(TRIM(froot)//'_condensed.cdf', NF90_NOWRITE, ncid)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF

  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, "vmax_mps", vmax)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  PRINT *,'vmax = ',vmax

  ierr = NF90_INQ_VARID(ncid, "mass", mid)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  ierr = NF90_GET_VAR(ncid, mid, m_kg)  !Reads just 1st value by default
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  PRINT *,'Species 1 mass in kg = ',m_kg

  ierr = NF90_INQ_VARID(ncid, "charge", mid)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  ierr = NF90_GET_VAR(ncid, mid, e_C)  !Reads just 1st value by default
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  PRINT *,'Species 1 charge in C = ',e_C

  ierr = nf90_close(ncid)
END SUBROUTINE my_nc_read

!--------------------------------------------------------------------
SUBROUTINE ckerr(ierr,sbrtn)
  integer, intent(in) :: ierr
  character*(*), intent(in) :: sbrtn

  IF(ierr.NE.0) then
     write(6,*) ' ?plasma_state_test: error in call: '//trim(sbrtn)
     stop
  ENDIF
END SUBROUTINE ckerr

!--------------------------------------------------------------------
!Compute a particle distribution that is uniform in x,v to get
! numerical Jacobian.
SUBROUTINE writejac_omp(froot, nparts, nthreads, &
                        Rmin, Rmax, zmin, zmax, &
                        vlmin, vlmax, vpsmax, ioncharge, ionmass)
  USE plasma_state_mod
  USE netcdf
  IMPLICIT NONE
  INTRINSIC SQRT, TRIM

  CHARACTER*96, INTENT(IN)     :: froot
  INTEGER, INTENT(IN)          :: nparts, nthreads
  REAL(KIND=rspec), INTENT(IN) :: Rmin, Rmax, zmin, zmax
  REAL(KIND=rspec), INTENT(IN) :: vlmin, vlmax, vpsmax
  REAL(KIND=rspec), INTENT(IN) :: ioncharge, ionmass

  REAL(KIND=rspec), PARAMETER  :: one=1.0
  REAL, PARAMETER    :: interval = 5.0 ! wallclock seconds
  INTEGER, PARAMETER :: nfbufsize = 10080

  REAL(KIND=rspec), ALLOCATABLE, DIMENSION(:) :: pphi, mu, Etot
  REAL(KIND=rspec) buffer(3)
  REAL(KIND=rspec) zh, dz, rh, modB, vlh, dvl, bphi
  REAL(KIND=rspec) mot, pmag, tmp
  REAL(KIND=rspec) rhs, drs, vphs, psi_lcfs, psiw
  REAL x, lasttime, thistime
  INTEGER ids(3), dimids(2), strt(2), cnt(2)
  INTEGER ierr, ipart, jbuf, minseed, lparts
  INTEGER ncid, ptcdim, specdim, ppid, muid, Eid, wtid, lid
  INTEGER, DIMENSION(:), ALLOCATABLE :: isv, seed

  WRITE(*,*)
  IF (nthreads.EQ.1) THEN
     WRITE(*,'(A)')'STARTING JAC_OMP WITH ONE THREAD'
  ELSE
     WRITE(*,'(A,I4,A)')'STARTING JAC_OMP WITH',nthreads,' THREADS'
  ENDIF

  ids(1) = ps%id_BRRZ
  ids(2) = ps%id_BZRZ
  ids(3) = ps%id_BphiRZ

  dz = zmax - zmin
  drs = Rmax**2 - Rmin**2
  dvl = vlmax - vlmin
  mot = 0.5*ionmass

  ! Get max flux (at lcfs)
  CALL ps_intrp_2d(ps%R_min_lcfs, ps%z_axis, ps%id_PsiRZ, tmp, ierr)
  CALL ckerr(ierr,'ps_intrp_2d')
  CALL ps_intrp_2d(ps%R_max_lcfs, ps%z_axis, ps%id_PsiRZ, psi_lcfs, ierr)
  CALL ckerr(ierr,'ps_intrp_2d')
  IF (tmp.GT.psi_lcfs) psi_lcfs = tmp
  PRINT *,'psi_lcfs = ',psi_lcfs
  CALL ps_intrp_2d(ps%R_axis, ps%z_axis, ps%id_PsiRZ, tmp, ierr)
  CALL ckerr(ierr,'ps_intrp_2d')
  PRINT *,'psi_axis = ',tmp
  CALL ps_intrp_2d(ps%R_axis, ps%z_axis, ps%id_BphiRZ, tmp, ierr)
  CALL ckerr(ierr,'ps_intrp_2d')
  PRINT *,'B_phi_axis = ',tmp
  psiw = ioncharge*psi_lcfs

  ! Create NetCDF output file, define variables
  ierr = nf90_create(TRIM(froot)//'_jacobian.cdf', NF90_NETCDF4, ncid)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  ierr = nf90_def_dim(ncid, 'nptcl', NF90_UNLIMITED, ptcdim)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_def_dim(ncid, 'nspec', 1, specdim)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  dimids = (/ specdim, ptcdim /)
  strt = (/ 1, 1 /);  cnt = (/ 1, nfbufsize /)
  ierr = nf90_def_var(ncid, 'pphi', NF90_DOUBLE, dimids, ppid, CHUNKSIZES=cnt)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_put_att(ncid, ppid, 'units', 'kg m^2/s')
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_def_var(ncid, 'mu', NF90_DOUBLE, dimids, muid, CHUNKSIZES=cnt)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_put_att(ncid, muid, 'units', 'A m^2')
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_def_var(ncid, 'E', NF90_DOUBLE, dimids, Eid, CHUNKSIZES=cnt)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_put_att(ncid, Eid, 'units', 'J')
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_def_var(ncid, 'weight', NF90_DOUBLE, dimids, wtid, CHUNKSIZES=cnt)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_def_var_fill(ncid, wtid, 0, one)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_def_var(ncid, 'sgn_v', NF90_BYTE, dimids, lid, CHUNKSIZES=cnt)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
  ierr = nf90_enddef(ncid) !Leave define mode, enter data mode
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))

  ALLOCATE(pphi(nfbufsize), mu(nfbufsize), Etot(nfbufsize), isv(nfbufsize))
  WRITE(*,'(A,I7)')' I/O buffer size =',nfbufsize

  ! Initialize seed for random number generator
  CALL RANDOM_SEED(SIZE=minseed) !Returns minimum array length for random seed
  ALLOCATE(seed(minseed))
  seed = 1
  CALL RANDOM_SEED(PUT=seed)

  ! Main loop
  PRINT *,'Generating',nparts,'random particles within LCFS...'
  PRINT *,' Inside     /     Total'
  CALL cpu_time(lasttime)
  ipart = 0
  DO
     lparts =0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP             PRIVATE(seed, x, zh, rhs, rh, pmag, ierr, vphs, vlh, buffer, modB, bphi) &
!$OMP             SCHEDULE(GUIDED, nfbufsize/nthreads) &
!$OMP             REDUCTION(+:lparts)
     DO jbuf=1,nfbufsize
        !Choose random real space coordinates
        CALL RANDOM_NUMBER(x)
        zh = zmin + x*dz
        CALL RANDOM_NUMBER(x)
        rhs = Rmin**2 + x*drs
        rh = SQRT(rhs)

        !Determine if it falls within the lcfs
        CALL ps_intrp_2d(rh, zh, ps%id_PsiRZ, pmag, ierr)
        CALL ckerr(ierr,'ps_intrp_2d')
        IF (pmag.LE.psi_lcfs) lparts = lparts + 1

        !Choose random velocity space coordinates
        CALL RANDOM_NUMBER(x)
        vphs = x*vpsmax
        CALL RANDOM_NUMBER(x)
        vlh = vlmin + x*dvl

        !Convert to constants of motion
        CALL ps_intrp_2d(rh, zh, ids, buffer, ierr) !Retrieve magnetic field
        CALL ckerr(ierr,'ps_intrp_2d')
        modB = SQRT(buffer(1)**2 + buffer(2)**2 + buffer(3)**2)
        bphi = buffer(3)/modB

        Etot(jbuf) = mot*(vphs + vlh**2)
        mu(jbuf) = mot*vphs/modB
        pphi(jbuf) = ioncharge*pmag - ionmass*rh*vlh*bphi  !J.Breslau sign of vlh term
        !pphi(jbuf) = ioncharge*pmag + ionmass*rh*vlh*bphi  !D.Liu sign of vlh term

        !Classify orbit
        if (vlh.GT.0.0) THEN
           isv(jbuf) = 1
        ELSE
           isv(jbuf) = -1
        ENDIF
     ENDDO !jbuf
!$OMP END PARALLEL DO

     ! Flush buffers to NetCDF file
     ierr = nf90_put_var(ncid, ppid, pphi, start=strt, count=cnt)
     IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
     ierr = nf90_put_var(ncid, muid, mu,   start=strt, count=cnt)
     IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
     ierr = nf90_put_var(ncid, Eid,  Etot, start=strt, count=cnt)
     IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
     ierr = nf90_put_var(ncid, lid,  isv,  start=strt, count=cnt)
     IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
     strt(2) = strt(2) + cnt(2)

     ipart = ipart + lparts
     IF (ipart.GE.nparts) EXIT

     CALL cpu_time(thistime)
     IF (thistime - lasttime .GE. interval*nthreads) THEN
        PRINT *,ipart,'/',strt(2)-1
        lasttime = thistime
     ENDIF
  END DO
  PRINT *,strt(2)-1,' total particles written.'

  DEALLOCATE(pphi, mu, Etot, isv, seed)
  ierr = nf90_close(ncid)
  IF (ierr.NE.NF90_NOERR) PRINT *,TRIM(nf90_strerror(ierr))
END SUBROUTINE writejac_omp

!--------------------------------------------------------------------
SUBROUTINE sanity(froot, ioncharge, ionmass)
  USE plasma_state_mod
  USE netcdf
  IMPLICIT NONE
  INTRINSIC TRIM
  CHARACTER*96, INTENT(IN) :: froot
  REAL(KIND=rspec), INTENT(IN) :: ioncharge, ionmass

  REAL(KIND=rspec) buffer(3)
  REAL(KIND=rspec) rh, zh, vlh, vph, modB, bphi, Etot, mu, pphi
  REAL(KIND=rspec) E_cond, mu_cond, pp_cond, pmag
  INTEGER ids(3)
  INTEGER ierr, ncid, vid

  ids(1) = ps%id_BRRZ
  ids(2) = ps%id_BZRZ
  ids(3) = ps%id_BphiRZ

  ierr = nf90_open(TRIM(froot)//'_condensed.cdf', NF90_NOWRITE, ncid)

  !Read phase space coordinates of test particle
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "R1_diag", rh)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "z1_diag", zh)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "vl1_diag", vlh)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "vr1_diag", vph)
  IF (ierr.NE.NF90_NOERR) THEN
     PRINT *,TRIM(nf90_strerror(ierr))
     STOP "Stopped"
  ENDIF

  !Read constants-of-motion coords of test particle
  ierr = NF90_INQ_VARID(ncid, "pphi", vid)
  ierr = NF90_GET_VAR(ncid, vid, pp_cond)
  ierr = NF90_INQ_VARID(ncid, "mu", vid)
  ierr = NF90_GET_VAR(ncid, vid, mu_cond)
  ierr = NF90_INQ_VARID(ncid, "E", vid)
  ierr = NF90_GET_VAR(ncid, vid, E_cond)

  !Convert to constants of motion
  rh = 0.01*rh;  zh = 0.01*zh   !cm to m unit conversion
  CALL ps_intrp_2d(rh, zh, ids, buffer, ierr)
  CALL ckerr(ierr,'ps_intrp_2d')
  modB = SQRT(buffer(1)**2 + buffer(2)**2 + buffer(3)**2)
  bphi = buffer(3)/modB
  write(0,*)'At ',rh,', ',zh,' modB = ',modB

  vlh = 0.01*vlh;  vph = 0.01*vph  !cm/s to m/s unit conversion
  Etot = 0.5*ionmass*(vph**2 + vlh**2)
  write(0,*)'Kinetic energy = ',Etot,' J vs ',E_cond
  write(0,*)' % error = ',100.0*(E_cond - Etot)/E_cond

  mu = 0.5*ionmass*vph*vph/modB
  write(0,*)'Magnetic moment = ',mu,' A m^2 vs ',mu_cond
  write(0,*)' % error = ',100.0*(mu_cond - mu)/mu_cond

  CALL ps_intrp_2d(rh, zh, ps%id_PsiRZ, pmag, ierr)
  CALL ckerr(ierr,'ps_intrp_2d')
  pphi = ioncharge*pmag - ionmass*rh*vlh*bphi
  write(0,*)'P terms = ',ioncharge*pmag,', ',-ionmass*rh*vlh*bphi
  write(0,*)'Angular momentum = ',pphi,' kg m^2/s vs ',pp_cond
  write(0,*)' % error = ',100.0*(pp_cond - pphi)/pp_cond

  ierr = nf90_close(ncid)
END SUBROUTINE sanity
