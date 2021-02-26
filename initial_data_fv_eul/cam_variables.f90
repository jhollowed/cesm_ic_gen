 module cam_variables

!=======================================================================
!     CAM3 grid and resolutions
!=======================================================================
      use cam_grid
      
      implicit none

!=======================================================================
!     CAM3 model variables (here with constant initial values)
!=======================================================================
      real*8, parameter ::  ts1 = 288.d0,                              &
                            ts2 = ts1,                                 &
                            ts3 = ts1,                                 &
                            ts4 = ts1,                                 &
                            snowh = 0.d0,                              &
                            oro   = 0.d0,                              &
                            tsice = 0.d0,                              &
                            sgh   = 0.d0,                              &
                            fland = 0.d0,                              &
                            sicthk = 0.d0,                             &
                            cwat   = 0.d0,                             &
                            landm  = 0.d0

!=======================================================================
!     allocatable arrays
!=======================================================================
      real*8, allocatable :: u(:,:,:),                                 & ! zonal wind
                             v(:,:,:),                                 & ! meridional wind
                             us(:,:,:),                                & ! staggered zonal wind
                             vs(:,:,:),                                & ! staggered meridional wind
                             t(:,:,:),                                 & ! temperature
                             q(:,:,:),                                 & ! specific humidity
                             tt_lw(:,:,:),                             & ! tracer 1
                             tt_md(:,:,:),                             & ! tracer 2
                             tt_hi(:,:,:),                             & ! tracer 3
                             ttrmd(:,:,:),                             & ! tracer 4
!                            tt_un(:,:,:),                             & ! tracer
                             psi(:,:),                                 & ! streamfunction
                             ps(:,:),                                  & ! surface pressure
                             phis(:,:),                                & ! surface geopotential
                             output_1D(:),                             & ! 1D output stream 
                             output_1Ds(:)                               ! 1D output strean for staggered data

contains

!-----------------------------------------------------------------------
!  allocate data arrays
!-----------------------------------------------------------------------
  subroutine allocate_data_arrays
  
    implicit none

      if (model_version == 2) then                                       ! FV, staggered grid
        allocate (us(nlon,nlev,nslat),                                 & ! FV, staggered grid
                  vs(nlon,nlev,nlat),                                  &
                  output_1Ds(nstotal))
      else
        allocate (u(nlon,nlev,nlat),                                   & ! Eul or SLD, Gaussian grid
                  v(nlon,nlev,nlat))
      endif
                
      allocate (t (nlon,nlev,nlat),                                    &
                psi  (nlon,nlat),                                       &
                ps  (nlon,nlat),                                       &
                phis(nlon,nlat),                                       &
                output_1D(ntotal),                                     &
                q (nlon,nlev,nlat))
                
      if (choice_testcase == 11 .or. choice_testcase == 13) then
        allocate (tt_lw(nlon,nlev,nlat),                               &
                  tt_md(nlon,nlev,nlat),                               &
                  tt_hi(nlon,nlev,nlat),                               &
                  ttrmd(nlon,nlev,nlat))
!-----------------------------------
!       set tracers to zero, default
!-----------------------------------
        tt_lw = 0._r8
        tt_md = 0._r8
        tt_hi = 0._r8
        ttrmd = 0._r8
      endif
      
      if (choice_testcase == 12) then
        allocate (tt_lw(nlon,nlev,nlat))
!-----------------------------------
!       set tracers to zero, default
!-----------------------------------
        tt_lw = 0._r8
      endif

      if (choice_testcase == 400 .or. choice_testcase == 410 .or. choice_testcase == 411  .or. &
          choice_testcase == 412 .or. choice_testcase == 413 .or. &
          choice_testcase == 42  .or. choice_testcase == 43  .or. &
          choice_testcase == 42  .or. choice_testcase == 43  .or. &
          choice_testcase == 60  .or. choice_testcase == 61) then
        allocate (tt_lw(nlon,nlev,nlat),                               &
                  tt_md(nlon,nlev,nlat))
!-----------------------------------
!       set tracers to zero, default
!-----------------------------------
        tt_lw = 0._r8
        tt_md = 0._r8
      endif      

      q = 0._r8
                
   end subroutine allocate_data_arrays


!=======================================================================
!  subroutine write_to_file
!=======================================================================
   subroutine write_to_file

      implicit none
      integer :: i, j, ii
      character*250 line

!-----------------------------------------------------------------------
!     open input and output files
!-----------------------------------------------------------------------
      if (model_version == 1) then
!-----------------------------------------------------------------------
!       the resolutions correspond to the spectral resolutions T...
!-----------------------------------------------------------------------
        write (*,*) 'template.cam.'//trim(model)//'.T'//resol//'L'//c_level
        open (1, file='template.cam.'//trim(model)//'.T'//resol//'L'//c_level)
        open (2, file='initial_data.cam.'//trim(model)//'.T'//resol//'L'//c_level//  &
             '.'//trim(tc_nr)//'.'//trim(tc_name)//'.nc.data')
      else
!-----------------------------------------------------------------------
!       the resolution is arbitrary (FV only)
!-----------------------------------------------------------------------
        write (cnx,'(I4.4)') nlon
        write (cny,'(I4.4)') nlat
        open (1, file='template.cam.'//trim(model)//'.'//cny//'x'//cnx//'L'//c_level)
        open (2, file='initial_data.cam.'//trim(model)//'.'//cny//'x'//cnx//'L'//c_level// &
             '.'//trim(tc_nr)//'.'//trim(tc_name)//'.nc.data')
      endif        

!-----------------------------------------------------------------------
!     copy the netcdf template from a file
!-----------------------------------------------------------------------
   25 read (1,'(A)',end=30) line
      write (2,'(A)') trim(line)
      goto 25
   30 close (1)

!=======================================================================
!     PS surface pressure
!=======================================================================
      call concatenate_2D (ps, nlat, output_1D, nlayer)
      write (2,*) 'PS = '
      do i = 1, nlayer/4 - 1
         write (2,150) (output_1D(ii), ii=(i-1)*4+1,i*4)
      enddo
      write (2,155) (output_1D(ii), ii=nlayer-3,nlayer)
!=======================================================================
!     TS Surface temperature
!=======================================================================
      write (2,*) 'TS = '
      do i = 1, nlayer/8 - 1
         write (2,130) (ts1, j=1,8)
      enddo
      write (2,135) (ts1, j=1,8)
!=======================================================================
!     TS1 Surface temperature (level 1)
!=======================================================================
      write (2,*) 'TS1 = '
      do i = 1, nlayer/8 - 1
         write (2,130) (ts1, j=1,8)
      enddo
      write (2,135) (ts1, j=1,8)
!=======================================================================
!     TS2 Surface temperature (level 2)
!=======================================================================
      write (2,*) 'TS2 = '
      do i = 1, nlayer/8 - 1
         write (2,130) (ts2, j=1,8)
      enddo
      write (2,135) (ts2, j=1,8)
!=======================================================================
!     TS3 Surface temperature (level 3)
!=======================================================================
      write (2,*) 'TS3 = '
      do i = 1, nlayer/8 - 1
         write (2,130) (ts3, j=1,8)
      enddo
      write (2,135) (ts3, j=1,8)
!=======================================================================
!     TS4 Surface temperature (level 4)
!=======================================================================
      write (2,*) 'TS4 = '
      do i = 1, nlayer/8 - 1
         write (2,130) (ts4, j=1,8)
      enddo
      write (2,135) (ts4, j=1,8)
!=======================================================================
!     SNOWHICE Water equivalent snow depth
!=======================================================================
      write (2,*) 'SNOWHICE = '
      do i = 1, nlayer/8 - 1
         write (2,120) (snowh, j=1,8)
      enddo
      write (2,125) (snowh, j=1,8)
!=======================================================================
!     ORO ocean(1), land(2), sea ice(3) flag
!=======================================================================
!      write (2,*) 'ORO = '
!      do i = 1, nlayer/8 - 1
!         write (2,120) (oro, j=1,8)
!      enddo
!      write (2,125) (oro, j=1,8)
!=======================================================================
!     TSICE sea-ice model snow/ice surface temperature
!=======================================================================
      write (2,*) 'TSICE = '
      do i = 1, nlayer/8 - 1
         write (2,120) (tsice, j=1,8)
      enddo
      write (2,125) (tsice, j=1,8)
!=======================================================================
!     SGH orography standard deviation
!=======================================================================
      write (2,*) 'SGH = '
      do i = 1, nlayer/8 - 1
         write (2,120) (sgh, j=1,8)
      enddo
      write (2,125) (sgh, j=1,8)
!=======================================================================
!     LANDFRAC Fractional land in gridbox
!=======================================================================
      write (2,*) 'LANDFRAC = '
      do i = 1, nlayer/8 - 1
         write (2,120) (fland, j=1,8)
      enddo
      write (2,125) (fland, j=1,8)
!=======================================================================
!     LANDM land-ocean transition mask
!=======================================================================
      write (2,*) 'LANDM = '
      do i = 1, nlayer/8 - 1
         write (2,120) (landm, j=1,8)
      enddo
      write (2,125) (landm, j=1,8)
!=======================================================================
!     LANDM_COSLAT land-ocean transition mask
!=======================================================================
      write (2,*) 'LANDM_COSLAT = '
      do i = 1, nlayer/8 - 1
         write (2,120) (landm, j=1,8)
      enddo
      write (2,125) (landm, j=1,8)
!=======================================================================
!     SICTHK sea ice thickness
!=======================================================================
!      write (2,*) 'SICTHK = '
!      do i = 1, nlayer/8 - 1
!         write (2,120) (sicthk, j=1,8)
!      enddo
!      write (2,125) (sicthk, j=1,8)
!=======================================================================
!     CWAT
!=======================================================================
      write (2,*) 'CWAT = '
      do i = 1, ntotal/8 - 1
         write (2,120) (cwat, j=1,8)
      enddo
      write (2,125) (cwat, j=1,8)
!=======================================================================
!     Q specific humidity
!=======================================================================
      call concatenate (q, nlat, output_1D, ntotal)
      write (2,*) 'Q = '
      do i = 1, ntotal/4 - 1
         write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
      enddo
      write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)

   if (choice_testcase == 11 .or. choice_testcase == 12 .or. choice_testcase == 13 .or.  choice_testcase == 400 .or. &
        choice_testcase == 410 .or. choice_testcase == 411 .or. choice_testcase == 412 .or. choice_testcase == 413 .or. &
        choice_testcase == 42 .or. choice_testcase == 43 .or. choice_testcase == 60 .or. choice_testcase == 61) then

!=======================================================================
!     Tracers q1
!=======================================================================
      call concatenate (tt_lw, nlat, output_1D, ntotal)
      write (2,*) 'TT_LW = '
      j = 1
      DO i=1,ntotal-1
         IF (j.NE.4) THEN
            WRITE(2,200,ADVANCE='no') output_1D(i)
            j=j+1
         ELSE
            WRITE(2,200) output_1D(i)
            j = 1
         ENDIF
      ENDDO
      WRITE(2,201) output_1D(ntotal)
  200 format ('  ',(f22.15,', '))
  201 format ('  ',(f22.15,'; '))

    endif


  if (choice_testcase == 11 .or. choice_testcase == 13 .or. choice_testcase == 400 .or. &
        choice_testcase == 410 .or. choice_testcase == 411 .or. choice_testcase == 412 .or. choice_testcase == 413 .or. &
        choice_testcase == 42 .or. choice_testcase == 43 .or. choice_testcase == 60 .or. choice_testcase == 61) then

!=======================================================================
!     Tracers q2
!=======================================================================
      call concatenate (tt_md, nlat, output_1D, ntotal)
      write (2,*) 'TT_MD = '
      j = 1
      DO i=1,ntotal-1
         IF (j.NE.4) THEN
            WRITE(2,200,ADVANCE='no') output_1D(i)
            j=j+1
         ELSE
            WRITE(2,200) output_1D(i)
            j = 1
         ENDIF
      ENDDO
      WRITE(2,201) output_1D(ntotal)

  endif


   if (choice_testcase == 11 .or. choice_testcase == 13) then
!=======================================================================
!     Tracers q3
!=======================================================================
      call concatenate (tt_hi, nlat, output_1D, ntotal)
      write (2,*) 'TT_HI = '
      j = 1
      DO i=1,ntotal-1
         IF (j.NE.4) THEN
            WRITE(2,200,ADVANCE='no') output_1D(i)
            j=j+1
         ELSE
            WRITE(2,200) output_1D(i)
            j = 1
         ENDIF
      ENDDO
      WRITE(2,201) output_1D(ntotal)

!=======================================================================
!     Tracers q4
!=======================================================================
      call concatenate (ttrmd, nlat, output_1D, ntotal)
      write (2,*) 'TTRMD = '
      j = 1
      DO i=1,ntotal-1
         IF (j.NE.4) THEN
            WRITE(2,200,ADVANCE='no') output_1D(i)
            j=j+1
         ELSE
            WRITE(2,200) output_1D(i)
            j = 1
         ENDIF
      ENDDO
      WRITE(2,201) output_1D(ntotal)
!=======================================================================
!     offline tracer
!=======================================================================
      if (ABS(dble(ntotal)/dble(4)-ntotal/4)>0.001) then
         WRITE(*,*) 'ntotal not divisable by 4 - edit write in cam_variables'
         stop
      endif
    endif

      if (model_version.eq.1) then
!=======================================================================
!     U zonal velocity
!=======================================================================
        call concatenate (u, nlat, output_1D, ntotal)
        write (2,*) 'U = '
        do i = 1, ntotal/4 - 1
           write (2,140) (output_1D(ii), ii=(i-1)*4+1,i*4 )
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)

      else
!=======================================================================
!       US zonal velocity staggered
!=======================================================================
        call concatenate (us, nslat, output_1Ds, nstotal)
        write (2,*) 'US = '
        do i = 1, nstotal/4 - 1
           write (2,140) (output_1Ds(ii), ii=(i-1)*4+1,i*4 )
        enddo
        write (2,145) (output_1Ds(ii), ii=nstotal-3,nstotal)
      endif

      if (model_version.eq.1) then
!=======================================================================
!     V meridional velocity
!=======================================================================
        call concatenate (v, nlat, output_1D, ntotal)
        write (2,*) 'V = '
        do i = 1, ntotal/4 - 1
           write (2,140) (output_1D(ii), ii=(i-1)*4+1,i*4 )
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
!=======================================================================
!       VS zonal velocity staggered
!=======================================================================
        call concatenate (vs, nlat, output_1D, ntotal)
        write (2,*) 'VS = '
        do i = 1, ntotal/4 - 1
           write (2,140) (output_1D(ii), ii=(i-1)*4+1,i*4 )
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
!=======================================================================
!     T temperature
!=======================================================================
      call concatenate (t, nlat, output_1D, ntotal)
      write (2,*) 'T = '
      do i = 1, ntotal/4 - 1
         write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
      enddo
      write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
!=======================================================================
!     PHIS surface geopotential
!=======================================================================
      call concatenate_2D (phis, nlat, output_1D, nlayer)
      write (2,*) 'PHIS = '
      do i = 1, nlayer/4 - 1
         write (2,150) (output_1D(ii), ii=(i-1)*4+1,i*4)
      enddo
      write (2,155) (output_1D(ii), ii=nlayer-3,nlayer)
!=======================================================================
!     close file
!=======================================================================
      write (2,'(''}'')')
      close (2)

!-----------------------------------------------------------------------
!     output formats
!-----------------------------------------------------------------------
  120 format ('  ', 8(f2.0,', '))
  125 format ('  ', 7(f2.0,', '),f2.0,'; ',/)
  130 format ('  ', 8(f4.0,', '))
  135 format ('  ', 7(f4.0,', '),f4.0,'; ',/)
  140 format (('  ',4(e22.15,', ')))
  145 format ('  ', 3(e22.15,', '),e22.15,'; ',/)
  150 format (('  ',4(e22.15,', ')))
  155 format ('  ', 3(e22.15,', '),e22.15,'; ',/)
  160 format ('  ', 8(e6.1,', '))
  165 format ('  ', 7(e6.1,', '),e6.1,'; ',/) !6.1

  end subroutine write_to_file

!=======================================================================
! convert 3D fields into 1D output streams
!=======================================================================
  subroutine concatenate (var, nlimit, output, n)

      implicit none
      integer nlimit, n
      real*8 :: var (nlon, nlev, nlimit), &
                output (n)
      integer index
      integer i,j ,k

      index = 0
      do j = 1, nlimit
        do k = 1, nlev
           do i = 1, nlon
             index = index + 1
             output(index) = REAL(var(i,k,j))
           enddo
        enddo
      enddo
   end subroutine concatenate
!=======================================================================
! convert 3D fields into 1D output streams
!=======================================================================
  subroutine concatenate_2D (var, nlimit, output, n)

      implicit none
      integer nlimit, n
      real*8 :: var (nlon, nlimit), &
                output (n)
      integer index
      integer i, j

      index = 0
      do j = 1, nlimit
           do i = 1, nlon
             index = index + 1
             output(index) = real(var(i,j))
           enddo
      enddo
   end subroutine concatenate_2D

!-----------------------------------------------------------------------
!    deallocate arrays
!-----------------------------------------------------------------------
   subroutine deallocate_variables
     implicit none
     
     deallocate (q, t, phis, ps, output_1D)
     if (model_version == 2) then 
        deallocate (us, vs, output_1Ds)
     else 
        deallocate (u, v)
     endif
     if (choice_testcase == 11 .or. choice_testcase == 13) &
       deallocate (tt_lw, tt_md, tt_hi, ttrmd)
     if (choice_testcase == 12) &
       deallocate (tt_lw)
     if (choice_testcase == 410 .or. choice_testcase == 411 .or. choice_testcase == 412 .or. choice_testcase == 413 .or. &
         choice_testcase == 42  .or. choice_testcase == 43  .or. choice_testcase == 60  .or. choice_testcase == 61)      &
       deallocate (tt_lw, tt_md)
   end subroutine deallocate_variables

 end module cam_variables
