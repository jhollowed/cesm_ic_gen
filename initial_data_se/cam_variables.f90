 module cam_variables

      use cam_grid
      implicit none

!=======================================================================
!     CAM3 model variables (here with constant initial values)
!=======================================================================
      real*8, parameter ::  ts  = 302.15d0,                            &
                            ts1 = 271.36d0,                            &
                            ts2 = ts1,                                 &
                            ts3 = ts1,                                 &
                            ts4 = ts1,                                 &
                            snowh = 0.d0,                              &
                            oro   = 0.d0,                              &
                            tsice = 271.36d0,                          &
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
                             phis(:,:),                                & ! surface geopotential
                             ps(:,:),                                  & ! surface pressure
                             q(:,:,:),                                 & ! specific humidity
                             q1(:,:,:),                                 & ! tracer 1 
                             q2(:,:,:),                                 & ! tracer 2
                             q3(:,:,:),                                 & ! tracer 3 
                             q4(:,:,:),                                 & ! tracer 4
                             p(:,:,:),                                 & ! pressure
                             rho(:,:,:),                               & ! density
                             output_1D(:),                             & ! 1D output stream
                             output_1Ds(:)                               ! 1D output strean for staggered data

contains

!-----------------------------------------------------------------------
!  allocate data arrays
!-----------------------------------------------------------------------
   subroutine allocate_data_arrays
      use cam_grid
      implicit none
      if (model_version.eq.1 .or.model_version.eq.2)                   & ! EUL, SLD, & FV
      allocate (u(nlon,nlev,nlat),                                     &
                v(nlon,nlev,nlat),                                     &
                t(nlon,nlev,nlat),                                     &
                q(nlon,nlev,nlat),                                     &
                q1(nlon,nlev,nlat),                                     &
                q2(nlon,nlev,nlat),                                     &
                q3(nlon,nlev,nlat),                                     &
                q4(nlon,nlev,nlat),                                     &
                ps(nlon,nlat),                                         &
                phis(nlon,nlat),                                       &
                p(nlon,nlev,nlat),                                     &
                rho(nlon,nlev,nlat),                                   &
                output_1D(ntotal))

 
      if (model_version.eq.2)                                          & ! FV
        allocate (us(nlon,nlev,nslat),                                 & 
                  vs(nlon,nlev,nlat),                                  &
                  output_1Ds(nstotal))

      if (model_version.eq.3)                                          & ! HOMME
      allocate (u(nlev,ncol,1),                                          &
                v(nlev,ncol,1),                                          &
                t(nlev,ncol,1),                                          &
                q(nlev,ncol,1),                                          &
                q1(nlev,ncol,1),                                          &
                q2(nlev,ncol,1),                                          &
                q3(nlev,ncol,1),                                          &
                q4(nlev,ncol,1),                                          &
                ps(ncol,1),                                              &
                phis(ncol,1),                                            &
                p(nlev,ncol,1),                                          &
                rho(nlev,ncol,1),                                          &
                output_1D(ntotal))

   end subroutine allocate_data_arrays


!=======================================================================
!  subroutine write_to_file
!=======================================================================
   subroutine write_to_file

      use cam_grid

      implicit none
      integer :: i, j, ii
 
!-----------------------------------------------------------------------
!     open input and output files
!-----------------------------------------------------------------------
      if (model_version == 1) then
!-----------------------------------------------------------------------
!       the resolutions correspond to the spectral resolutions T...
!-----------------------------------------------------------------------
        open (1, file='template.T'//resol//'.L'//c_level//'.cam1.'//        &
             trim(model))
        open (2, file='initial_data.cam.T'//resol//'.L'//c_level//'.'//  &
             trim(model)//'.nc.data')
      else if (model_version == 2) then
!-----------------------------------------------------------------------
!       the resolution is arbitrary (FV only)
!-----------------------------------------------------------------------
        write (cnx,'(I4.4)') nlon
        write (cny,'(I4.4)') nlat
        open (1, file='template.'//cnx//'x'//cny//'.L'//   &
             c_level//'.cam1.'//trim(model))
        open (2, file='initial_data.cam.'//cnx//'x'//cny//'.L'//c_level//'.'// &
             trim(model)//'.nc.data')
      else
        write (cn,'(I4.4)') ne
        open (1, file='template.'//cn//'.L'//   &
             c_level//'.cam1.'//trim(model))
        open (2, file='initial_data.cam.'//trim(model)//'.ne'//cn//'L'//c_level//'.'// &
              trim(tc_nr)//'.'//trim(tc_name)//'.nc.data') 
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
      if (model_version.eq.1 .or.model_version.eq.2) then
      call concatenate_2D (ps, nlat, output_1D, nlayer)
      write (2,*) 'PS = '
      do i = 1, nlayer/4 - 1
         write (2,150) (output_1D(ii), ii=(i-1)*4+1,i*4)
      enddo 
      write (2,155) (output_1D(ii), ii=nlayer-3,nlayer)
      else
      write (2,*) 'PS = '
      do i = 1, (nlayer-2)/4
         write (2,150) (ps(ii,1), ii=(i-1)*4+1,i*4)
      enddo
      write (2,170) (ps(ii,1), ii=nlayer-1,nlayer)
      end if
 
   !   write (2,*) 'PS = '
   !   do i = 1, nlayer/8 - 1
   !      write (2,160) (p0, j=1,8)
   !   enddo
   !   write (2,165) (p0, j=1,8)
!=======================================================================
!     TS Surface temperature
!=======================================================================
      write (2,*) 'TS = '
      do i = 1, nlayer/8 - 1
         write (2,130) (ts, j=1,8)
      enddo
      write (2,135) (ts, j=1,8)
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
         write (2,130) (tsice, j=1,8)
      enddo
      write (2,135) (tsice, j=1,8)
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
      if (model_version.eq.1 .or.model_version.eq.2) then
        call concatenate (q, nlat, output_1D, ntotal)
        write (2,*) 'Q = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
        call concatenate_2DHOMME (q, output_1D, ntotal)
        write (2,*) 'Q = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
!=======================================================================
!     Q1 Tracer (TT_LW)
!=======================================================================
      if (choice_testcase .ne. 200) then
      if (model_version.eq.1 .or.model_version.eq.2) then
        call concatenate (q1, nlat, output_1D, ntotal)
        write (2,*) 'TT_LW = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
        call concatenate_2DHOMME (q1, output_1D, ntotal)
        write (2,*) 'TT_LW = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
      endif
!=======================================================================
!     Q2 Tracer (TT_MD)
!=======================================================================
      if (choice_testcase .ne. 200) then
      if (model_version.eq.1 .or.model_version.eq.2) then
        call concatenate (q2, nlat, output_1D, ntotal)
        write (2,*) 'TT_MD = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
        call concatenate_2DHOMME (q2, output_1D, ntotal)
        write (2,*) 'TT_MD = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
      endif
!=======================================================================
!     Q3 Tracer (TT_HI)
!=======================================================================
      if (choice_testcase .eq. 11) then
      if (model_version.eq.1 .or.model_version.eq.2) then
        call concatenate (q3, nlat, output_1D, ntotal)
        write (2,*) 'TT_HI = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
        call concatenate_2DHOMME (q3, output_1D, ntotal)
        write (2,*) 'TT_HI = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
      endif
!=======================================================================
!     Q4 Tracer (TTRMD)
!=======================================================================
      if (choice_testcase .eq. 11) then
      if (model_version.eq.1 .or.model_version.eq.2) then
        call concatenate (q4, nlat, output_1D, ntotal)
        write (2,*) 'TTRMD = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
        call concatenate_2DHOMME (q4, output_1D, ntotal)
        write (2,*) 'TTRMD = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
      endif
!=======================================================================
!     U zonal velocity
!=======================================================================
      if (model_version.eq.1) then
        call concatenate (u, nlat, output_1D, ntotal)
        write (2,*) 'U = '
        do i = 1, ntotal/4 - 1
           write (2,140) (output_1D(ii), ii=(i-1)*4+1,i*4 )
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else if (model_version.eq.2) then
!=======================================================================
!       US zonal velocity staggered
!=======================================================================
        call concatenate (us, nslat, output_1Ds, nstotal)
        write (2,*) 'US = '
        do i = 1, nstotal/4 - 1
           write (2,140) (output_1Ds(ii), ii=(i-1)*4+1,i*4 )
        enddo
        write (2,145) (output_1Ds(ii), ii=nstotal-3,nstotal)
      else
        call concatenate_2DHOMME (u, output_1D, ntotal)
        write (2,*) 'U = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
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
      else if (model_version.eq.2) then
!=======================================================================
!       VS zonal velocity staggered
!=======================================================================
        call concatenate (vs, nlat, output_1D, ntotal)
        write (2,*) 'VS = '
        do i = 1, ntotal/4 - 1
           write (2,140) (output_1D(ii), ii=(i-1)*4+1,i*4 )
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
!=======================================================================
!     V meridional velocity (HOMME)
!=======================================================================
        call concatenate_2DHOMME (v, output_1D, ntotal)
        write (2,*) 'V = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
!=======================================================================
!     T temperature
!=======================================================================
      if (model_version.eq.1 .or.model_version.eq.2) then
        call concatenate (t, nlat, output_1D, ntotal)
        write (2,*) 'T = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      else
        call concatenate_2DHOMME (t, output_1D, ntotal)
        write (2,*) 'T = '
        do i = 1, ntotal/4 - 1
           write(2,140) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,145) (output_1D(ii), ii=ntotal-3,ntotal)
      endif
!=======================================================================
!     PHIS surface geopotential
!=======================================================================
      if (model_version.eq.1 .or.model_version.eq.2) then
        call concatenate_2D (phis, nlat, output_1D, nlayer)
        write (2,*) 'PHIS = '
        do i = 1, nlayer/4 - 1
           write (2,150) (output_1D(ii), ii=(i-1)*4+1,i*4)
        enddo
        write (2,155) (output_1D(ii), ii=nlayer-3,nlayer)
      else
      write (2,*) 'PHIS = '
      do i = 1, (nlayer-2)/4 
         write (2,150) (phis(ii,1), ii=(i-1)*4+1,i*4)
      enddo
      write (2,170) (phis(ii,1), ii=nlayer-1,nlayer)
      endif
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
  165 format ('  ', 7(e6.1,', '),e6.1,'; ',/)
  170 format ('  ', 1(e22.15,', '),e22.15,'; ',/)

  end subroutine write_to_file

!=======================================================================
! convert 3D fields into 1D output streams
!=======================================================================
  subroutine concatenate (var, nlimit, output, n)
      use cam_grid, only: nlev, nlon
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
             output(index) = var(i,k,j)
           enddo
        enddo
      enddo
   end subroutine concatenate
!=======================================================================
! convert 2D fields into 1D output streams for EUL, SLD, & FV
!=======================================================================
  subroutine concatenate_2D (var, nlimit, output, n)
      use cam_grid, only: nlon
      implicit none
      integer nlimit, n
      real*8 :: var (nlon, nlimit), &
                output (n)
      integer index
      integer i,j ,k

      index = 0
      do j = 1, nlimit
           do i = 1, nlon
             index = index + 1
             output(index) = var(i,j)
           enddo
      enddo
   end subroutine concatenate_2D
!=======================================================================
! convert 2D fields into 1D output streams for HOMME
!=======================================================================
  subroutine concatenate_2DHOMME (var, output, n)
      use cam_grid, only: ncol, nlev
      implicit none
      integer n
      real*8 :: var (nlev,ncol,1), &
                output (n)
      integer index
      integer j ,k

      index = 0
      do k = 1, nlev
        do j = 1, ncol
           index = index + 1
           output(index) = var(k,j,1)
        enddo
      enddo
   end subroutine concatenate_2DHOMME

!-----------------------------------------------------------------------
!    deallocate arrays
!-----------------------------------------------------------------------
   subroutine deallocate_variables
     implicit none
     integer :: model_version     

     deallocate (u, v, t, phis, ps, q, output_1D)
     if (model_version.eq.2) deallocate (us, vs, output_1Ds)

   end subroutine deallocate_variables
 
 end module cam_variables
