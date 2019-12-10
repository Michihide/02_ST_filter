! Filtereing data of ST by Gaussian filter
program main
  use filter

  call input_prmtr
  call read_fln
  call read_array_and_xy_data

  do k = 1, kend, 1
     call read_vtk(5, 1, wtr_mdn) ! ---- deviation_of_water_level ----
     call read_vtk(5, 2, wl_mdn)  ! ---- water_level ----
     call read_vtk(5, 3, bed_mdn) ! ---- deviation_of_bed_level ----
     call read_vtk(5, 4, bel_mdn) ! ---- bed_level ----     
     call read_vtk(5, 5, dep_mdn) ! ---- depth_of_water ----

     call intrplltn_mssng_pnts(iend, jend, wtr_mdn)
     call intrplltn_mssng_pnts(iend, jend, wl_mdn)     
     call intrplltn_mssng_pnts(iend, jend, bed_mdn)
     call intrplltn_mssng_pnts(iend, jend, bel_mdn)     
     
!!$     call input_scaling_structure
!!$     call interpolate_around_scaling_structure

     call gaussian_filter(iend, jend, wtr_mdn, wtr_gaus, smthng)
     call gaussian_filter(iend, jend, wl_mdn , wl_gaus , smthng)
     call gaussian_filter(iend, jend, bed_mdn, bed_gaus, smthng)
     call gaussian_filter(iend, jend, bel_mdn, bel_gaus, smthng)

     call add_grad(iend, jend, x, y, wl_gaus , ib, wl_gaus)
     call add_grad(iend, jend, x, y, bel_gaus, ib, bel_gaus)     

     call cal_grad_xy(iend, jend, x, y, wl_mdn  , dwldx_mdn,  dwldy_mdn) 
     call cal_grad_xy(iend, jend, x, y, wl_gaus , dwldx_gaus, dwldy_gaus)
     call cal_grad_xy(iend, jend, x, y, bed_mdn, beddx_mdn, beddy_mdn)
     call cal_grad_xy(iend, jend, x, y, bed_gaus, beddx_gaus, beddy_gaus)

     wtr_dpth(:,:) = wl_gaus(:,:) - bel_gaus(:,:)
     
     call cal_tau(iend, jend, x, y, wtr_dpth, dwldx_gaus, dwldy_gaus, taus, tausx, tausy)

     call output_schalar_vtk
     call output_vector_vtk

!!$     do j = 1, jend
!!$        call drw_cntr
!!$     enddo

  enddo

end program main



! Input parameter
subroutine input_prmtr
  use filter  

  open(100,file = '../Target_file.txt')
  read(100,*)case
  close(100)
  
  ib       = 0
  intrvl_i = 1
  intrvl_j = 1
  smthng   = 1
    
  call system('mkdir ../00_input_vtk/'//trim(case)//'/vtk_schalar_gaussian')
  call system('mkdir ../00_input_vtk/'//trim(case)//'/vtk_vector_gaussian')  

end subroutine input_prmtr



! Read file name
subroutine read_fln
  use filter

  call system('ls -1 ../00_input_vtk/'//trim(case)//'/vtk_schalar_median/ > ./file.txt')

  open(100,file = 'file.txt')  

  k = 0
  do
     k = k + 1
     read(100, *, end = 998)
  enddo
998 write(*,*) 'Number of file =', k - 1
  kend = k - 1
  close(100)
  
  allocate(fln(kend))
  
  open(100,file = 'file.txt')  
  read(100,*)(fln(k), k = 1, kend) 
  close(100)

end subroutine read_fln



! Read array
subroutine read_array_and_xy_data
  use filter
  
  open(100,file = '../00_input_vtk/'//trim(case)//'/vtk_schalar_median/'//trim(fln(1)))
  read(100,'(a)')(chr, i = 1, 4)
  read(100,*)chr, iend, jend
  call set_array
  read(100,*)  
  do j = 1, jend
     do i = 1, iend
        read(100,*)x(i,j), y(i,j)
     enddo
  enddo
  close(100)
  
end subroutine read_array_and_xy_data



! Set array
subroutine set_array
  use filter
  lend = iend*jend
  
  allocate(x(iend,jend)         , y(iend,jend))
  allocate(wtr_mdn(iend,jend)   , wl_mdn(iend,jend) , dep_mdn(iend,jend))
  allocate(wtr_gaus(iend,jend)  , wl_gaus(iend,jend), dev_wl_gaus(iend,jend))
  allocate(bed_mdn(iend,jend)   , bel_mdn(iend,jend))
  allocate(bed_gaus(iend,jend)  , bel_gaus(iend,jend))
  allocate(dwldx_mdn(iend,jend) , dwldy_mdn(iend,jend))
  allocate(beddx_mdn(iend,jend) , beddy_mdn(iend,jend))
  
  allocate(dwldx_gaus(iend,jend), dwldy_gaus(iend,jend))
  allocate(beldx_gaus(iend,jend), beldy_gaus(iend,jend))
  allocate(beddx_gaus(iend,jend), beddy_gaus(iend,jend))
  
  allocate(tausx(iend,jend), tausy(iend,jend), taus(iend,jend), wtr_dpth(iend,jend))
  allocate(dis(lend), x_itp(lend), y_itp(lend), bel_itp(lend), bed_itp(lend), wl_itp(lend))

end subroutine set_array



! Read vtk
subroutine read_vtk(nend, nnm, f)
  use filter
  integer            :: n  
  integer,intent(in) :: nend, nnm
  real*8,intent(out) :: f(iend,jend)  

  open(100,file = '../00_input_vtk/'//trim(case)//'/vtk_schalar_median/'//trim(fln(k)))  

  read(100,'(a)')(chr, i = 1, 6)
  do j = 1, jend
     do i = 1, iend
        read(100,*)
     enddo
  enddo

  read(100,'(a)')(chr, i = 1, 3)  

  do n = 1, nend
     read(100,*)   
     do j = 1, jend
        do i = 1, iend
           if(n.eq.nnm)then
              read(100,*)f(i,j)
           elseif(n.ne.nnm)then              
              read(100,*)
           endif
        end do
     end do
  end do
  close(100)

end subroutine read_vtk



subroutine input_scaling_structure
  use filter
  do j = 1, jend
     do i = 1, iend
        if(wl_mdn(i,j) .lt. 0.0 .or. bel_mdn(i,j) .lt. 0.0 )then
           wl_mdn  (i,j) = 0.0
           wtr_mdn (i,j) = 0.0
           bed_mdn (i,j) = 0.1
           bel_mdn (i,j) = 0.1
           wl_gaus (i,j) = 0.0
           wtr_gaus(i,j) = 0.0
           bed_gaus(i,j) = 0.1
           bel_gaus(i,j) = 0.1
        endif
     enddo
  enddo
end subroutine input_scaling_structure



subroutine interpolate_around_scaling_structure
  use filter
  mend = 5
  
  ! ---- interpolate_around_scaling_structure ----
  do j = 1, jend
     do i = 1, iend
        if(i.ne.1)then
           if(bel_mdn(i,j).eq.0.1 .and. bel_mdn(i-1,j).ne.0.1)then
              bed_gaus(i-mend:i-1,j) = bed_gaus(i-(mend+1),j)
              bel_gaus(i-mend:i-1,j) = bel_gaus(i-(mend+1),j)
           endif
        endif
        if(i.ne.iend)then
           if(bel_mdn(i,j).eq.0.1 .and. bel_mdn(i+1,j).ne.0.1)then
              bed_gaus(i+1:i+mend,j) = bed_gaus(i+(mend+1),j)
              bel_gaus(i+1:i+mend,j) = bel_gaus(i+(mend+1),j)
           endif
        endif
     enddo
  enddo
  
  do i = 1, iend
     do j = 1, jend
        if(j.ne.jend)then
           if(bel_mdn(i,j).eq.0.1 .and. bel_mdn(i,j+1).ne.0.1)then
              bed_gaus(i,j+1:j+mend) = bed_gaus(i,j+(mend+1))
              bel_gaus(i,j+1:j+mend) = bel_gaus(i,j+(mend+1))
           endif
        endif
        if(j.ne.1)then
           if(bel_mdn(i,j).eq.0.1 .and. bel_mdn(i,j-1).ne.0.1)then
              bed_gaus(i,j-mend:j-1) = bed_gaus(i,j-(mend+1))
              bel_gaus(i,j-mend:j-1) = bel_gaus(i,j-(mend+1))
           endif
        endif
     enddo
  enddo
  
  do i = 1, iend
     do j = 1, jend
        if(bel_gaus(i,j) .le. 0.02)then
           bed_gaus(i,j) = 0.05
           bel_gaus(i,j) = 0.05
        endif
     enddo
  enddo
  
  ! ! ! ---- moving_average_interpolate_z_around_scaling_structure ----
  ! do j = 1, jend
  !    do i = 1, iend
  !       if(i.ne.1)then
  !          if(bel_gaus(i,j).eq.0.1 .and. bel_gaus(i-1,j).ne.0.1)then
  !             do m = mend, 1, -1
  !                bel_gaus(i-m,j) = sum(bel_gaus((i-m-mend):(i-m),j)) / mend
  !             enddo
  !          endif
  !       endif
  !       if(i.ne.iend)then
  !          if(bel_gaus(i,j).eq.0.1 .and. bel_gaus(i+1,j).ne.0.1)then
  !             do m = mend, 1, -1
  !                bel_gaus(i+m,j) = sum(bel_gaus((i+m):(i+m+mend),j)) / mend
  !             enddo
  !          endif
  !       endif
  !    enddo
  ! enddo

  ! do j = 1, jend
  !    do i = 1, iend
  !       if(i.ne.1)then
  !          if(bel_mdn(i,j).eq.0.1 .and. bel_mdn(i-1,j).ne.0.1)then
  !             do m = mend, 1, -1
  !                bed_gaus(i-m,j) = sum(bed_gaus(i-m-mend:i-m,j)) / mend
  !                bel_gaus(i-m,j) = sum(bel_gaus(i-m-mend:i-m,j)) / mend
  !                write(*,*)i-m-mend,i-m,i,m
  !             enddo
  !          endif
  !       endif
  !       if(i.ne.iend)then
  !          if(bel_mdn(i,j).eq.0.1 .and. bel_mdn(i+1,j).ne.0.1)then
  !             do m = 1, mend
  !                bed_gaus(i+m,j) = bed_gaus(i+(m+1),j)
  !                bel_gaus(i+m,j) = bel_gaus(i+(m+1),j)
  !             enddo
  !          endif
  !       endif
  !    enddo
  ! enddo

end subroutine interpolate_around_scaling_structure



! subroutine interpolate_z
!   use filter
!   l = 0
!   do i = 1, iend
!      do j = 1, jend
!         l = l + 1
!         x_itp(l)   = x(i,j)
!         y_itp(l)   = y(i,j)
!         wl_itp(l)  = wl_mdn(i,j)
!         bel_itp(l) = bel_mdn(i,j)
!         bed_itp(l) = bed_mdn(i,j)
!      enddo
!   enddo

!   do i = 1, iend
!      do j = 1, jend
!         do l = 1, lend
!            if(bed_mdn(i,j) .lt. 0.0d0)then
!               dis(l) = sqrt( (x_itp(l) - x(i,j))**2 + (y_itp(l) - y(i,j))**2 )
!            endif
!         enddo
!         if(bed_mdn(i,j) .lt. 0.0d0)then
!            wl_mdn(i,j)  = wl_itp( minloc(dis(:),1) )
!            bed_mdn(i,j) = bed_itp( minloc(dis(:),1) )
!            bel_mdn(i,j) = bel_itp( minloc(dis(:),1) )
!         endif
!      enddo
!   enddo
    
! end subroutine interpolate_z



! Add gradient
subroutine add_grad(iend, jend, x, y, z, ib, zg)
  implicit none   
  integer            :: i, j
  integer,intent(in) :: iend, jend
  real*8,intent(in)  :: x(iend,jend), y(iend,jend), z(iend,jend), ib
  real*8,intent(out) :: zg(iend,jend)

  do j = 1, jend
     do i = 1, iend
        zg(i,j) = z(i,j) - (x(i,j)*ib)
     enddo
  enddo

end subroutine add_grad



! Interpolation of missing point
subroutine intrplltn_mssng_pnts(iend, jend, f)
  implicit none
  integer              :: i, j, n
  integer,parameter    :: nend=1
  integer,intent(in)   :: iend, jend
  real*8,intent(inout) :: f(iend,jend)  
  real*8               :: ff(iend,jend)

  ff = f
  do n = 1, nend
     do j = 2, jend-1
        do i = 2, iend-1
           if(f(i,j).le.0.1)then
              ff(i,j) = (sum(f(i-1:i+1,j-1:j+1)) - f(i,j))/8.d0
           endif
        enddo
     enddo
     
     do i = 2, iend-1
        if(f(i,1).le.0.1)then
           ff(i,1) = (sum(f(i-1:i+1,1:2))-f(i,1))/5.
        elseif(f(i,jend).le.0.1)then
           ff(i,jend) = (sum(f(i-1:i+1,jend-1:jend))-f(i,jend))/5.
        endif
     enddo

     do j = 2, jend-1
        if(f(1,j).le.0.1)then        
           ff(1,j) = (sum(f(1:2,j-1:j+1))-f(1,j))/5.
        elseif(f(iend,j).le.0.1)then        
           ff(iend,j) = (sum(f(iend-1:iend,j-1:j+1))-f(iend,j))/5.
        endif
     enddo

     if(f(1,1).le.0.1)then       
        ff(1,1) = (sum(f(1:2,1:2))-f(1,1))/3.
     elseif(f(iend,1).le.0.)then       
        ff(iend,1) = (sum(f(iend-1:iend,1:2))-f(iend,1))/3.
     elseif(f(iend,jend).le.0.1)then       
        ff(iend,jend) = (sum(f(iend-1:iend,jend-1:jend))-f(iend,jend))/3.
     elseif(f(1,jend).le.0.1)then       
        ff(1,jend) = (sum(f(1:2,jend-1:jend))-f(1,jend))/3.
     endif
     
  enddo
  f = ff
  
end subroutine intrplltn_mssng_pnts



! Gaussian filter
subroutine gaussian_filter(iend, jend, f, gf, nend)
  implicit none
  integer            :: i, j, n
  integer,intent(in) :: iend, jend, nend
  real*8,intent(in)  :: f(iend,jend)
  real*8,intent(out) :: gf(iend,jend)
  real*8             :: ff(iend,jend)  

  ff = f
  gf = f

  do n = 1, nend
     do j = 2, jend-1
        do i = 2, iend-1
           if(sum(f(i-1:i+1,j-1:j+1)).ne.0.)then
              gf(i,j) = 1./16.*(4.*ff(i,j)+2.*ff(i-1,j)+2.*ff(i+1,j)+2.*ff(i,j-1)+2.*ff(i,j+1)&
                   +ff(i-1,j-1)+ff(i-1,j+1)+ff(i+1,j-1)+ff(i+1,j+1))
           else
              gf(i,j) = ff(i,j)
           endif
        enddo
     enddo

     do i = 2, iend-1
        gf(i,1)    = sum(ff(i-1:i+1,1:2))/6.
        gf(i,jend) = sum(ff(i-1:i+1,jend-1:jend))/6.        
     enddo

     do j = 2, jend-1
        gf(1,j)    = sum(ff(1:2,j-1:j+1))/6.
        gf(iend,j) = sum(ff(iend-1:iend,j-1:j+1))/6.        
     enddo

     gf(1,1)       = sum(ff(1:2,1:2))/4.
     gf(iend,1)    = sum(ff(iend-1:iend,1:2))/4.
     gf(iend,jend) = sum(ff(iend-1:iend,jend-1:jend))/4.
     gf(1,jend)    = sum(ff(1:2,jend-1:jend))/4.
     
     ff = gf
  enddo
  
end subroutine gaussian_filter



! Calculate of gradient
subroutine cal_grad_xy(iend, jend, x, y, f, dfdx, dfdy)
  implicit none   
  integer            :: i, j
  integer,intent(in) :: iend, jend
  real*8,intent(in)  :: x(1:iend,1:jend), y(1:iend,1:jend), f(1:iend,1:jend)
  real*8,intent(out) :: dfdx(1:iend,1:jend), dfdy(1:iend,1:jend)
  real*8             :: dfdx_ave, dfdy_ave
  real*8             :: dfdx_lim, dfdy_lim

  do j = 1, jend
     do i = 1, iend
        if(i .eq. 1)then
           dfdx(i,j) = (f(i,j) - f(i+1,j)) / abs((x(i,j) - x(i+1,j)))
        elseif(i .eq. iend)then            
           dfdx(i,j) = (f(i-1,j) - f(i,j)) / abs((x(i-1,j) - x(i,j)))
        else
           dfdx(i,j) = (f(i-1,j) - f(i+1,j)) / abs((x(i-1,j) - x(i+1,j)))
        endif
     enddo
  enddo

  do i = 1, iend
     do j = 1, jend
        if(j .eq. 1)then
           dfdy(i,j) = (f(i,j) - f(i,j+1)) / abs(y(i,j) - y(i,j+1))
        elseif(j .eq. jend)then            
           dfdy(i,j) = (f(i,j-1) - f(i,j)) / abs(y(i,j-1) - y(i,j))
        else
           dfdy(i,j) = (f(i,j-1) - f(i,j+1)) / abs(y(i,j-1) - y(i,j+1))
        endif
     enddo
  enddo

  dfdx_ave = sum(dfdx(1:iend,1:jend)) / (iend*jend)
  dfdy_ave = sum(dfdx(1:iend,1:jend)) / (iend*jend)

  dfdx_lim = 0.01
  dfdy_lim = 0.01
  
  do i = 1, iend
     do j = 1, jend
        if(dfdx(i,j).lt.-1.*dfdx_lim .or. dfdx(i,j).gt.dfdx_lim)then
           dfdx(i,j) = dfdx(i,j) * 0.01
        endif
        if(dfdy(i,j).lt.-1.*dfdy_lim .or. dfdy(i,j).gt.dfdy_lim)then
           dfdy(i,j) = dfdy(i,j) * 0.01
        endif
     enddo
  enddo
  
end subroutine cal_grad_xy



! Output vtk data
subroutine cal_tau(iend, jend, x, y, dep, dfdx, dfdy, tau, taux, tauy)
  implicit none   
  integer            :: i, j
  real*8 :: d,s
  integer,intent(in) :: iend, jend
  real*8,intent(in)  :: x(1:iend,1:jend), y(1:iend,1:jend)
  real*8,intent(in)  :: dep(1:iend,1:jend), dfdx(1:iend,1:jend), dfdy(1:iend,1:jend)
  real*8,intent(out) :: tau(1:iend,1:jend), taux(1:iend,1:jend), tauy(1:iend,1:jend)
  s = 1.65
  d = 0.00076

  do i = 1, iend
     do j = 1, jend
        taux(i,j) = (dep(i,j) * dfdx(i,j) / (s * d)) / 0.034
        tauy(i,j) = (dep(i,j) * dfdy(i,j) / (s * d)) / 0.034
        tau(i,j)  = sqrt(taux(i,j)**2 + tauy(i,j)**2) / 0.034
     enddo
  enddo
  
end subroutine cal_tau



! Output vtk data
subroutine output_schalar_vtk
  use filter  

  ic = 0
  do i = 1, iend, intrvl_i
!!$     do i = 1, iend-120, intrvl_i
     ic = ic + 1
  enddo
  
  jc = 0
  do j = 1, jend, intrvl_j
     jc = jc + 1
  enddo

  ijc = ic * jc

  open( 100,file='../00_input_vtk/'//trim(case)//'/vtk_schalar_gaussian/'//trim(fln(k)))
  write(100,'(a)') '# vtk DataFile Version 3.1'
  write(100,'(a,i4.4,a)') 'out.vtk'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET STRUCTURED_GRID'
  write(100,'(a,5x,i4,5x,i4,5x,i4)')'DIMENSIONS', ic, jc, 1
  write(100,'(a,5x,i8,5x,a)')'POINTS', ijc,'float'

  do j = 1, jend, intrvl_j          
     do i = 1, iend, intrvl_i
        write(100,*)x(i,j), y(i,j), 0.d0
     enddo
  enddo

  write(100,*)
  write(100,'(a,5x,i8)')'POINT_DATA', ijc
  write(100,'(a,5x,i8)')'FIELD FieldData', 5
 
  call out_scalar(wtr_gaus, 'deviation_of_water_level')
  call out_scalar(wl_gaus,  'water_level')
  call out_scalar(bed_gaus, 'deviation_of_bed_level')
  call out_scalar(bel_gaus, 'bed_level')
  call out_scalar(wtr_gaus(:,:) - bed_gaus(:,:), 'depth_of_water')
  
!!$  call read_vtk(5, 1, wtr_mdn) ! ---- deviation_of_water_level ----
!!$  call read_vtk(5, 2, wl_mdn)  ! ---- water_level ----
!!$  call read_vtk(5, 3, bed_mdn) ! ---- deviation_of_bed_level ----
!!$  call read_vtk(5, 4, bel_mdn) ! ---- deviation_of_water_level ----     
!!$  call read_vtk(5, 5, dep_mdn) ! ---- depth_of_water ----

  close(100)
end subroutine output_schalar_vtk



! Output scalar data
subroutine out_scalar(f, name)
  use filter  
  real*8,intent(in)       :: f(iend, jend)
  character(*),intent(in) :: name

  write(100,'(a,5x,i4,5x,i8,5x,a)')name, 1, ijc, 'double'

  do j = 1, jend , intrvl_j          
     do i = 1, iend, intrvl_i
        write(100,*)f(i,j)
     end do
  end do

end subroutine out_scalar




! Output vector data to vtk
subroutine output_vector_vtk
  use filter  

  ic = 0
  do i = 1, iend, intrvl_i
     ic = ic + 1
  enddo

  jc = 0
  do j = 1, jend, intrvl_j
     jc = jc + 1
  enddo

  ijc = ic * jc

  open( 100,file='../00_input_vtk/'//trim(case)//'/vtk_vector_gaussian/'//trim(fln(k)))  
  write(100,'(a)') '# vtk DataFile Version 3.1'
  write(100,'(a,i4.4,a)') 'out.vtk'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET STRUCTURED_GRID'
  write(100,'(a,5x,i4,5x,i4,5x,i4)')'DIMENSIONS', ic, jc, 1
  write(100,'(a,5x,i8,5x,a)')'POINTS', ijc,'float'

  do j = 1, jend, intrvl_j
     do i = 1, iend, intrvl_i
        write(100,*)x(i,j), y(i,j), 0.d0
     enddo
  enddo

  write(100,*)
  write(100,'(a,5x,i8)')'POINT_DATA', ijc
  write(100,'(a,5x,i8)')'FIELD FieldData', 1
 
  call out_vector(dwldx_gaus, dwldy_gaus, 'gradient_of_water_level')
  ! call out_vector(dbldx_gaus, dbldy_gaus, 'gradient_of_bed_level')  
  
  close(100)

end subroutine output_vector_vtk



! Output_vector
subroutine out_vector(fx, fy, name)
  use filter
  real*8,intent(in)       :: fx(iend,jend), fy(iend,jend)
  character(*),intent(in) :: name

  write(100,'(a,5x,i4,5x,i8,5x,a)')name, 3, ijc, 'double'
  ! do j = 1, jend, intrvl_j
  !    do i = 1, iend, intrvl_i
  !       if( fx(i,j) .lt. 0 )then
  !          write(100,*)0, fy(i,j), 0
  !       else
  !          write(100,*)fx(i,j), fy(i,j), 0
  !       endif
  !    end do
  ! end do
  
  do j = 1, jend, intrvl_j
     do i = 1, iend, intrvl_i
        write(100,*)fx(i,j), fy(i,j), 0
     end do
  end do

end subroutine out_vector
