! Structuring and Filtering data
program main
  use filter

  call input_prmtr
  call read_fln
  call structuring

  do k = 1, kend
     call set_txt_file
     call read_arry(wend, pss_wtr)
     call read_arry(bend, pss_bed)
     call set_arry
     call read_txt(wend, xw, yw, zw, pss_wtr)
     call read_txt(bend, xb, yb, zb, pss_bed)
     call count_arry(wend, wend_trm, xw, yw)
     call count_arry(bend, bend_trm, xb, yb)
     call set_arry2     
     call trim_txt(wend, xw, yw, zw, xw_trm, yw_trm, zw_trm)
     call trim_txt(bend, xb, yb, zb, xb_trm, yb_trm, zb_trm)    

     call median_filter(xw, yw, zw, wend, x, y, wtr_mdn)  !!wtr_mdn : h + dev_bl
     call median_filter(xb, yb, zb, bend, x, y, bed_mdn)  !!bed_mdn : dev_bl
     
!#######################################################################
     do i = 1, iend
        do j = 1, jend
           if(bed_mdn(i,1).lt.0)then
              bed_mdn(i,1) = bed_mdn(i+1,1)
           elseif(bed_mdn(i,2).lt.0)then
              bed_mdn(i,2) = bed_mdn(i+1,2)
           elseif(bed_mdn(i,j).lt.0)then
              bed_mdn(i,j) = bed_mdn(i,j-1) + (bed_mdn(i,j-1)-bed_mdn(i,j-2))
           endif
        end do
     end do
     
     do i = 1, iend
        do j = 1, jend
           if(wtr_mdn(i,1).lt.0)then
              wtr_mdn(i,1) = wtr_mdn(i+1,1)
           elseif(wtr_mdn(i,2).lt.0)then
              wtr_mdn(i,2) = wtr_mdn(i+1,2)
           elseif(wtr_mdn(i,j).lt.0)then
              wtr_mdn(i,j) = wtr_mdn(i,j-1) + (wtr_mdn(i,j-1)-wtr_mdn(i,j-2))
           endif
        end do
     end do


     do i = 1, iend
        do j = 1, jend
           if(wtr_mdn(i,j).lt.0)then
              wtr_mdn(i,j) = wtr_mdn(i+1,j)
           endif
        end do
     end do

     do i = 1, iend
        do j = 1, jend
           if(bed_mdn(i,j).lt.0)then
              bed_mdn(i,j) = bed_mdn(i+1,j)
           endif
        end do
     end do

! ######################################################################

     call add_grad(iend, jend, x, y, wtr_mdn, ib, wl_mdn)  !!wl_mdn : h + z
     call add_grad(iend, jend, x, y, bed_mdn, ib, bel_mdn) !!bel_mdn : z
     
     call cal_grad_xy(iend, jend, x, y, wl_mdn, dwldx, dwldy)     
     call cal_grad_xy(iend, jend, x, y, bel_mdn, dbldx, dbldy)
     
!!$     call output_csv(wtr_mdn)
     call output_schalar_vtk
     call output_vector_vtk
     
     if(k.eq.1)then
        call drw_cntr
     endif
     
     call delete_arry     
  enddo
     
end program main



! Input parameter
subroutine input_prmtr
  use filter

  open(100,file = '../Target_file.txt')
  read(100,*)case
  close(100)

  ib       = 1./180.

  xtrm_min = 0.2d0  
!!$  xtrm_max =  7.85d0
  xtrm_max = 8.2d0
  
!!$ ---- 20cm ----
!!$ ytrm_min = -0.04d0
!!$ ytrm_max =  0.161d0

!!$ ---- 35cm ----
  ytrm_min = -0.121d0
  ytrm_max =  0.301d0

!!$ ---- 40cm ----
!!$  ytrm_min  = -0.111d0
!!$  ytrm_max  =  0.311d0    

!!$ ---- 45cm ----
!!$  ytrm_min  = -0.141d0
!!$  ytrm_max  =  0.311d0    
  
  intrvl_i = 1                !size of 1 section (cm) 
  intrvl_j = 1  
  
  dx = 0.02d0
  dy = 0.02d0

  iend = int((xtrm_max-xtrm_min)/dx) + 1 
  jend = int((ytrm_max-ytrm_min)/dy) + 1

  call system('mkdir ../00_input_vtk')
  call system('mkdir ../00_input_vtk/'//trim(case))
  call system('mkdir ../00_input_vtk/'//trim(case)//'/vtk_schalar_median')
  call system('mkdir ../00_input_vtk/'//trim(case)//'/vtk_vector_median')    

end subroutine input_prmtr



! Read file name
subroutine read_fln
  use filter

  call system('ls -1 ../mocm_vtk/'//trim(case)//'/ > ./file.txt')

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



! Structuring
subroutine structuring
  use filter

  allocate(x(iend,jend), y(iend,jend))
  allocate(wtr_mdn(iend,jend), bed_mdn(iend,jend))  
  allocate(wl_mdn(iend,jend), bel_mdn(iend,jend))
  allocate(dwldx(iend,jend), dwldy(iend,jend))
  allocate(dbldx(iend,jend), dbldy(iend,jend))
  
  do i = 1, iend
     do j = 1, jend
        x(i,j) = ((i-1) * dx) + xtrm_min
        y(i,j) = ((j-1) * dy) + ytrm_min
     enddo
  enddo  
  
end subroutine structuring



! Set txt file
subroutine set_txt_file
  use filter  

  call system('mv ../mocm_vtk/'//trim(case)//'/'//trim(fln(k))//&
       '/bed_srfc/* ../mocm_vtk/'//trim(case)//'/'//trim(fln(k))//'/bed_srfc/bed.txt')
  call system('mv ../mocm_vtk/'//trim(case)//'/'//trim(fln(k))//&
       '/wtr_srfc/* ../mocm_vtk/'//trim(case)//'/'//trim(fln(k))//'/wtr_srfc/wtr.txt')

  pss_wtr = '../mocm_vtk/'//trim(case)//'/'//trim(fln(k))//'/wtr_srfc/wtr.txt'  
  pss_bed = '../mocm_vtk/'//trim(case)//'/'//trim(fln(k))//'/bed_srfc/bed.txt'

!!$  pss_csv = '../input_vtk/'//trim(case)//'/csv/'//trim(fln(k))//'.csv'
  pss_shr = '../00_input_vtk/'//trim(case)//'/vtk_schalar_median'
  pss_vct = '../00_input_vtk/'//trim(case)//'/vtk_vector_median'  
  
end subroutine set_txt_file



! Read array
subroutine read_arry(nend, pss)
  use filter
  integer                 :: n
  integer,intent(out)     :: nend
  character(*),intent(in) :: pss
  
  open(100,file = trim(pss))  
  read(100,*)nend
  write(*,*)'File name       : ', trim(pss)
!!$  write(*,*)'Number of point : ', nend
  close(100)

end subroutine read_arry



! Set array
subroutine set_arry
  use filter

  allocate(xw(wend), yw(wend), zw(wend))
  allocate(xb(bend), yb(bend), zb(bend))
  
end subroutine set_arry



! Delete array
subroutine delete_arry
  use filter

  deallocate(xw, yw, zw)
  deallocate(xb, yb, zb)
  deallocate(xw_trm, yw_trm, zw_trm)
  deallocate(xb_trm, yb_trm, zb_trm)
  
end subroutine delete_arry



! Read txt
subroutine read_txt(nend, xx, yy, zz, pss)
  integer                 :: n
  integer,intent(in)      :: nend
  real*8,intent(out)      :: xx(nend), yy(nend), zz(nend)
  character(*),intent(in) :: pss

  open(100,file = trim(pss))    
  read(100,*)
  do n = 1, nend
     read(100,*)xx(n), yy(n), zz(n)
  end do
  xx = xx * 0.01d0
  yy = yy * 0.01d0
  zz = zz * 0.01d0
  close(100)

end subroutine read_txt



! Count array
subroutine count_arry(nend, nend_trm, xn, yn)
  use filter
  integer                 :: n, ncount  
  integer,intent(in)      :: nend
  integer,intent(out)     :: nend_trm
  real*8,intent(in)       :: xn(nend), yn(nend)

  ncount = 0
  do n = 1, nend
     if(xn(n).ge.xtrm_min .and. xn(n).le.xtrm_max .and. yn(n).ge.ytrm_min .and. yn(n).le.ytrm_max)then
        ncount = ncount + 1
     endif
  enddo
  nend_trm = ncount

end subroutine count_arry



! Set array 2
subroutine set_arry2
  use filter

  allocate(xw_trm(wend_trm), yw_trm(wend_trm), zw_trm(wend_trm))
  allocate(xb_trm(bend_trm), yb_trm(bend_trm), zb_trm(bend_trm))

end subroutine set_arry2



! Trim txt data
subroutine trim_txt(nend, xn, yn, zn, xn_trm, yn_trm, zn_trm)
  use filter
  integer            :: n, ncount
  integer,intent(in) :: nend
  real*8,intent(in)  :: xn(nend), yn(nend), zn(nend)
  real*8,intent(out) :: xn_trm(nend), yn_trm(nend), zn_trm(nend)  

  ncount = 0
  do n = 1, nend
     if(xn(n).ge.xtrm_min .and. xn(n).le.xtrm_max .and. yn(n).ge.ytrm_min .and. yn(n).le.ytrm_max)then
        ncount = ncount + 1
        xn_trm(ncount) = xn(n)
        yn_trm(ncount) = yn(n)
        zn_trm(ncount) = zn(n)                
     endif
  enddo
  
end subroutine trim_txt



! Median filter
subroutine median_filter(xn, yn, zn, nend, xx, yy, zz)
  use filter
  integer            :: n, ncount
  integer,intent(in) :: nend
  real*8             :: x1, x2, y1, y2, mdn
  real*8,intent(in)  :: xn(nend), yn(nend), zn(nend)  
  real*8,intent(in)  :: xx(iend, jend), yy(iend, jend)
  real*8,intent(out) :: zz(iend, jend)
  real*8,allocatable :: za(:)  
  
  do i = 1, iend
     do j = 1, jend
        x1 = x(i,j) - (dx/2.)
        x2 = x(i,j) + (dx/2.)        
        y1 = y(i,j) - (dy/2.)
        y2 = y(i,j) + (dy/2.)
       
        ncount = 0
        do n = 1, nend
           if(xn(n).gt.x1.and.xn(n).lt.x2.and.yn(n).gt.y1.and.yn(n).lt.y2)then
              ncount = ncount + 1
           endif
        enddo

        allocate(za(0:ncount-1))

        ncount = 0        
        do n = 1, nend
           if(xn(n).gt.x1.and.xn(n).lt.x2.and.yn(n).gt.y1.and.yn(n).lt.y2)then
              za(ncount) = zn(n)
              ncount = ncount + 1              
           endif
        enddo

        zz(i,j) = 0.d0        
        if(ncount.eq.0) then
           zz(i,j) = -0.1d0
        else
           call sort_ascndng_ordr_and_cal_median(ncount, za, mdn)           
           zz(i,j) = mdn
        endif
        
        deallocate(za)      
     enddo
  enddo

end subroutine median_filter



! Sort ascending order
subroutine sort_ascndng_ordr_and_cal_median(iend, val, mdn)
  implicit none
  integer              :: i, j
  integer,intent(in)   :: iend  
  real*8               :: temp, mdn
  real*8,intent(inout) :: val(0:iend-1)

  do i = iend-1, 1, -1 
     do j = 0, i -1
        if(val(j) > val(j+1)) then
           temp = val(j)
           val(j) = val(j+1)
           val(j+1) = temp
        endif
     enddo
  enddo
  mdn = val(int(iend/2))

end subroutine sort_ascndng_ordr_and_cal_median



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



! Calculate of gradient
subroutine cal_grad_xy(iend, jend, x, y, f, dfdx, dfdy)
  implicit none   
  integer            :: i, j
  integer,intent(in) :: iend, jend
  real*8,intent(in)  :: x(1:iend,1:jend), y(1:iend,1:jend), f(1:iend,1:jend)
  real*8,intent(out) :: dfdx(1:iend,1:jend), dfdy(1:iend,1:jend)

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

end subroutine cal_grad_xy




!!$! Output csv data
!!$subroutine output_csv(f)
!!$  use filter
!!$  real*8,intent(in) :: f(iend,jend)
!!$
!!$  open(200, file=trim(pss_csv))    
!!$  write(200,*)iend, ',', jend
!!$  do j = 1, jend     
!!$     do i = 1, iend
!!$        write(200,'(f0.5, a, f0.5, a, f0.5)')x(i,j), ',', y(i,j), ',', f(i,j)
!!$     enddo
!!$  enddo
!!$  close(200)
!!$
!!$end subroutine output_csv



! Output schalar data to vtk
subroutine output_schalar_vtk
  use filter  

  open( 100,file = trim(pss_shr)//'/'//trim(fln(k))//'_trm.vtk')
  write(100,'(a)') '# vtk DataFile Version 3.1'
  write(100,'(a,i4.4,a)') 'out.vtk'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET STRUCTURED_GRID'
  write(100,'(a,5x,i4,5x,i4,5x,i4)')'DIMENSIONS', iend, jend, 1
  write(100,'(a,5x,i8,5x,a)')'POINTS', iend * jend,'float'

  do j = 1, jend
     do i = 1, iend
        write(100,*)x(i,j), y(i,j), 0.d0
     enddo
  enddo

  write(100,*)
  write(100,'(a,5x,i8)')'POINT_DATA', iend * jend
  write(100,'(a,5x,i8)')'FIELD FieldData', 5
 
  call out_scalar(wtr_mdn(:,:), 'deviation_of_water_level')
  call out_scalar(wl_mdn,  'water_level')
  call out_scalar(bed_mdn, 'deviation_of_bed_level')
  call out_scalar(bel_mdn, 'bed_level')
  call out_scalar(wtr_mdn(:,:) - bed_mdn(:,:), 'depth_of_water')  

  
  close(100)
  
end subroutine output_schalar_vtk



! Output scalar data
subroutine out_scalar(f, name)
  use filter  
  real*8,intent(in)       :: f(iend, jend)
  character(*),intent(in) :: name

  write(100,'(a,5x,i4,5x,i8,5x,a)')name, 1, iend*jend, 'double'
!!$  do j = 1, jend
  do j = jend,1,-1
     do i = 1, iend
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
  
  open( 100,file = trim(pss_vct)//'/'//trim(fln(k))//'_trm.vtk')
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
  
!!$  do j = 1, jend
!!$     do i = 1, iend
!!$        write(*,*)wtr_mdn(i,j)
!!$     enddo
!!$  enddo

  write(100,*)
  write(100,'(a,5x,i8)')'POINT_DATA', ijc
  write(100,'(a,5x,i8)')'FIELD FieldData', 1
 
  call out_vector(dwldx, dwldy, 'gradient_of_water_level')
  call out_vector(dbldx, dbldy, 'gradient_of_bed_level')  
  
  close(100)

end subroutine output_vector_vtk



! Output_vector
subroutine out_vector(fx, fy, name)
  use filter
  real*8,intent(in)       :: fx(iend,jend), fy(iend,jend)
  character(*),intent(in) :: name

  write(100,'(a,5x,i4,5x,i8,5x,a)')name, 3, ijc, 'double'
!!$  do j = 1, jend, intrvl_j
  do j = jend, 1, -1*intrvl_j
     do i = 1, iend, intrvl_i
        if( fx(i,j) .lt. 0 )then
           write(100,*)0, fy(i,j), 0
        else
           write(100,*)fx(i,j), fy(i,j), 0
        endif
     end do
  end do

end subroutine out_vector

