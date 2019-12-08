! Draw 2d contour
subroutine drw_cntr
  use filter
  use plplot

  plparseopts_rc = plparseopts(PL_PARSE_FULL)
  if(plparseopts_rc .ne. 0) stop "plparseopt error"
  call plscol0( 0,255,255,255)
  call plscol0(15,  0,  0,  0)
  call plscmap0n(0)
  call plinit     

  barname = 'dev_depth(m)'  
  kend = wend
  allocate(xk(kend), yk(kend), zk(kend))       
  xk = xw
  yk = yw
  zk = zw

  kend_trm = wend_trm
  allocate(xk_trm(kend_trm), yk_trm(kend_trm), zk_trm(kend_trm))
  xk_trm = xw_trm
  yk_trm = yw_trm
  zk_trm = zw_trm


!!$  barname = 'dev_bed(m)'
!!$  kend = bend
!!$  allocate(xk(kend), yk(kend), zk(kend))       
!!$  xk = xb
!!$  yk = yb
!!$  zk = zb    
!!$
!!$  kend_trm = bend_trm
!!$  allocate(xk_trm(kend_trm), yk_trm(kend_trm), zk_trm(kend_trm))
!!$  xk_trm = xb_trm
!!$  yk_trm = yb_trm
!!$  zk_trm = zb_trm    

  xmin = minval(xk) - (maxval(xk)-minval(xk))*0.02d0
  xmax = maxval(xk) + (maxval(xk)-minval(xk))*0.02d0
  ymin = minval(yk) - (maxval(yk)-minval(yk))*0.10d0
  ymax = maxval(yk) + (maxval(yk)-minval(yk))*0.10d0
  zmin = minval(zk) - (maxval(zk)-minval(zk))*0.10d0
  zmax = maxval(zk) + (maxval(zk)-minval(zk))*0.10d0

!!$  zmin  = 0.d0
!!$  zmax  = 5.d0     
  
  pdy   = 0.30d0
  
  write(*,*)minval(yk),maxval(yk)
  
  axis_l = 'Width(m)'
  axis_t = ''
  axis_b = ''

  clrnm = 3

  pxmin = 0.15d0
  pxmax = 0.85d0
  pymin = 0.70d0
  pymax = 0.90d0

  call pladv(0)
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_point_color(kend, xk, yk, zk, 1)
  call plot_rctngl(minval(x), maxval(x), minval(y), maxval(y), 3)
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.05d0, pxmax+0.06d0, pymin, pymax, trim(barname))      

  axis_t   = ''

  pymin = pymin - pdy
  pymax = pymax - pdy

  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_point_color(kend_trm, xk_trm, yk_trm, zk_trm, 1)
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.05d0, pxmax+0.06d0, pymin, pymax, trim(barname))      


  axis_t   = ''

  pymin = pymin - pdy
  pymax = pymax - pdy

  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, dx, dy, wtr_mdn)
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.05d0, pxmax+0.06d0, pymin, pymax, trim(barname))      
  
  call plend

  
!!$  ps2pdf = 'ps2pdf figs/'//trim(case)//'_'//trim(fln(n))//'.ps &
!!$       figs/'//trim(case)//'_'//trim(fln(n))//'.pdf'
!!$  rm = 'rm figs/'//trim(case)//'_'//trim(fln(n))//'.ps'  
!!$  call system(trim(ps2pdf)//' && '//trim(rm))
  
end subroutine drw_cntr



! Set window
subroutine set_window(axmin, axmax, aymin, aymax, axis_b, axis_l, axis_t, pxmin, pxmax, pymin, pymax)
  use plplot
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: axis_b, axis_l, axis_t

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax ,aymin, aymax)
  call plcol0(15)
  call plbox('bcnst', 0.0d0, 5, 'bcnstv', 0.0d0, 0)
  call plmtex('b', 3.0d0, 0.5d0, 0.5d0, trim(axis_b))
  call plmtex('l', 4.5d0, 0.5d0, 0.5d0, trim(axis_l))
  call plmtex('t', 1.d0, 0.5d0, 0.5d0, trim(axis_t))

end subroutine set_window



! Plot rectanguler
subroutine plot_rctngl(xmin, xmax, ymin, ymax, icol)
  use plplot
  integer,intent(in)    :: icol
  real*8,intent(in)     :: xmin, xmax, ymin, ymax
  real*8,dimension(0:4) :: xp, yp

  xp(0) = xmin
  yp(0) = ymin
  xp(1) = xmax
  yp(1) = yp(0)
  xp(2) = xp(1)
  yp(2) = ymax
  xp(3) = xp(0)
  yp(3) = yp(2)
  xp(4) = xp(0)
  yp(4) = yp(0)  
  call plcol0(icol)
  call plline(xp, yp)
  call plpoin(xp, yp, 5)  

end subroutine plot_rctngl



! Plot grid
subroutine plot_grid(iend, jend, xmin, ymin, dx, dy)
  use plplot
  implicit none  
  integer             :: i, j
  integer,intent(in)  :: iend, jend  
  real*8,intent(in)   :: dx, dy, xmin, ymin
  real*8              :: xp(0:iend,0:jend), yp(0:iend,0:jend)    

  do i = 0, iend
     do j = 0, jend
        xp(i,j) = (i * dx) + xmin
        yp(i,j) = (j * dy) + ymin
     enddo
  enddo  

  xp(:,:) = xp(:,:) - (dx/2.d0)
  yp(:,:) = yp(:,:) - (dy/2.d0)
  
  do i = 0, iend
     call plot_line(jend+1, xp(i,0:jend), yp(i,0:jend), 9)
  enddo

  do j = 0, jend
     call plot_line(iend+1, xp(0:iend,j), yp(0:iend,j), 9)
  enddo
  
end subroutine plot_grid



! Set color of contour
subroutine set_clr_cntr(Clr_nm, zmin, zmax)
  use Val_Cmap
  integer,intent(in) :: Clr_nm
  real*8,intent(in)  :: zmin, zmax

  iclrDst = 40
  cqmin = zmin
  cqmax = zmax
  
  if(Clr_nm.eq.3)then
     iru(1) = 0
     igu(1) = 0
     ibu(1) = 255

     iru(2) = 255
     igu(2) = 255
     ibu(2) = 255

     iru(3) = 255
     igu(3) = 0
     ibu(3) = 0
  elseif(Clr_nm.eq.5)then
     iru(1) = 0
     igu(1) = 0
     ibu(1) = 255
     iru(2) = 0
     igu(2) = 255
     ibu(2) = 255
     iru(3) = 0
     igu(3) = 255
     ibu(3) = 20
     iru(4) = 255
     igu(4) = 255
     ibu(4) = 0
     iru(5) = 255
     igu(5) = 0
     ibu(5) = 0
  endif

  if(Clr_nm .eq. 2)call remaprgb_TwoClr
  if(Clr_nm .eq. 3)call remaprgb_ThreClr
  if(Clr_nm .eq. 4)call remaprgb_FourClr
  if(Clr_nm .eq. 5)call remaprgb_FiveClr

end subroutine set_clr_cntr



! Draw contour
subroutine plot_cntr(imin, imax, jmin, jmax, x, y, dx, dy, z)
  use Val_Cmap
  integer               :: i, j
  integer,intent(in)    :: imin, imax, jmin, jmax
  real*8                :: zave
  real*8,intent(in)     :: dx, dy  
  real*8,dimension(0:4) :: xp, yp
  real*8,dimension(imin:imax,jmin:jmax),intent(in) :: x, y, z   

  do j = jmin, jmax
     do i = imin, imax
        xp(0) = x(i,j) - (dx/2.d0)
        yp(0) = y(i,j) - (dy/2.d0)
        xp(1) = x(i,j) + (dx/2.d0)
        yp(1) = yp(0)
        xp(2) = xp(1)
        yp(2) = y(i,j) + (dy/2.d0)
        xp(3) = xp(0)
        yp(3) = yp(2)
        xp(4) = xp(0)
        yp(4) = yp(0)
        call clrdqnty(z(i,j))
        if(z(i,j).ne.0)then        
           call plfill(xp, yp)
        else
        endif
     enddo
  enddo

!!$  do j = jmin, jmax-1
!!$     do i = imin, imax-1
!!$        xp(0) = x(i,j)
!!$        yp(0) = y(i,j)
!!$        xp(1) = x(i+1,j)
!!$        yp(1) = y(i+1,j)
!!$        xp(2) = x(i+1,j+1)
!!$        yp(2) = y(i+1,j+1)
!!$        xp(3) = x(i,j+1)
!!$        yp(3) = y(i,j+1)
!!$        xp(4) = x(i,j)
!!$        yp(4) = y(i,j)
!!$        zave = sum(z(i:i+1,j:j+1)) / 4.d0
!!$        call clrdqnty(zave)        
!!$        call plfill(xp, yp)
!!$     enddo
!!$  enddo

end subroutine plot_cntr



! Draw color bars
subroutine colbar(axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax, barname)
  use Val_Cmap
  integer                 :: i, j 
  real*8                  :: cqt, xp(0:4), yp(0:4)
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: barname

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax, aymin, aymax)
  call plcol0(15)
  call plbox('bc', 0.0d0, 0, 'bcmst', 0.0d0, 0)      
  call plmtex('l', 1.d0, 0.5d0, 0.5d0, barname)   

  dpRng   = aymax - aymin
  RngDlt  = dpRng / iclrDst
  
  do i = 0, iclrDst-1
     xp(0) = axmin
     yp(0) = aymin + RngDlt * (i)
     xp(1) = axmax
     yp(1) = yp(0)
     xp(2) = xp(1)
     yp(2) = aymax + RngDlt * (i)
     xp(3) = xp(0)
     yp(3) = yp(2)
     xp(4) = xp(0)
     yp(4) = yp(0)
     cqt   = rng(i)
     call clrdqnty(cqt)
     call plfill(xp, yp)
  enddo

end subroutine colbar



! Plot point
subroutine plot_point(ni, xp, yp, icol, ptype)
  use plplot
  implicit none
  integer,intent(in) :: icol, ni, ptype
  real*8,intent(in)  :: xp(0:ni-1), yp(0:ni-1)

  call plcol0(icol)
  call plpoin(xp, yp, ptype)

end subroutine plot_point



! Plot line
subroutine plot_line(ni, xp, yp, icol)
  use plplot
  implicit none
  integer,intent(in) :: icol, ni
  real*8,intent(in)  :: xp(0:ni-1), yp(0:ni-1)

  call plcol0(icol)
  call plline(xp, yp)

end subroutine plot_line



! Plot point
subroutine plot_point_color(ni, x, y, z, ptype)
  use plplot
  implicit none
  integer            :: i
  integer,intent(in) :: ni, ptype
  real*8,intent(in)  :: x(0:ni-1), y(0:ni-1), z(0:ni-1)
  real*8             :: xp(0:1), yp(0:1)

  do i = 0, ni-1
     call clrdqnty(z(i))
     xp(0) = x(i)
     yp(0) = y(i)     
     call plpoin(xp, yp, ptype)
  enddo
     
end subroutine plot_point_color
