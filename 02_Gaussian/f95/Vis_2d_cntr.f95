! Draw 2d contour
subroutine drw_cntr
  use filter
  use plplot

!!$  if(k .eq. 1)then
  if(k .eq. 1 .and. j .eq. 1)then
     plparseopts_rc = plparseopts(PL_PARSE_FULL)
     if(plparseopts_rc .ne. 0) stop "plparseopt error"
     call plscol0( 0,255,255,255)
     call plscol0(15,  0,  0,  0)
     call plscmap0n(0)
     call plinit
  endif
     call pladv(0)
     call plschr(2.0d0, 2.0d0)
     call plwidth(1.d0)     

  pdy      = 0.23d0     


  ! ---- cntr_of_Deviation_of_bed_level_median ----
  axis_t   = trim(fln(k))//' Median'
  axis_b   = ''
  axis_l   = 'Width [m]'     
  clrnm    = 3
  barsname = 'B.L [m]'

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ymin = minval(y) - (maxval(y)-minval(y))*0.10d0
  ymax = maxval(y) + (maxval(y)-minval(y))*0.10d0
!!$  zmin = minval(bed_mdn) - abs(maxval(bed_mdn) - minval(bed_mdn))*0.1d0
!!$  zmax = maxval(bed_mdn) + abs(maxval(bed_mdn) - minval(bed_mdn))*0.1d0
  zmin = 0.03d0 ; zmax = 0.06d0

  pxmin = 0.10d0    ;    pxmax = 0.90d0
  pymin = 0.78d0    ;    pymax = 0.95d0

  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, bed_mdn)
  call plcol0(9)
  call pljoin(x(1,j), y(1,j), x(iend,j), y(iend,j))
!!$  call plot_rctngl(x(1,1), x(iend,1), y(1,1), y(1,jend), 3)
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.03d0, pxmax+0.04d0, pymin, pymax, trim(barsname))    


  ! ---- cntr_of_Deviation_of_bed_level_gaussian ----
  axis_t   = trim(fln(k))//' Gaussian'

  pymin = pymin - pdy   ;   pymax = pymax - pdy
  
!!$  zmin = minval(bed_gaus) - abs(maxval(bed_gaus) - minval(bed_gaus))*0.1d0
!!$  zmax = maxval(bed_gaus) + abs(maxval(bed_gaus) - minval(bed_gaus))*0.1d0
  zmin = 0.03d0 ; zmax = 0.06d0
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, bed_gaus)
  call plcol0(9)
  call pljoin(x(1,j), y(1,j), x(iend,j), y(iend,j))
  write(*,*)y(iend,j)
!!$  call plot_rctngl(x(1,1), x(iend-120,1), y(1,1), y(1,jend), 3)
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.03d0, pxmax+0.04d0, pymin, pymax, trim(barsname))   



  ! ---- lng_of_water_level_gaussian ----
  axis_l   = 'D.B.L [m]'
  axis_t   = 'j number : '//trim(j_number)
  axis_b   = ''
  pymin = pymin - pdy ; pymax = pymax - pdy
  
!!$  ymin = minval(bed_gaus) - abs(maxval(bed_gaus) - minval(bed_gaus))*0.1d0
!!$  ymax = maxval(bed_gaus) + abs(maxval(bed_gaus) - minval(bed_gaus))*0.1d0
  ymin = 0.03d0 ; ymax = 0.08d0
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call plot_line(iend, x(:,j), wtr_gaus(:,j), 9)
  call plot_line(iend, x(:,j), bed_gaus(:,j), 9)
  call plot_line(iend, x(:,j), wtr_mdn (:,j), 7)
  call plot_line(iend, x(:,j), bed_mdn (:,j), 7)
  
!!$  call plot_line(iend, x(:,j), wtr_mdn (:,jend/2), 7)
!!$  call plot_line(iend, x(:,j), bed_mdn (:,jend/2), 7)
!!$  call plot_line(iend, x(:,j), wtr_gaus(:,jend/2), 9)
!!$  call plot_line(iend, x(:,j), bed_gaus(:,jend/2), 9)
!!$  call plot_line(iend, x(:,j), wtr_mdn (:,jend), 3)
!!$  call plot_line(iend, x(:,j), bed_mdn (:,jend), 3)
!!$  call plot_line(iend, x(:,j), wtr_gaus(:,j), 9)
!!$  call plot_line(iend, x(:,j), bed_gaus(:,j), 9)
!!$
!!$  call plot_line(iend, x(:,j), wtr_mdn (:,j), 3)
!!$  call plot_line(iend, x(:,j), bed_mdn (:,j), 3)
!!$  call plot_line(iend, x(:,j), wtr_gaus(:,j), 9)
!!$  call plot_line(iend, x(:,j), bed_gaus(:,j), 9)  



  ! ---- lng_water_surface_slope ----
  axis_l   = 'dZdx'
  axis_t   = 'Grey:Median, Red:Guassian'     
  axis_b   = 'Distance from upstream(m)'

  pymin = pymin - pdy ; pymax = pymax - pdy
 
  ymin = minval(beddx_mdn(:,1)) - abs(maxval(beddx_mdn(:,1)) - minval(beddx_mdn(:,j)))*0.1d0
  ymax = maxval(beddx_mdn(:,1)) + abs(maxval(beddx_mdn(:,j)) - minval(beddx_mdn(:,j)))*0.1d0
  ymin = 0.0d0 ; ymax = 0.05
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call plot_line(iend, x(:,j), dwldx_mdn (:,j), 7)  
  call plot_line(iend, x(:,j), abs(beddx_gaus(:,j)), 9)
  call plot_line(iend, x(:,j), abs(beddx_mdn(:,j)), 7)
!!$  call plot_line(iend, x(:,j), dwldx_gaus(:,j), 9)

  
!!$  if(k .eq. nend)then
  if(k .eq. nend .and. j.eq.jend)then  
     call plend
  endif

end subroutine drw_cntr




! Set window
subroutine set_window(axmin, axmax, aymin, aymax, axis_b, axis_l, axis_t, pxmin, pxmax, pymin, pymax)
  use plplot
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: axis_b, axis_l, axis_t

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax ,aymin, aymax)
  call plcol0(15)
  call plbox ('bcnst', 0.0d0, 5    , 'bcnstv', 0.0d0, 0)
  call plmtex('b'    , 3.0d0, 0.5d0, 0.5d0   , trim(axis_b))
  call plmtex('l'    , 4.5d0, 0.5d0, 0.5d0   , trim(axis_l))
  call plmtex('t'    , 1.d0 , 0.5d0, 0.5d0   , trim(axis_t))

end subroutine set_window



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
subroutine plot_cntr(imin, imax, jmin, jmax, x, y, z)
  use Val_Cmap
  integer               :: i, j
  integer,intent(in)    :: imin, imax, jmin, jmax
  real*8                :: zave
  real*8,dimension(0:4) :: xp, yp
  real*8,dimension(imin:imax,jmin:jmax),intent(in) :: x, y, z   

  do j = jmin, jmax-1
     do i = imin, imax-1
        xp(0) = x(i,j)
        yp(0) = y(i,j)
        xp(1) = x(i+1,j)
        yp(1) = y(i+1,j)
        xp(2) = x(i+1,j+1)
        yp(2) = y(i+1,j+1)
        xp(3) = x(i,j+1)
        yp(3) = y(i,j+1)
        xp(4) = x(i,j)
        yp(4) = y(i,j)
        zave = sum(z(i:i+1,j:j+1)) / 4.d0
        call clrdqnty(zave)        
        call plfill(xp, yp)
     enddo
  enddo

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
  call plwidth(2.d0)
  call plline(xp, yp)
  call plwidth(1.d0)  

end subroutine plot_line




! Plot line
subroutine plot_line_xave(ni, nj, x, y, icol)
  use plplot
  implicit none
  integer            :: i, j 
  integer,intent(in) :: icol, ni, nj
  real*8,intent(in)  :: x(0:ni-1,0:nj-1), y(0:ni-1,0:nj-1)
  real*8             :: xp(0:nj-1), yp(0:nj-1)  

  do j = 0, nj-1
     xp(j) = x(1,j)     
     yp(j) = sum(y(:,j))/ni
  enddo

  yp(0:nj-1) = yp(nj-1:0:-1)
  ! yp(:) = yp(:) + abs(yp(0))
  
  call plcol0(icol)
  call plwidth(2.d0)
  call plline(xp, yp)
  call plwidth(1.d0)  

end subroutine plot_line_xave



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
