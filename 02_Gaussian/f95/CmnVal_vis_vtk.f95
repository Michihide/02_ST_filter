Module filter
  implicit none
  
  integer :: i, j, k, l,m
  integer :: ic, jc, ijc, intrvl_i, intrvl_j  
  integer :: iend, jend, kend, mend
  integer :: smthng, clrnm

  real*8  :: ib, br
  real*8  :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8  :: pxmin, pxmax, pymin, pymax, pdy

  real*8,dimension(:),allocatable   :: dis, x_itp, y_itp, bel_itp, bed_itp, wl_itp
  real*8,dimension(:,:),allocatable :: x, y
  real*8,dimension(:,:),allocatable :: wtr_mdn, wl_mdn, dep_mdn
  real*8,dimension(:,:),allocatable :: wtr_gaus, wl_gaus, dev_wl_gaus
  real*8,dimension(:,:),allocatable :: wtr_dpth
  real*8,dimension(:,:),allocatable :: bed_mdn, bel_mdn
  real*8,dimension(:,:),allocatable :: bed_gaus, bel_gaus
  real*8,dimension(:,:),allocatable :: dwldx_mdn, dwldy_mdn
  real*8,dimension(:,:),allocatable :: beldx_mdn, beldy_mdn
  real*8,dimension(:,:),allocatable :: beddx_mdn, beddy_mdn 
  real*8,dimension(:,:),allocatable :: dwldx_gaus, dwldy_gaus
  real*8,dimension(:,:),allocatable :: beldx_gaus, beldy_gaus
  real*8,dimension(:,:),allocatable :: u, v
  real*8,dimension(:,:),allocatable :: tausx, tausy,taus
  real*8,dimension(:,:),allocatable :: an_x, an_y
  
  character*50  :: chr, case, j_number
  character*100 :: barsname, axis_t, axis_l, axis_b    

  character*50,allocatable :: fln(:)  

End Module filter
