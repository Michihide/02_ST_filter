Module filter
  implicit none
  
  integer :: i, j, k
  integer :: ic, jc, ijc, intrvl_i, intrvl_j
  integer :: iend, jend, wend, bend, kend
  integer :: wend_trm, bend_trm, kend_trm
  integer :: smthng, clrnm

  real*8  :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8  :: pxmin, pxmax, pymin, pymax
  real*8  :: dx, dy, ib, br  
  real*8  :: xtrm_min, xtrm_max, ytrm_min, ytrm_max

  real*8,dimension(:),allocatable   :: xw, yw, zw, xb, yb, zb, xk, yk, zk
  real*8,dimension(:),allocatable   :: xw_trm, yw_trm, zw_trm, xb_trm, yb_trm, zb_trm
  real*8,dimension(:),allocatable   :: xk_trm, yk_trm, zk_trm
  real*8,dimension(:,:),allocatable :: x, y
  real*8,dimension(:,:),allocatable :: wtr_mdn, bed_mdn
  real*8,dimension(:,:),allocatable :: wl_mdn, bel_mdn
  real*8,dimension(:,:),allocatable :: dwldx, dwldy
  real*8,dimension(:,:),allocatable :: dbldx, dbldy  
  
  character*100 :: pss_wtr, pss_bed, pss_wtr_trm, pss_bed_trm, pss_csv, pss_shr, pss_vct
  character*50  :: chr, case, stage, smthng_chr, j_number  
  character*50  :: wtr_number, bed_number
  character*50  :: barname, axis_t, axis_l, axis_b
  character*150 :: ifile, ps2pdf, mv, rm , rot
  integer       :: plparseopts_rc
  character*50,allocatable :: fln(:)  

End Module filter
