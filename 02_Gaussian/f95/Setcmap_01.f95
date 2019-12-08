!===================================================================================================================================
subroutine SetClrRng
!===================================================================================================================================
  use Val_Cmap
!-----------------------------------------------------------------------------------------------------------------------------------
  dpRng   = cqmax - cqmin   ! Colored Data Range
  RngDlt  = dpRng / iclrDst
  k = 0
  rng(k) = cqmin
  do k = 1, iclrDst
     rng(k) = rng(k-1) + RngDlt
  enddo
  ! write(*,*)  
  ! write(*,*) 'cqmax, cqmin, dpRng, iclrDst:'
  ! write(*,*)  cqmax, cqmin, dpRng, iclrDst
  ! write(*,*)
!-----------------------------------------------------------------------------------------------------------------------------------
end subroutine SetClrRng
!===================================================================================================================================
subroutine clrdqnty(cqt)
!===================================================================================================================================
  use Val_Cmap
  real*8 :: cqt
!-----------------------------------------------------------------------------------------------------------------------------------
  ic = idf+1
  do k = 1, iclrDst-1
     if(cqt.ge.rng(k) .and. cqt.lt.rng(k+1))then
        ic = k + idf !+ 1
        goto 220
     endif
  enddo
220 if(cqt.gt.rng(iclrDst)) ic = idf + iclrDst
  icc = idf + iclrDst
  if(ic.ge.icc) ic = icc - 1
  call plcol0(ic)
!-----------------------------------------------------------------------------------------------------------------------------------
 end subroutine clrdqnty
!===================================================================================================================================
 subroutine remaprgb_TwoClr
!===================================================================================================================================
   use Val_Cmap
!-----------------------------------------------------------------------------------------------------------------------------------
   open(90, file='SetGradtnColor.dat')
   call SetClrRng
   is  = 0
   idf = 16
   rc = 1. / iclrDst
   rr = abs(iru(1) - iru(2)) * rc
   gg = abs(igu(1) - igu(2)) * rc
   bb = abs(ibu(1) - ibu(2)) * rc
   if(iru(1) .gt. iru(2)) rr = -rr
   if(iru(1) .lt. iru(2)) rr =  rr
   if(igu(1) .gt. igu(2)) gg = -gg
   if(igu(1) .lt. igu(2)) gg =  gg
   if(ibu(1) .gt. ibu(2)) bb = -bb
   if(ibu(1) .lt. ibu(2)) bb =  bb
   write(90, *) is, ie, ig, iclrDst, rc
   write(90, *) iru(1), igu(1), ibu(1)
   write(90, *) iru(2), igu(2), ibu(2)
   write(90, *) rr, gg, bb
   write(90, *)
   do i = is, iclrDst
      r(i+idf) = iru(1) + rr * i
      g(i+idf) = igu(1) + gg * i
      b(i+idf) = ibu(1) + bb * i
      write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
   enddo
   do i = 1, idf
      call plgcol0(i-1, r(i), g(i), b(i))
   enddo
   nc = iclrDst + idf
   call plscmap0(r, g, b)
   close(90)
!-----------------------------------------------------------------------------------------------------------------------------------
end subroutine remaprgb_TwoClr
!===================================================================================================================================
subroutine remaprgb_ThreClr
!===================================================================================================================================
  use Val_Cmap
!-----------------------------------------------------------------------------------------------------------------------------------
  open(90, file='SetGradtnColor.dat')
  call SetClrRng
  is  = 0
  idf = 16
  ig = iclrDst * 0.5
  ig = int(iclrDst / 2)
  rc = 1. / iclrDst
  icg = ig * 2
  if(iclrDst.ne.icg) then
     write(*,*) 'Warning!'
     write(*,*) 'iclrDst must set to (mod(iclrDst,2).eq.0.'
     write(*,*) 'Warning!'
     write(*,*) 
  endif
  rr = abs(iru(1) - iru(2)) * rc * 2.
  gg = abs(igu(1) - igu(2)) * rc * 2.
  bb = abs(ibu(1) - ibu(2)) * rc * 2.
  if(iru(1) .gt. iru(2)) rr = -rr
  if(iru(1) .lt. iru(2)) rr =  rr
  if(igu(1) .gt. igu(2)) gg = -gg
  if(igu(1) .lt. igu(2)) gg =  gg
  if(ibu(1) .gt. ibu(2)) bb = -bb
  if(ibu(1) .lt. ibu(2)) bb =  bb
  write(90, *) is, ie, ig, iclrDst, rc
  write(90, *) iru(1), igu(1), ibu(1)
  write(90, *) iru(2), igu(2), ibu(2)
  write(90, *) iru(3), igu(3), ibu(3)
  write(90, *) rr, gg, bb
  write(90, *)
  do i = is, ig
     r(i+idf) = iru(1) + rr * i
     g(i+idf) = igu(1) + gg * i
     b(i+idf) = ibu(1) + bb * i
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  rr = abs(iru(2) - iru(3)) * rc * 2.
  gg = abs(igu(2) - igu(3)) * rc * 2.
  bb = abs(ibu(2) - ibu(3)) * rc * 2.
  if(iru(2) .gt. iru(3)) rr = -rr
  if(iru(2) .lt. iru(3)) rr =  rr
  if(igu(2) .gt. igu(3)) gg = -gg
  if(igu(2) .lt. igu(3)) gg =  gg
  if(ibu(2) .gt. ibu(3)) bb = -bb
  if(ibu(2) .lt. ibu(3)) bb =  bb
  write(90, *) rr, gg, bb
  write(90, *)
  do i = ig+1, iclrDst
     ii = i - ig
     r(i+idf) = iru(2) + rr * ii
     g(i+idf) = igu(2) + gg * ii
     b(i+idf) = ibu(2) + bb * ii
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  do i = 1, idf
     call plgcol0(i-1, r(i), g(i), b(i))
  enddo
  nc = iclrDst + idf
  call plscmap0(r, g, b)
  close(90)
!-----------------------------------------------------------------------------------------------------------------------------------
end subroutine remaprgb_ThreClr
!===================================================================================================================================
subroutine remaprgb_FourClr
!===================================================================================================================================
  use Val_Cmap
!-----------------------------------------------------------------------------------------------------------------------------------
  open(90, file='SetGradtnColor.dat')
  call SetClrRng
  is  = 0
  idf = 16
  ig = 33
  ig = int(iclrDst / 3)
  rc = 1. / iclrDst
  icg = ig * 3
  if(iclrDst.ne.icg) then
     write(*,*) 'Warning!'
     write(*,*) 'iclrDst must set to (mod(iclrDst,3).eq.0.'
     write(*,*) 'Warning!'
     write(*,*) 
  endif
  r = 0
  rr = abs(iru(1) - iru(2)) * rc * 3.
  gg = abs(igu(1) - igu(2)) * rc * 3.
  bb = abs(ibu(1) - ibu(2)) * rc * 3.
  if(iru(1) .gt. iru(2)) rr = -rr
  if(iru(1) .lt. iru(2)) rr =  rr
  if(igu(1) .gt. igu(2)) gg = -gg
  if(igu(1) .lt. igu(2)) gg =  gg
  if(ibu(1) .gt. ibu(2)) bb = -bb
  if(ibu(1) .lt. ibu(2)) bb =  bb
  write(90, *) is, ie, ig, iclrDst, rc
  write(90, *) iru(1), igu(1), ibu(1)
  write(90, *) iru(2), igu(2), ibu(2)
  write(90, *) iru(3), igu(3), ibu(3)
  write(90, *) rr, gg, bb
  write(90, *)
  do i = is, ig
     r(i+idf) = iru(1) + rr * i
     g(i+idf) = igu(1) + gg * i
     b(i+idf) = ibu(1) + bb * i
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  rr = abs(iru(2) - iru(3)) * rc * 3.
  gg = abs(igu(2) - igu(3)) * rc * 3.
  bb = abs(ibu(2) - ibu(3)) * rc * 3.
  if(iru(2) .gt. iru(3)) rr = -rr
  if(iru(2) .lt. iru(3)) rr =  rr
  if(igu(2) .gt. igu(3)) gg = -gg
  if(igu(2) .lt. igu(3)) gg =  gg
  if(ibu(2) .gt. ibu(3)) bb = -bb
  if(ibu(2) .lt. ibu(3)) bb =  bb
  write(90, *) rr, gg, bb
  write(90, *)
  do i = ig+1, ig*2
     ii = i - ig
     r(i+idf) = iru(2) + rr * ii
     g(i+idf) = igu(2) + gg * ii
     b(i+idf) = ibu(2) + bb * ii
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  rr = abs(iru(3) - iru(4)) * rc * 3.
  gg = abs(igu(3) - igu(4)) * rc * 3.
  bb = abs(ibu(3) - ibu(4)) * rc * 3.
  if(iru(3) .gt. iru(4)) rr = -rr
  if(iru(3) .lt. iru(4)) rr =  rr
  if(igu(3) .gt. igu(4)) gg = -gg
  if(igu(3) .lt. igu(4)) gg =  gg
  if(ibu(3) .gt. ibu(4)) bb = -bb
  if(ibu(3) .lt. ibu(4)) bb =  bb
  write(90, *) rr, gg, bb
  write(90, *)
  do i = ig*2+1, iclrDst
     iii = i - ig*2
     r(i+idf) = iru(3) + rr * iii
     g(i+idf) = igu(3) + gg * iii
     b(i+idf) = ibu(3) + bb * iii
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  do i = 1, idf
     call plgcol0(i-1, r(i), g(i), b(i))
  enddo
  nc = iclrDst + idf
  call plscmap0(r, g, b)
  close(90)
!-----------------------------------------------------------------------------------------------------------------------------------
end subroutine remaprgb_FourClr
!===================================================================================================================================
subroutine remaprgb_FiveClr
!===================================================================================================================================
  use Val_Cmap
!-----------------------------------------------------------------------------------------------------------------------------------
  open(90, file='SetGradtnColor.dat')
  call SetClrRng
  is  = 0
  idf = 16
  ig = 50
  ig = int(iclrDst / 4)
  rc = 1. / iclrDst
  icg = ig * 4
  if(iclrDst.ne.icg) then
     write(*,*) 'Warning!'
     write(*,*) 'iclrDst must set to (mod(iclrDst,4).eq.0.'
     write(*,*) 'Warning!'
     write(*,*) 
  endif
  r = 0
  rr = abs(iru(1) - iru(2)) * rc * 4.
  gg = abs(igu(1) - igu(2)) * rc * 4.
  bb = abs(ibu(1) - ibu(2)) * rc * 4.
  if(iru(1) .gt. iru(2)) rr = -rr
  if(iru(1) .lt. iru(2)) rr =  rr
  if(igu(1) .gt. igu(2)) gg = -gg
  if(igu(1) .lt. igu(2)) gg =  gg
  if(ibu(1) .gt. ibu(2)) bb = -bb
  if(ibu(1) .lt. ibu(2)) bb =  bb
  write(90, *) is, ie, ig, iclrDst, rc
  write(90, *) iru(1), igu(1), ibu(1)
  write(90, *) iru(2), igu(2), ibu(2)
  write(90, *) iru(3), igu(3), ibu(3)
  write(90, *) rr, gg, bb
  write(90, *)
  do i = is, ig
     r(i+idf) = iru(1) + rr * i
     g(i+idf) = igu(1) + gg * i
     b(i+idf) = ibu(1) + bb * i
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  rr = abs(iru(2) - iru(3)) * rc * 4.
  gg = abs(igu(2) - igu(3)) * rc * 4.
  bb = abs(ibu(2) - ibu(3)) * rc * 4.
  if(iru(2) .gt. iru(3)) rr = -rr
  if(iru(2) .lt. iru(3)) rr =  rr
  if(igu(2) .gt. igu(3)) gg = -gg
  if(igu(2) .lt. igu(3)) gg =  gg
  if(ibu(2) .gt. ibu(3)) bb = -bb
  if(ibu(2) .lt. ibu(3)) bb =  bb
  write(90, *) rr, gg, bb
  write(90, *)
  do i = ig+1, ig*2
     ii = i - ig
     r(i+idf) = iru(2) + rr * ii
     g(i+idf) = igu(2) + gg * ii
     b(i+idf) = ibu(2) + bb * ii
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  rr = abs(iru(3) - iru(4)) * rc * 4.
  gg = abs(igu(3) - igu(4)) * rc * 4.
  bb = abs(ibu(3) - ibu(4)) * rc * 4.
  if(iru(3) .gt. iru(4)) rr = -rr
  if(iru(3) .lt. iru(4)) rr =  rr
  if(igu(3) .gt. igu(4)) gg = -gg
  if(igu(3) .lt. igu(4)) gg =  gg
  if(ibu(3) .gt. ibu(4)) bb = -bb
  if(ibu(3) .lt. ibu(4)) bb =  bb
  write(90, *) rr, gg, bb
  write(90, *)
  do i = ig*2+1, ig*3
     iii = i - ig*2
     r(i+idf) = iru(3) + rr * iii
     g(i+idf) = igu(3) + gg * iii
     b(i+idf) = ibu(3) + bb * iii
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  rr = abs(iru(4) - iru(5)) * rc * 4.
  gg = abs(igu(4) - igu(5)) * rc * 4.
  bb = abs(ibu(4) - ibu(5)) * rc * 4.
  if(iru(4) .gt. iru(5)) rr = -rr
  if(iru(4) .lt. iru(5)) rr =  rr
  if(igu(4) .gt. igu(5)) gg = -gg
  if(igu(4) .lt. igu(5)) gg =  gg
  if(ibu(4) .gt. ibu(5)) bb = -bb
  if(ibu(4) .lt. ibu(5)) bb =  bb
  write(90, *) rr, gg, bb
  write(90, *)
  do i = ig*3+1, iclrDst
     iii = i - ig*3
     r(i+idf) = iru(4) + rr * iii
     g(i+idf) = igu(4) + gg * iii
     b(i+idf) = ibu(4) + bb * iii
     write(90, *) i, r(i+idf), g(i+idf), b(i+idf), rng(i)
  enddo
  do i = 1, idf
     call plgcol0(i-1, r(i), g(i), b(i))
  enddo
  nc = iclrDst + idf
  call plscmap0(r, g, b)
  close(90)
!-----------------------------------------------------------------------------------------------------------------------------------
end subroutine remaprgb_FiveClr
! !===================================================================================================================================
! subroutine cmap1_init(gray) 
! !===================================================================================================================================
! ! For gray.eq.1, basic grayscale variation from half-dark 
! ! to light.  Otherwise, hue variations around the front of the 
! ! colour wheel from blue to green to red with constant lightness 
! ! and saturation. 
!   implicit none 
!   integer gray 
!   real*8 i(0:1), h(0:1), l(0:1), s(0:1) 
!   integer rev(0:1) 
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! left boundary 
!   i(0) = 0.d0 
! ! right boundary 
!   i(1) = 1.d0 
!   if (gray.eq.1) then 
! ! hue -- low: red (arbitrary if s=0) 
!      h(0) = 0.0d0 
! ! hue -- high: red (arbitrary if s=0) 
!      h(1) = 0.0d0 
! ! lightness -- low: half-dark 
!      l(0) = 0.5d0 
! ! lightness -- high: light 
!      l(1) = 1.0d0 
! ! minimum saturation 
!      s(0) = 0.0d0 
! ! minimum saturation 
!      s(1) = 0.0d0 
!   else 
! ! This combination of hues ranges from blue to cyan to green to yellow to red (front of colour wheel) with constant lightness = 0.6 
! ! and saturation = 0.8. 
! ! hue -- low: blue 
!      h(0) = 240.d0 
! ! hue -- high: red 
!      h(1) = 0.0d0 
! ! lightness -- low: 
!      l(0) = 0.6d0 
! ! lightness -- high:
!      l(1) = 0.6d0 
! ! saturation 
!      s(0) = 0.8d0 
! ! minimum saturation 
!      s(1) = 0.8d0 
!   endif
!   rev(0) = 0 
!   rev(1) = 0 
!   call plscmap1n(256) 
!   call plscmap1l(0, 2, i, h, l, s, rev) 
! !-----------------------------------------------------------------------------------------------------------------------------------
!     end subroutine cmap1_init
! !===================================================================================================================================
