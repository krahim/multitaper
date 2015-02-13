C$$$    The multitaper R package
C$$$    Multitaper and spectral analysis package for R
C$$$    Copyright (C) 2011 Wesley S. Burr, Karim J. Rahim, David J. Thomson 
C$$$
C$$$    This file is part of the multitaper package for R.
C$$$
C$$$    The multitaper package is free software: you can redistribute it and
C$$$    or modify
C$$$    it under the terms of the GNU General Public License as published by
C$$$    the Free Software Foundation, either version 2 of the License, or
C$$$    any later version.
C$$$
C$$$    The multitaper package is distributed in the hope that it will be 
C$$$    useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
C$$$    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C$$$    GNU General Public License for more details.
C$$$
C$$$    You should have received a copy of the GNU General Public License
C$$$    along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
C$$$
C$$$    If you wish to report bugs please contact the author. 
C$$$    wesley.burr@gmail.com

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc This files contains modified spectral estimation code adapted from 
cc Robert Parker's psd.f 

c*********************************************************************
cc 
cc  quickSineF
cc
cc  Simple non-adaptive (possibly weighted) sine taper multitaper
cc  computation program. Explicitly for calling from within 
cc  the adaptive loop of spec.mtm.sine. The R-specific version 
cc  of this (quickSine) runs quickly on its own; it is the adaptive
cc  loops that need speeding up.
cc
cc  Adapted from Robert Parker's 'psd.f'.
cc
c*********************************************************************
      subroutine quickSineF(nFreqs,nFFT,k,cft,useAdapt,kadapt,spec)
      
      implicit none
      
      integer nFreqs,nFFT,k, ks, i, j, i2, j1, j2
      logical useAdapt
      complex*16 cft(nFFT), zz
      real*8 spec(1:nFreqs), ck, wt, kadapt(1:nFreqs)
      
      do 5 j=1,nFreqs
         spec = 0.0d+00
 5    continue   
      
      do 6 i=1,nFreqs
         i2 = 2*(i-1)
         if(useAdapt) then
            ks = int(kadapt(i))
         else
            ks = k
         endif
         
         ck = 1/(real(ks)**2)
         
         do 7 j=1,ks
            j1 = mod(i2+nFFT-j,nFFT)
            j2 = mod(i2+j,nFFT)
            zz = cft(j1+1) - cft(j2+1)
            wt = 1.0d+00 - ck*(j-1)**2
            spec(i) = spec(i) + (dble(zz)**2 + dimag(zz)**2)*wt
 7       continue
         
         spec(i) = spec(i)*(6.0d+00 *dble(ks))/(4*(dble(ks)**2) + 
     *        (3*dble(ks)) -1)
 6    continue
      
      return
      end subroutine
      
c*********************************************************************
cc 
cc  curbF
cc
cc  Rewrites input vector so all points lie below a piecewise
cc  linear function v(k) + abs(j-k); clips strong peaks, keeps
cc  slopes below 1 in magnitude. Based on Robert Parker's 'psd.f'.
cc
c*********************************************************************
 
      subroutine curbF(n, v)
      implicit none
      integer n, j, k
      real*8 v(n), vloc
      
      do 1500 j=2, n-1
         if (v(j).lt.v(j+1) .and. v(j).lt.v(j-1)) then
            vloc=v(j)
            do 1200 k=1, n
               v(k)=min(v(k), vloc+abs(j-k))
 1200       continue
         endif
 1500 continue
      
      return
      end
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc northF
cc
cc Performs quadratically-weighted LS fit to some function 's' by
cc a degree-two polynomial in an orthogonal basis; returns
cc d1 and d2, estimates of 1st and 2nd derivatives at center of record
cc
cc  Taken directly from Robert Parker's 'psd.f'.
cc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine northF(n, i1, i2, s, ds, dds)
      implicit none
      integer i1, i2, n, el, L, kk, i, u0sq
      real*8 gamma, s(n), ds, dds, amid, u1sq, u2sq, dot0, dot1, dot2
     *     ,   ssq
      
      L = i2 - i1 + 1
      el=L
      gamma = (el**2 - 1.0)/12.0
      u0sq = el
      u1sq = el*(el**2 - 1.0)/12.0
      u2sq = (el*(el**2 - 1.0)*(el**2- 4.0))/180.0
      amid= 0.5*(el + 1.0)
      dot0=0.0
      dot1=0.0
      dot2=0.0
      ssq=0.0
      do 1100 kk=1, L
         i=kk + i1 - 1
c     Negative or excessive index uses even function assumption
         if (i.le. 0) i=2 - i
         if (i.gt. n) i=2*n - i
         dot0 = dot0 + s(i)
         dot1 = dot1 + (kk - amid) * s(i)
         dot2 = dot2 + ((kk - amid)**2 - gamma)*s(i)
 1100 continue
      ds = dot1/u1sq
      dds = 2.0*dot2/u2sq
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc adapt
cc
cc Performs adaptive spectral estimation via sine taper approach.
cc From pilot estimate of spectrum, computes estimates of S'' to be used
cc in Eq. (13) of Riedel & Sidorenko (1995).
cc
cc  Adapted somewhat from Robert Parker's 'psd.f'.
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine adapt(ntimes, k, nFreqs, sx, nFFT, cft, df,kopt,fact)
      implicit none
      integer k, ntimes, nFreqs, ispan, iter, nFFT, j
      real*8 sx(nFreqs), kopt(nFreqs), y(nFreqs), dy, ddy, R, ak, phi
     *     , sigR, opt(nFreqs), fact, df, c1, c2
      complex*16 cft(nFFT)
      data c1/1.2000/, c2/3.437/
      
      do 5 j=1,nFreqs
         kopt(j) = k
 5    continue
      
c  Adaptive iteration for MSE spectrum
      do 1600 iter=1, ntimes
c     
         do 1100 j=1, nFreqs
            y(j)= log(sx(j))
 1100    continue
         
c  Estimate K, number of tapers at each freq for MSE spectrum
c  R = S"/S -- use R = Y" + (Y')**2 , Y=ln S.
         do 1200 j=1, nFreqs
c     
            ispan = int(kopt(j)*1.4)
            call northF(nFreqs, j-ispan, j+ispan, y, dy, ddy)
            R = (ddy  + dy**2)/df**2
            ak=kopt(j)/dble(2*ispan)
            phi=720.0*ak**5*(1.0 - 1.286*ak + 
     $           ak**3 - 0.0909*ak**5)
            
            sigR= sqrt(phi/dble(kopt(j))**5) / df**2
            opt(j)=c2/(df**4 *( R**2 + 1.4*sigR**2) /fact**2)** 0.2
 1200    continue
         
         call curbF(nFreqs, opt)
         do 1550 j=1, nFreqs
            kopt(j)=max(3.0, opt(j))
 1550    continue
c     Recompute spectrum with optimal variable taper numbers
         call quickSineF(nFreqs,nFFT,1,cft,.true.,kopt,sx)
 1600 continue
      return
      end
      
