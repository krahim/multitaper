!  slp: Slepian Regression Smoothers
!   
!  Based on LISP code accompanying Percival and Walden 1993, the work of David Thomson and David Slepian.
!  Code: Karim Rahim 2010
!  Modifications: Wesley Burr, small modifications,  2013
!  Maintainer: Karim Rahim 2017
!
!  This file is part of the slp package for R.
!
!  If you wish to report bugs please contact the maintainer:
!  
!
!  dpss.f90 calculates Discrete Prolate Spheroidal Sequences using 
!  LAPACK routines 'dstebz' and 'dstein', using the tridiagonal method.
!
!  Note the subroutine expects the memory for the matrix 'v' and the vector 
!  'ev' to be allocated by the calling program. This is done by the 
!  '.dpss' subroutine in /R.
!
! name change to avoid namespace issue (2017)

subroutine fdpss (n, k, nw, v, ev)
  ! Calculate dpss using the tridiagonal formulation given in 
  ! Percival and Walden (1993), Chapter 8.3, using LAPACK functions
  ! in place of EISPACK.
  !
  ! Also reduces the problem using the trick of reducing the symmetric 
  ! tridiagonal matrix to two half sized matrices: one for the even 
  ! eigenfunctions, and one for the odd. This trick was mentioned in a 
  ! Bell Labs technical memo by Slepian in 1977.
  !

  implicit none
  integer :: n, k,  oddK, evenK, nOdd, nEven, i, j, &
       oddK_2, evenK_2M1, iTest
  double precision :: nw, w, ctpw, pi, twopi, sr2, &
       dlamch, abstol, sqrtsumsq
  parameter(pi=3.141592653589793d0,twopi=2.0d0*pi)
  double precision, pointer :: d(:), e(:), work(:), &
       blockDbleMem(:), evLocal(:)!, vlocal(:,:)
  double precision :: ev(k), v(n,k)

  integer, pointer :: blockIntMem(:)
  logical :: is_evenN
  character :: cmach

  !
  ! interface block added to conform with f90 standard.
  ! In response to bug reported by Brian Ripley, July 25, 2010
  ! - krahim
  ! 
  interface
     subroutine tridiagMatrixEigen(n, k, d, e, v, ldv, ev, &
          abstol, blockIntMem, work)
       implicit none
       character :: range, order, cmach
       integer :: nsplit, m, info
       integer :: k, ldv, n, il
       integer, target :: blockIntMem(5*n+k)
       integer,  pointer :: iblock(:), isplit(:), iwork(:), ifail(:)
       double precision :: vl, vu, abstol
       double precision :: d(n), e(n-1), v(ldv,k),  ev(n), work(5*n)
     end subroutine tridiagMatrixEigen
  end interface
  
  w = nw/dble(n)
  ctpw =  dcos(twopi*w)
  oddK = k/2
  evenK = k - oddK
  sr2 = dsqrt(2.0d0)

  nOdd = n/2
  nEven = n - nOdd
  
  ! Allocate memory and use pointers to reduce malloc calls
  ! The values '5' and '8' provide the required memory space for the two 
  ! LAPACK calls. See the documentation for 'dstebz' and 'dstein'.
  ! Memory used in this procedure is allocated in the following two 
  ! calls and then appropriate blocks of memory are accessed using 
  ! pointers. Note: memory space for the matrix used by the calling 
  ! program must be allocated in the calling program

  allocate(blockIntMem(5*nEven+k))
  allocate(blockDbleMem(8*nEven-1))

  d => blockDbleMem(1:nEven)
  e => blockDbleMem((nEven+1):(2*nEven-1))
  work => blockDbleMem((2*nEven):(7*nEven-1))
  evLocal => blockDbleMem((7*nEven):(8*nEven-1))

  i = 0
  d = (/ (((n-1-2*i) / 2.0d0)**2 * ctpw, i=0, nEven-1) /)

  ! 
  ! convert n and i to double before the multiplication, in response to 
  ! reported bug that occurs for large N 
  ! - wburr, 2012
  !
  e = (/ ((dble(i) * dble(n - i)) / 2.0d0, i = 1, nEven-1) /)

  is_evenN  = (modulo(n, 2) .eq. 0)
  
  ! ensure integer values are converted to double before mult to avoid overflow
  if(is_evenN)  then 
     d(nEven) = ((n+1-2*nEven)/2.0d0)**2 * ctpw + dble(nEven) * dble(nOdd)/2.0d0 
  else 
     e(nEven -1) =  e(nEven -1) * sr2
  end if
  
  cmach = 'S'
  abstol = 2.0d0*dlamch(cmach)
  
  ! set ldv (ldz) to 2*n to force 'dstein' to skip odd column which
  ! will be used for odd eigenvectors
  ! The matrix v is considered to be of a different shape in the 
  ! call to 'tridiagMatrixEigen2'
  !
  call tridiagMatrixEigen(nEven, evenK, d, e, v, 2*n, evlocal, &
       abstol, blockIntMem, work)
  
  if(.not. is_evenN) then
    do i = 1, k, 2
       v(nEven, i) =  v(nEven, i) * sr2
    end do
  end if
  
  do i = 1, k, 2
    v((nEven+1):n ,i) = v(nOdd:1:-1, i)
  end do
  
  ! reorder eigenvalues and eigenfunction columns for even
  evenK_2M1 = evenK * 2 - 1
  j = 1
  ev((/ (i, i=evenK_2M1, 1, -2) /)) = evLocal((/ (j, j=1, evenK, 1) /)) 
  v(:, (/ (i, i=1, evenK_2M1, 2) /)) =  & 
       v(:,(/ (j, j=evenK_2M1, 1, -2) /))
  
  if(k > 1)  then
     
     ! odd eigenfunctions (if any)
     ! similar procedure to the even functions
     ! if n is odd, nOdd < nEven
     ! ensure integer values are converted to double before mult to avoid overflow

     d(1:nOdd) = (/ (((n-1-2*i) / 2.0d0)**2 * ctpw, i=0, nOdd-1) /)
     
     ! convert n and i to double before the multiplication, in response to bug when n is large
     e(1:(nOdd-1)) = (/ ((dble(i) * dble(n - i)) / 2.0d0, i = 1, nOdd-1) /)
     
     if(is_evenN) then 
        ! convert n and i to double before the multiplication, in response to bug when n is large
        d(nOdd) = ((n+1-2*nEven)/2.0d0)**2 * ctpw - dble(nEven) * dble(nOdd)/2.0d0  
     end if
     
     call tridiagMatrixEigen(nOdd, oddK, d, e, v(1,2), 2*n, evlocal, &
          abstol, blockIntMem, work)
     
     if(.not. is_evenN) then
        do i= 2, k, 2
           v(nEven,i) = 0.0d0
        end do
     end if
     
     do i = 2, k, 2
        v((nEven+1):n ,i) = - v(nOdd:1:-1, i)
     end do
     
     ! reorder eigenvalues and eigenfunction columns for odd
     oddK_2 = oddK * 2
     ev((/ (i, i=oddK_2, 2, -2) /)) = evLocal((/ (j, j=1, oddK, 1) /))
     v(:, (/ (i, i=2, oddK_2, 2) /)) =  & 
          v(:,(/ (j, j=oddK_2, 2, -2) /))
     
  end if
  
  ! normalize
  iTest = 0
  if (mod(n,2) .eq. 0) then
     iTest = n/2 +1
  else 
     iTest = n/2 +2
  end if

  do j = 1 , k
     sqrtsumsq = dsqrt(sum( v(:,j)**2 ))
     ! set polarity to agree with Slepian (1978)
     ! differs from Percival and Walden's formulation
     ! dpss slope up at centre, this agrees with Thomson (1982)
     if(v(iTest,j) < 0.0d0) then
        sqrtsumsq = -1.0d0 * sqrtsumsq
     end if
     v(:,j) = v(:,j) / sqrtsumsq
  end do

  nullify(d,e,work,evLocal)
  deallocate(blockDbleMem)
  deallocate(blockIntMem)
  
end subroutine fdpss

subroutine tridiagMatrixEigen(n, k, d, e, v, ldv, ev, &
     abstol, blockIntMem, work)
  !assumes memory is allocated by the calling subroutine.
  implicit none
  
  character :: range, order, cmach
  integer :: nsplit, m, info
  integer :: k, ldv, n, il
  integer, target :: blockIntMem(5*n+k)
  integer,  pointer :: iblock(:), isplit(:), iwork(:), ifail(:)
  double precision :: vl, vu, abstol
  double precision :: d(n), e(n-1), v(ldv,k),  ev(n), work(5*n)

  range = 'I'
  order = 'E'
  cmach = 'S'

  il = n - k + 1
  m = k
     
  iblock => blockIntMem(1:n)
  isplit => blockIntMem((n+1):(2*n))
  iwork => blockIntMem((2*n+1):(5*n))
  ifail => blockIntMem((5*n+1):(5*n+k))

  call dstebz(range, order, n, vl, vu, il, n, &   
       abstol, d, e, m, nsplit, ev, iblock, isplit, &
       work, iwork, info)

  call dstein(n, d, e, m, ev, iblock, isplit, v, ldv, &
           work, iwork, ifail, info)

  nullify(iblock, isplit, iwork, ifail)

end subroutine tridiagMatrixEigen
