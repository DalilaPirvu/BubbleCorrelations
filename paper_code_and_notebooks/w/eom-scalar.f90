!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
  use fftw3
  implicit none

  integer, parameter :: nFld = 1, nLat = 8192, nyq = nLat/2+1, nVar = 2*nFld*nLat+1
  real(dl), dimension(1:nVar), target :: yvec
  real(dl), parameter :: nu = 2.e-3
  real(dl), parameter :: omega = 0.25*50._dl*2._dl*nu**0.5
  real(dl), parameter :: del = (nu/2._dl)**0.5*(1._dl+0.5_dl)
  real(dl), parameter :: rho = 200._dl*2._dl*(nu)**0.5*2.**(-3)
  real(dl), parameter :: lambda = del*(2._dl/nu)**0.5
  real(dl), parameter :: m2eff = 4._dl*nu*(-1._dl+lambda**2)
  real(dl), parameter :: len = 8*0.5*50._dl / (2.*nu)**0.5
  real(dl), parameter :: dx = len/dble(nLat), dk = twopi/len
  type(transformPair1D) :: tPair

contains
  
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -4._dl*nu*( sin(yc(FLD)) + 0.5_dl*lambda**2*sin(2._dl*yc(FLD)) )
!    yp(DFLD) = - m2eff*yc(FLD)

    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
  end subroutine derivs

end module eom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS
