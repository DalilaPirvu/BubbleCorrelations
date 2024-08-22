!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS

#include "macros.h"
#define SMOOTH 1

program Gross_Pitaevskii_1d
  ! the following modules are incorporated:
  use, intrinsic :: iso_c_binding
  use gaussianRandomField
  use integrator
  use constants
  use eom
  implicit none

  real(dl), pointer :: time
  integer :: nTime = 1
  integer :: sim, nSims = 5000, minSim = 0
  integer :: spec = nyq/16
  real(dl), parameter :: alph = 16._dl
  real(dl), parameter :: phi0 = twopi / 7._dl
  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call initialize_rand(93286123,12)
  call setup(nVar)

  do sim = 0, nSims-1 ! run nSims simulations for the parameters, each with its output files
      call initialise_fields(fld, nyq, phi0, spec)
      if (sim >= minSim) then
          call time_evolve(dx/alph, nLat, sim)
          print*, "Simulation ", sim+1, " in ", nSims , " done!"
      endif
  enddo
  print*, 'm2eff = ', m2eff, 'len = ', len, 'dt = ', dx/alph, 'dx = ', dx, 'dk = ', dk, 'phi0 = ', phi0, 'alph = ', alph, 'spec = ', spec

contains

  subroutine initialise_fields(fld,kmax,phi,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax
    real(dl), intent(in), optional :: phi
    integer, intent(in), optional :: klat
    integer :: kc
    real(dl) :: phiL

    kc = nLat/2+1; if (present(klat)) kc = klat
    phiL = 0.5_dl*twopi; if (present(phi)) phiL = phi

    ! initialize mean fields
    fld(:,1) = 0._dl!5_dl*twopi
    fld(:,2) = 0._dl

    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here
    call initialize_vacuum_fluctuations(fld,len,m2eff,kmax,phiL,kc) ! Change this call as necessary
  end subroutine initialise_fields

!!!!!!!!!!!!!!!!!!
! Time Evolution !
!!!!!!!!!!!!!!!!!!

  subroutine time_evolve(dt, ns, sim)
    real(dl) :: dt, dtout; integer :: ns, i, j, outsize, sim
    if (dt > dx) print*, "Warning, violating Courant condition" !i.e. alph > 1
    outsize = ns/nTime ! how much time is needed to compute ncross field intersections or sth like that
    dtout = dt*outsize

    call output_fields(fld, dt, dtout, sim)
    do i = 1, nTime-1 ! this loops over time slices
       do j = 1, outsize ! how many integrations are needed to evolve by one time slice; depends on how many crosses
          call gl10(yvec,dt)
       enddo
       call output_fields(fld, dt, dtout, sim)
    enddo
  end subroutine time_evolve

!!!!!!!!!!!!!!!!
! Fluctuations !
!!!!!!!!!!!!!!!!

  subroutine initialize_vacuum_fluctuations(fld,len,m2,kspec,phi,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len, m2
    integer, intent(in), optional :: kspec, klat
    real(dl), intent(in), optional :: phi

    real(dl), dimension(1:size(fld(:,1)/2+1)) :: spec, w2eff  ! remove w2eff here, it's unneeded
    real(dl), dimension(1:size(fld(:,1))) :: df
    integer :: i,km,kc
    real(dl) :: phiL, norm

    integer :: n, nn; real(dl) :: dk
    dk = twopi / len; n = size(fld(:,1)); nn = n/2+1

    km = size(spec); if (present(kspec)) km = kspec
    kc = size(spec); if (present(klat))  kc = klat

    phiL = 0.5_dl*twopi; if (present(phi)) phiL = phi

    norm = (0.5_dl)**0.5 / phiL / sqrt(2._dl * len)

    do i=1,nn
       w2eff(i) = m2 + dk**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:km) = norm / w2eff(2:km)**0.25
    call generate_1dGRF(df,spec(:kc))
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kc))
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_vacuum_fluctuations

!!!!!!!!!!
! Output !
!!!!!!!!!! 

  subroutine setup(nVar)
    integer, intent(in) :: nVar
    call init_integrator(nVar)
    call initialize_transform_1d(tPair,nLat)
  end subroutine setup

  character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

  character(len=20) function real_str(k)
    real(dl), intent(in) :: k
    write (real_str, '(f12.4)') k
    real_str = adjustl(real_str)
  end function real_str


  subroutine output_fields(fld, dt, dtout, sim)
    real(dl), dimension(1:nLat, 1:2) :: fld
    real(dl) :: dt, dtout
    logical :: o; integer :: m, sim
    integer, parameter :: oFile = 98
    inquire(file='/gpfs/dpirvu/bubble_correlations/free_x'//trim(str(nLat))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lambda))//'_spec'//trim(str(spec))//'_sim'//trim(str(sim))//'_fields.dat', opened=o)
    if (.not.o) then
       open(unit=oFile,file='/gpfs/dpirvu/bubble_correlations/free_x'//trim(str(nLat))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lambda))//'_spec'//trim(str(spec))//'_sim'//trim(str(sim))//'_fields.dat')
       write(oFile,*) "# Lattice Parameters dx = ", dx, 'len = ', len
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters m2eff = ", m2eff, "dk = ", dk
    endif

    do m = 1, nLat
       write(oFile,*) fld(m,1) !:), 0.5_dl*fld(m,2)**2._dl + 0.5_dl*gsq(m) + 4._dl*nu*( - cos(fld(m,1)) + 0.5_dl*lambda**2._dl * sin(fld(m,1))**2._dl )
    enddo
  end subroutine output_fields

end program Gross_Pitaevskii_1d
