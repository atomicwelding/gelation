module math
  use iso_fortran_env, only: real32
  implicit none

  private
  
  public :: slpl
  public :: sid
  

contains
  ! single precision, procedures acting on square matrices only

  subroutine slpl(N, A)
    ! Compute in-place the laplacian of a matrix A of dimension N
    integer, intent(in) :: N
    real(kind=real32), dimension(N,N), intent(inout) :: A

    real(kind=real32), dimension(N, N) :: DEG ! Degree matrix
    integer :: i

    DEG = 0.0_real32
    do i=1,N
       DEG(i,i) = sum(A(i,:))
    end do

    A = DEG - A
  end subroutine slpl
  
  subroutine sid(N, A)
    ! Fill a matrix A  as an identity matrix of dimension N
    integer, intent(in) :: N
    real(kind=real32), dimension(N,N), intent(inout) :: A
    integer :: i

    A = 0.0_real32
    do i=1,N
       A(i,i) = 1.0_real32
    end do
  end subroutine sid

  subroutine spinv(N, A, A_PINV)
    ! Compute Moore-Penrose's  pseudo-inverse of a matrix A of dimension N
    ! Stores the result in A_PINV
    integer, intent(in) :: N
    real(kind=real32), dimension(N,N), intent(in) :: A
    real(kind=real32), dimension(N,N), intent(out) :: A_PINV

    A_PINV = 0.0_real32

    call handle_lapack_error(1, "hello, world!")
  end subroutine spinv


  ! helper
  subroutine handle_lapack_error(info, msg)
    integer, intent(in) :: info
    character(len=*), intent(in) :: msg

    if(info /= 0) then
       print *, " INFO = ", info
       print *, "ERROR : ", msg
       stop
    end if
  end subroutine handle_lapack_error
  
end module math
