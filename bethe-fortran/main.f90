program Main
  use iso_fortran_env, only: int8
  use utils, only : display, writefile, loadfile, &
       laplacian

  implicit none
  integer, parameter :: Nmonomers = 10000, functionality = 3
  real, parameter :: prob = 0.8

  integer(kind=int8), dimension(Nmonomers, Nmonomers) :: lattice, lapl

  ! init the lattice
  lattice = 0
  call rand_adj(Nmonomers, lattice, prob, functionality)

  ! compute the laplacian and display
  lapl = laplacian(Nmonomers, lattice)
  call display(Nmonomers, lapl)

contains
  subroutine rand_adj(size, matrix, p, max_degree)
    integer, intent(in) :: size, max_degree
    real, intent(in) :: p
    integer(kind=int8), dimension(size, size), intent(inout) :: matrix

    integer :: i, j
    integer, dimension(size) :: deg
    real :: r

    call random_seed()

    deg = 0
    do i=1, size - 1
       if(deg(i) >= max_degree) cycle
       do j=i+1, size
          if(deg(j) >= max_degree) cycle

          call random_number(r)
          if(r < p) then
             matrix(i,j) = 1
             matrix(j,i) = 1

             deg(i) = deg(i) + 1
             deg(j) = deg(j) + 1

             if(deg(i) >= max_degree) exit
          end if
          
       end do
    end do
  end subroutine rand_adj
end program Main
