module utils
  use iso_fortran_env, only: int8
  implicit none
  !! TODO
  !!! Calculer la distance de resistance
contains

  !! IO related
  subroutine display(size, matrix)
    ! displaying an integer matrix in a compact format
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit
    integer, intent(in) :: size
    integer(kind=int8), dimension(size, size), intent(in) :: matrix

    integer :: i, j
    do i=1,size
       do j=1,size
          write(stdout, "(I3)", advance="no") matrix(i,j)
       end do
       write(stdout,*) ! new line
    end do
  end subroutine display


  subroutine loadfile(size, matrix, filepath)
    integer, intent(in) :: size
    integer(kind=int8), dimension(size, size), intent(inout) :: matrix
    character(len=*), intent(in) :: filepath

    integer :: i, j

    integer :: io, iostat
    logical :: exists

    inquire(file=filepath, exist=exists)
    if(.not. exists) then
       print *, "Error : File ", trim(filepath), " does not exist."
       STOP 1
    end if


    print *, "Loading file ..."
    open(newunit=io, file=filepath)
    read(io, *, iostat = iostat) matrix
    if(iostat /= 0) then
       print *, "Error reading matrix from file '", trim(filepath), "'. IOSTAT = ", iostat
       close(io)
       stop iostat
    end if


    close(io)
  end subroutine loadfile

  subroutine writefile(size, matrix, filepath)
    integer, intent(in) :: size
    integer(kind=int8), dimension(size, size), intent(in) :: matrix
    character(len=*), intent(inout) :: filepath

    integer :: io, i, j, hash
    logical :: exists
    real :: r


    do
       inquire(file=trim(filepath), exist = exists)
       if(.not. exists) exit

       print *, "File ", trim(filepath), " already exists."
       print *, "Generating a new filename ..."

       call random_number(r)
       hash = int(r * 1000)
       write(filepath, "(A,I3,A)") "default", hash, ".txt"
    end do


    print *, "Writing to file: ", trim(filepath)

    open(newunit=io, file=trim(filepath))
    ! writing in a  compact format
    do i=1,size
       do j=1,size
          write(io, "(I2)", advance="no") matrix(i,j)
       end do
       write(io,*)
    end do
  end subroutine writefile

  !! Maths related
  pure function identity(dim)
    integer, intent(in) :: dim
    integer(kind=int8), dimension(dim, dim) :: identity
    integer :: i

    identity = 0
    do i=1,dim
       identity(i,i) = 1
    end do
  end function identity

  pure function laplacian(dim, adj_matrix)
    integer, intent(in) :: dim
    integer(kind=int8), dimension(dim, dim), intent(in) :: adj_matrix

    ! we assume such a low kind because polymers won't have a very high functionality
    integer(kind=int8), dimension(dim, dim) :: laplacian, degree_matrix
    integer :: i
    laplacian = 0

    do i=1,dim
       degree_matrix(i,i) = sum(adj_matrix(i,:))
    end do
    
    laplacian =  degree_matrix - adj_matrix   
  end function laplacian

end module utils
