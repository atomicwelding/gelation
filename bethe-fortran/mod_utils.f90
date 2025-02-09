module utils
  use iso_fortran_env, only: real32
  implicit none
contains

  !! IO related
  subroutine display(size, matrix)
    ! displaying an integer matrix in a compact format
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit
    integer, intent(in) :: size
    real(kind=real32), dimension(size, size), intent(in) :: matrix

    integer :: i, j
    do i=1,size
       do j=1,size
          write(stdout, "(F5.1)", advance="no") matrix(i,j)
       end do
       write(stdout,*) ! new line
    end do
  end subroutine display


  subroutine loadfile(size, matrix, filepath)
    integer, intent(in) :: size
    real(kind=real32), dimension(size, size), intent(inout) :: matrix
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
    real(kind=real32), dimension(size, size), intent(in) :: matrix
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
          write(io, "(F5.1)", advance="no") matrix(i,j)
       end do
       write(io,*)
    end do
  end subroutine writefile
end module utils
