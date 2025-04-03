program DLCA

    ! TODO
    !! Reread the article to determine the parameters to compute in order to reproduce the figures.
    implicit none

    integer, parameter :: N = 10000, L = 50, NUMBER_TRIALS = 100000
    real, parameter :: alpha = -0.55
    
    integer, dimension(0:L-1, 0:L-1, 0:L-1) :: grid
    type :: Particle
        integer :: id
        integer :: x, y, z
        integer :: m = 1
        integer :: cluster_id
    end type Particle
    type(Particle), dimension(0:N-1) :: particles
    integer, dimension(0:N-1) :: cluster_size

    integer :: t, fd = 10, io
    logical, parameter :: write_only_first_last = .false.
    integer, parameter :: freq_write = 1000
    character(len=*), parameter :: filepath = "data.xyz"
    character(len=100) :: errmsg

    ! ==== INIT ==== ! 
    grid = -1
    call init_particles

    ! ==== RUNNING THE SIMULATION ==== !
    open(fd, file=filepath, status="replace", iostat=io, iomsg=errmsg)
    if(io /= 0) then
        print *, "Error while opening file : ", trim(errmsg)
        stop 1
    end if

    do t=0,NUMBER_TRIALS-1
        call trial
        if(write_only_first_last .and. (t == 0 .or. t == NUMBER_TRIALS-1)) call write_frame(t)
        if (.not. write_only_first_last .and. modulo(t, freq_write) == 0) call write_frame(t)
    end do 
    close(fd)

contains
    subroutine init_particles
        integer :: i, x, y, z
        real :: r

        do i=0,N-1
            do
                call random_number(r)
                x = int(r * L)
                call random_number(r)
                y = int(r * L)
                call random_number(r)
                z = int(r * L)

                if (grid(x, y, z) == -1) exit
            end do

            particles(i)%id = i
            particles(i)%x = x 
            particles(i)%y = y
            particles(i)%z = z
            particles(i)%cluster_id = i
            cluster_size(i) = 1
            grid(x, y, z) = i
        end do 
    end subroutine init_particles

    subroutine flood_fill(id)
        integer, intent(in) :: id
        integer :: x, y, z, nx, ny, nz
        integer :: current_id, neighbor_id
        integer, dimension(:), allocatable :: queue
        integer :: qstart, qend, i
        integer, dimension(6,3) :: dirs

        dirs = reshape([1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], [6,3])
        allocate(queue(N))
        qstart = 1
        qend = 1
        queue(1) = id

        do while (qstart <= qend)
            current_id = queue(qstart)
            qstart = qstart + 1

            x = particles(current_id)%x
            y = particles(current_id)%y
            z = particles(current_id)%z

            do i = 1, 6
                nx = modulo(x + dirs(i,1),L)
                ny = modulo(y + dirs(i,2),L)
                nz = modulo(z + dirs(i,3),L)

                neighbor_id = grid(nx, ny, nz)
                if (neighbor_id /= -1) then
                    if (particles(neighbor_id)%cluster_id /= particles(id)%cluster_id) then
                        particles(neighbor_id)%cluster_id = particles(id)%cluster_id
                        qend = qend + 1
                        queue(qend) = neighbor_id
                    end if
                end if
            end do
        end do

        deallocate(queue)
    end subroutine flood_fill

    subroutine update_cluster_sizes()
        integer :: i
        cluster_size = 0
        do i = 0, N-1
            cluster_size(particles(i)%cluster_id) = cluster_size(particles(i)%cluster_id) + 1
        end do
    end subroutine update_cluster_sizes

    subroutine trial
        integer :: d, id, cid
        integer :: dx, dy, dz
        integer :: i
        integer :: x_old, y_old, z_old
        integer :: x_new, y_new, z_new
        real :: r, prob
        logical :: can_move

        ! choose a particle
        call random_number(r)
        id = int(r * N)
        cid = particles(id)%cluster_id

        ! choose sign of displacement
        d = 2
        call random_number(r)
        if (r <= 0.5) d = -1 * d

        dx = 0
        dy = 0
        dz = 0
        ! choose coordinates to move
        call random_number(r)
        select case ( int(r*3) )
        case (0)
            dx = d
        case (1)
            dy = d
        case (2)
            dz = d
        end select

        ! check if all particles in the cluster can move
        can_move = .true.
        do i = 0, N-1
            if (particles(i)%cluster_id == cid) then
                x_new = modulo(particles(i)%x + dx, L)
                y_new = modulo(particles(i)%y + dy, L)
                z_new = modulo(particles(i)%z + dz, L)
                if (grid(x_new, y_new, z_new) /= -1) then
                    ! if the target cell is occupied, check if the occupant belongs to the same cluster
                    ! if yes, then keep can_move as true as the occupant will also move.
                    if (particles(grid(x_new, y_new, z_new))%cluster_id /= cid) then
                        can_move = .false.
                        exit
                    end if
                end if
            end if
        end do

        if (.not. can_move) return

        ! check acceptance
        prob = real(cluster_size(cid)) ** alpha
        call random_number(r)
        if (r >= prob) return

        ! move the whole cluster
        do i = 0, N-1
            if (particles(i)%cluster_id == cid) then
                x_old = particles(i)%x
                y_old = particles(i)%y
                z_old = particles(i)%z
                grid(x_old, y_old, z_old) = -1

                x_new = modulo(x_old + dx, L)
                y_new = modulo(y_old + dy, L)
                z_new = modulo(z_old + dz, L)

                particles(i)%x = x_new
                particles(i)%y = y_new
                particles(i)%z = z_new
                grid(x_new, y_new, z_new) = particles(i)%id

            end if
        end do

        ! perform flood fill once after the entire cluster has moved
        call flood_fill(id)
        call update_cluster_sizes()
    end subroutine trial


    subroutine write_frame(current_frame)
        integer, intent(in) :: current_frame
        integer :: i
        character, parameter :: symbol = 'N'

        write(fd, *) N
        write(fd, *) "Frame ", current_frame
        do i = 0, N-1
            write(fd, *) symbol, real(particles(i)%x), real(particles(i)%y), real(particles(i)%z)
        end do
    end subroutine write_frame

    
    ! subroutine write_frame_in_construction(current_frame)
    !     integer, intent(in) :: current_frame
    !     integer :: i, largest_cluster, largest_size
    !     character :: symbol

    !     ! determine the largest cluster
    !     largest_cluster = -1
    !     largest_size = 0
    !     do i = 0, N-1
    !         if (cluster_size(i) > largest_size) then
    !             largest_size = cluster_size(i)
    !             largest_cluster = i
    !         end if
    !     end do

    !     print *, largest_cluster, " => ", cluster_size(largest_cluster), " particles"

    !     ! write all particles, with red for the largest cluster (symbol 'O'), blue otherwise (symbol 'N')
    !     write(fd, *) N
    !     write(fd, *) "Frame ", current_frame
    !     do i = 0, N-1
    !         if (particles(i)%cluster_id == largest_cluster) then
    !             symbol = 'O'
    !         else
    !             symbol = 'N'
    !         end if
    !         write(fd, *) symbol, real(particles(i)%x), real(particles(i)%y), real(particles(i)%z)
    !     end do
    ! end subroutine write_frame_in_construction
end program DLCA