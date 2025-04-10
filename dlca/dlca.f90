program DLCA

    ! TODO
    !! Reread the article to determine the parameters to compute in order to reproduce the figures.
    !! The simulation ends when Nc =1 or when a cluster joins two opposite sides of the system.

    !! Quand on finit la simu : on note snapshot du début et de fin
    !! on note temps de gel
    !! on note clusters
    implicit none

    logical, parameter :: write_trajectories = .true.
    integer, parameter :: freq_write = 10, freq_print = 5

    integer, parameter :: N = 100, L = 20
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
    
    integer :: tsim = 0, Nc = 0
    real :: tgel = 0
    logical :: hasEnded, hasPercolated

    integer :: fd = 10, io
    character(len=*), parameter :: filepath = "data.xyz"
    character(len=100) :: errmsg

    ! ==== INIT ==== !
    tsim = 0
    Nc = count_active_clusters()
    hasEnded = .false.
    grid = -1
    call init_particles

    ! ==== RUNNING THE SIMULATION ==== !
    open(fd, file=filepath, status="replace", iostat=io, iomsg=errmsg)
    if(io /= 0) then
        print *, "Error while opening file : ", trim(errmsg)
        stop 1
    end if

    ! write first frame
    call write_frame(0)
    ! dynamics
    do while(hasEnded .NEQV. .true.)
        
        tsim = tsim + 1
        Nc = count_active_clusters()
        tgel = tgel + 1.0 / Nc
    
        call trial

        if(modulo(tsim, freq_write) == 0 .AND. write_trajectories) call write_frame(tsim)
        if(modulo(tsim, freq_print) == 0) then
            print *, tsim, ": ", tgel
        end if 
    end do
    ! write last frame
    call write_frame(tsim+1)
    close(fd)

contains
    function count_active_clusters()
        integer :: count_active_clusters, i
        count_active_clusters = 0
        do i = 0, N-1
            if (cluster_size(i) > 0) count_active_clusters = count_active_clusters + 1
        end do
    end function count_active_clusters

    subroutine checkEndSimulation
        integer :: i, cid
        type FaceTouch
            logical :: min_x = .false.
            logical :: max_x = .false.
            logical :: min_y = .false.
            logical :: max_y = .false.
            logical :: min_z = .false.
            logical :: max_z = .false.
        end type FaceTouch
        type(FaceTouch), dimension(0:N-1) :: touches

        ! reset
        hasPercolated = .false.

        ! if all atoms are in one cluster, then it's over
        if (any(cluster_size == N)) then
            hasEnded = .true.
            return
        end if

        do i = 0, N-1
            cid = particles(i)%cluster_id
            if (particles(i)%x == 0)    touches(cid)%min_x = .true.
            if (particles(i)%x == L-1)  touches(cid)%max_x = .true.
            if (particles(i)%y == 0)    touches(cid)%min_y = .true.
            if (particles(i)%y == L-1)  touches(cid)%max_y = .true.
            if (particles(i)%z == 0)    touches(cid)%min_z = .true.
            if (particles(i)%z == L-1)  touches(cid)%max_z = .true.
        end do

        do cid = 0, N-1
            if ((touches(cid)%min_x .and. touches(cid)%max_x) .or. &
                (touches(cid)%min_y .and. touches(cid)%max_y) .or. &
                (touches(cid)%min_z .and. touches(cid)%max_z)) then
                print *, "Cluster", cid, " touches:"
                print *, "  X:", touches(cid)%min_x, touches(cid)%max_x
                print *, "  Y:", touches(cid)%min_y, touches(cid)%max_y
                print *, "  Z:", touches(cid)%min_z, touches(cid)%max_z
                hasPercolated = .true.
                hasEnded = .true.
                print *, "Terminated: cluster ", cid, " has percolated."
                call write_percolating_cluster(tsim, cid)
                call write_frame(tsim)
                return
            end if
        end do
    end subroutine checkEndSimulation

    subroutine flood_fill_all()
        integer :: i, cid
        logical, dimension(0:N-1) :: visited
        integer, dimension(:), allocatable :: queue
        integer :: qstart, qend
        integer :: current_id, neighbor_id
        integer :: x, y, z, nx, ny, nz, j
        integer, dimension(6,3) :: dirs

        dirs = reshape([1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], [6,3])
        visited = .false.
        cid = 0

        allocate(queue(N))

        do i = 0, N-1
            if (.not. visited(i)) then
                qstart = 1
                qend = 1
                queue(qstart) = i
                visited(i) = .true.
                particles(i)%cluster_id = cid

                do while (qstart <= qend)
                    current_id = queue(qstart)
                    qstart = qstart + 1

                    x = particles(current_id)%x
                    y = particles(current_id)%y
                    z = particles(current_id)%z

                    do j = 1, 6
                        nx = x + dirs(j,1)
                        ny = y + dirs(j,2)
                        nz = z + dirs(j,3)

                        if (nx < 0 .or. nx >= L .or. &
                            ny < 0 .or. ny >= L .or. &
                            nz < 0 .or. nz >= L) cycle

                        neighbor_id = grid(nx, ny, nz)
                        if (neighbor_id /= -1 .and. .not. visited(neighbor_id)) then
                            visited(neighbor_id) = .true.
                            particles(neighbor_id)%cluster_id = cid
                            qend = qend + 1
                            queue(qend) = neighbor_id
                        end if
                    end do
                end do

                cid = cid + 1
            end if
        end do

        deallocate(queue)
    end subroutine flood_fill_all

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

    subroutine update_cluster_sizes()
        integer :: i
        cluster_size = 0
        do i = 0, N-1
            cluster_size(particles(i)%cluster_id) = cluster_size(particles(i)%cluster_id) + 1
        end do
    end subroutine update_cluster_sizes

    subroutine trial
        integer :: d, cid
        integer :: dx, dy, dz
        integer :: i
        integer :: x_old, y_old, z_old
        integer :: x_new, y_new, z_new
        real :: r, prob
        logical :: can_move
        integer, allocatable :: active_clusters(:)
        integer :: num_clusters, cluster_idx, j

        ! choose a cluster
        ! count clusters
        num_clusters = count(cluster_size > 0)
        allocate(active_clusters(num_clusters))

        j = 1
        do i = 0, N-1
            if (cluster_size(i) > 0) then
                active_clusters(j) = i
                j = j + 1
            end if
        end do

        call random_number(r)
        cluster_idx = int(r * num_clusters)
        if (cluster_idx >= num_clusters) cluster_idx = num_clusters - 1
        cid = active_clusters(cluster_idx + 1)

        deallocate(active_clusters)

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
                x_new = particles(i)%x + dx
                y_new = particles(i)%y + dy
                z_new = particles(i)%z + dz

                if (x_new < 0 .or. x_new >= L .or. &
                    y_new < 0 .or. y_new >= L .or. &
                    z_new < 0 .or. z_new >= L) then
                    can_move = .false.
                    exit
                end if

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

                x_new = x_old + dx
                y_new = y_old + dy
                z_new = z_old + dz

                particles(i)%x = x_new
                particles(i)%y = y_new
                particles(i)%z = z_new
                grid(x_new, y_new, z_new) = particles(i)%id

            end if
        end do

        ! perform flood fill once after the entire cluster has moved
        call flood_fill_all()
        call update_cluster_sizes()
        call checkEndSimulation
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

    subroutine write_percolating_cluster(current_frame, cid)
        integer, intent(in) :: current_frame, cid
        integer :: i, fdc, iop
        character(len=100) :: errmsgp
        character(len=*), parameter :: path_percolating = "percolating_cluster.xyz"
        character :: symbol

        open(unit=fdc, file=path_percolating, status="replace", iostat=iop, iomsg=errmsgp)
        if (iop /= 0) then
            print *, "Error while opening percolating cluster file: ", trim(errmsgp)
            return
        end if

        ! Count the number of particles in the cluster
        write(fdc, *) cluster_size(cid)
        write(fdc, *) "Percolating Cluster at Frame ", current_frame

        symbol = 'O'
        do i = 0, N-1
            if (particles(i)%cluster_id == cid) then
                write(fdc, *) symbol, real(particles(i)%x), real(particles(i)%y), real(particles(i)%z)
            end if
        end do

        close(fdc)
    end subroutine write_percolating_cluster
end program DLCA