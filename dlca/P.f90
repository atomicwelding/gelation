program DLCA
  implicit none

  ! Parameters to tune
  integer, parameter :: Npts = 40 ! number of phis
  integer, parameter :: L = 10 ! size of the box
  integer, parameter :: MAX_STEPS_WITHOUT_AGGREGATION = 1e5  !!! Checking for no dynamics
  real, parameter :: alpha = 0.55, MAX_PHI = 0.3


  integer, dimension(0:L-1, 0:L-1, 0:L-1) :: grid
  type :: Particle
     integer :: id
     integer :: x, y, z
     integer :: m = 1
     integer :: cluster_id
  end type Particle
  type(Particle), allocatable :: particles(:)
  integer, allocatable :: cluster_size(:)

  integer :: tsim = 0, Nc = 0, N = 0
  real :: tgel = 0.0, tgel_report = 0.0
  logical :: hasEnded, hasPercolated

  integer :: run, runs = 100
  integer :: number_percolate, number_reach_max_aggreg = 0
  integer :: no_growth_steps, prev_max_cluster

  real, dimension(0:Npts-1) :: phis
  integer phi

  integer :: fd, io
  character(len=100) :: errmsg
  character(len=*), parameter :: filepath = "P.dat"

  ! Initialize phis from 0 to MAX_PHI in Npts steps -> Figure 2 of Gimel.
  do phi = 0, Npts-1
     phis(phi) = real(phi) * MAX_PHI / real(Npts-1)
  end do

  open(fd, file=filepath, status="replace", iostat=io, iomsg=errmsg)
  if(io /= 0) then
     print *, "Error while opening file : ", trim(errmsg)
     stop 1
  end if

  write(fd, *) "@ L = ", L
  call flush(fd)
  write(fd, *) "@ MAX_STEPS_WITHOUT_AGGREGATION = ", MAX_STEPS_WITHOUT_AGGREGATION
  write(fd, *) "@ phi0, P, tsim_avg, tgel_avg"
  call flush(fd)
  do phi=0, Npts-1
     N = int(phis(phi) * L**3)
     tsim = 0
     tgel_report = 0.0

     print *, "Running for ", N, " particles"
     allocate(particles(0:N-1))
     allocate(cluster_size(0:N-1))
     grid = -1
     hasEnded = .false.
     number_percolate = 0


     ! First stopping criterion
     !!! ARGUMENT TO FORCE CONVERGENCE : IF N < L, THEN THE SYSTEM CANNOT PERCOLATE.
     !!! AS SUCH, WE CAN DIRECTLY GO TO NEXT ITERATION 
     if(N < L) goto 82

     do run=0,runs-1
        
        tgel = 0.0
        hasEnded = .false.
        hasPercolated = .false.
        grid = -1
        no_growth_steps = 0
        prev_max_cluster = 1
        call init_particles

        ! Running the simulation 
        do while(hasEnded .NEQV. .true.)
           tsim = tsim + 1
           Nc = count_active_clusters()
           tgel = tgel + 1.0 / Nc
           call trial

           ! 2nd stopping criterion
           if (maxval(cluster_size) <= prev_max_cluster) then
              no_growth_steps = no_growth_steps + 1
           else
              no_growth_steps = 0
              prev_max_cluster = maxval(cluster_size)
           end if
           if (no_growth_steps >= MAX_STEPS_WITHOUT_AGGREGATION) then
              number_reach_max_aggreg = number_reach_max_aggreg + 1
              hasEnded = .true.
           end if
           
        end do
        
        if(hasPercolated) then
           number_percolate = number_percolate + 1
           tgel_report = tgel_report + tgel
        end if
        
     end do

82   print *, "Number of particles", N
     print *, "Number of percolated systems ", number_percolate
     print *, "Number of runs", runs
     print *, "Number of runs that reached MAX_STEPS_WITHOUT_AGGREGATION", number_reach_max_aggreg
     print *
     write(fd, *) phis(phi), dble(number_percolate)/dble(runs), dble(tsim)/dble(runs), dble(tgel_report)/dble(number_percolate)
     call flush(fd)

     deallocate(particles)
     deallocate(cluster_size)
  end do

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
            hasPercolated = .true.
            hasEnded = .true.
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

subroutine update_cluster_sizes
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
    d = 1
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
    prob = real(cluster_size(cid)) ** (-1.0 * alpha)
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
end program DLCA
