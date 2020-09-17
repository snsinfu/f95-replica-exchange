program main
    use rex

    implicit none

    real :: temps(3) = (/ 0.1, 1.0, 2.0 /)
    integer(8) :: seed
    integer :: rank
    integer :: replica_id
    integer :: step
    real :: energy
    logical :: ok

    real :: dt = 0.001
    real :: x = 0
    real :: v = 0
    real :: F = 0
    real :: kT

    seed = 1

    call rex_init(temps)
    call rex_seed(seed)
    call rex_rank(rank)
    call rex_id(replica_id)

    do step = 0, 10000
        ! Potential energy
        energy = (x + 3.5)**2 + sin(x - 2.5)
        kT = temps(replica_id)

        if (replica_id == 1 .and. mod(step, 100) == 0) then
            print "(i5,a,f8.3,a,f6.3)", step, " ", x, " ", energy
        end if

        if (mod(step, 100) == 0) then
            call rex_exchange(energy, ok)
            call rex_id(replica_id)
        end if

        v = v * sqrt(kT) / (abs(v) + 1e-6)

        x = x + dt * v
        F = -2 * (x + 3.5) - 2 * (x - 2.5)
        v = v + dt * F
    end do

    call rex_finalize()

end program
