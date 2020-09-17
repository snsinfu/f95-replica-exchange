program main
    use rex

    implicit none

    real :: temps(5) = (/ 0.1, 0.3, 0.8, 1.2, 1.7 /)
    integer(8) :: seed
    integer :: rank
    integer :: replica_id
    integer :: step
    logical :: ok

    real :: dt = 0.005
    real :: x = 0
    real :: v = 0
    real :: F = 0
    real :: R = 0
    real :: energy
    real :: kT
    real :: friction
    real :: damping
    real :: sigma

    seed = 1
    friction = 0.05

    call rex_init(temps)
    call rex_seed(seed)
    call rex_rank(rank)
    call rex_id(replica_id)

    do step = 0, 100000
        ! Potential energy
        energy = system_energy(x, v)
        kT = temps(replica_id)

        if (replica_id == 1 .and. mod(step, 10) == 0) then
            print "(i6,a,f8.3,a,f6.3)", step, " ", x, " ", energy
        end if

        if (mod(step, 100) == 0) then
            call rex_exchange(energy, ok)
            call rex_id(replica_id)
        end if

        damping = exp(-friction * dt)
        sigma = sqrt(kT * (1 - exp(-2 * friction * dt)))
        call random_normal(R)

        v = v + dt * F / 2
        x = x + dt * v / 2
        v = v * damping
        v = v + sigma * R
        x = x + dt * v / 2
        F = system_force(x)
        v = v + dt * F / 2
    end do

    call rex_finalize()

    stop

contains
    real function system_energy(x, v)
        real, intent(in) :: x
        real, intent(in) :: v

        system_energy = 0.5 * (x**2 + sin(5 * x)) + v * v / 2
    end function

    real function system_force(x)
        real, intent(in) :: x

        system_force = 0.5 * (-2 * x - 5 * cos(5 * x))
    end function

    subroutine random_normal(z)
        real, intent(out) :: z
        real :: u, v

        do
            call random_number(u)
            call random_number(v)

            if (u > 1e-6) then
                exit
            end if
        end do

        z = sqrt(-2 * log(u)) * cos(6.2831853 * v)

        return
    end subroutine
end program
