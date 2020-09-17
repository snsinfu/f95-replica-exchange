program main
    use rex
    implicit none

    real :: temps(2) = (/ 1.0, 1.2 /)
    integer(8) :: seed
    integer :: rank
    integer :: replica_id
    integer :: i
    real :: energy
    logical :: ok

    seed = 1

    call rex_init(temps)
    call rex_seed(seed)
    call rex_rank(rank)

    do i = 1, 100
        call random_number(energy)

        call rex_exchange(energy, ok)
        call rex_id(replica_id)
        print "(i1,a,i3,a,i1)", rank, ":", i, " -> ", replica_id
    end do

    call rex_finalize()

end program
