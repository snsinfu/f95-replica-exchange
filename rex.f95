module rex
    implicit none

    public :: rex_init
    public :: rex_seed
    public :: rex_exchange
    public :: rex_id
    public :: rex_rank
    public :: rex_finalize

private
    include "mpif.h"

    ! MPI message tag used for replica exchange messages. This is an arbitrary
    ! number that is not used for any other messages.
    integer, parameter :: TAG_EXCHANGE = 1

    ! Maximum number of expected replicas (MPI processes).
    integer, parameter :: MAX_REPLICAS = 256

    ! Temperature assigned to each replica.
    real :: replica_temps(MAX_REPLICAS)

    ! Number of MPI processes.
    integer :: processes

    ! Fixed MPI rank of this process.
    integer :: self_rank

    ! Replica ID currently owned by this process. This nunmber changes as
    ! replica exchange simulation progresses.
    integer :: self_id

    ! Number of replica exchange attempts so far.
    integer :: exchange_attempts

    ! SFC64 pseudorandom number generator state.
    integer(8) :: sfc_state(4)

contains

    ! Initializes MPI and replica configurations.
    subroutine rex_init(temps)
        real temps(:)
        integer err

        call mpi_init(err)
        call mpi_comm_size(MPI_COMM_WORLD, processes, err)
        call mpi_comm_rank(MPI_COMM_WORLD, self_rank, err)

        ! Set the temperature of each replica as specified. XXX: The caller
        ! needs to provide temperature values without knowing the number of
        ! replicas. Bad API design!
        replica_temps(:processes) = temps(:processes)

        ! Set initial replica ID of this process (which will be altered upon
        ! subsequent replica exchange events). Note that the MPI rank is
        ! zero-based. So, +1 to make replica ID one-based.
        self_id = self_rank + 1
        exchange_attempts = 0

        return
    end subroutine

    ! Gets the current replica ID owned by this process.
    subroutine rex_id(id)
        integer, intent(out) :: id
        id = self_id
        return
    end subroutine

    ! Gets the MPI process rank.
    subroutine rex_rank(rank)
        integer, intent(out) :: rank
        rank = self_rank
        return
    end subroutine

    ! Attempts to replica-exchange state with adjacent process.
    subroutine rex_exchange(energy, ok)
        real, intent(in) :: energy
        logical, intent(out) :: ok
        integer peer_rank
        integer last_id

        if (mod(exchange_attempts, 2) == 0) then
            if (mod(self_rank, 2) == 0) then
                peer_rank = self_rank + 1
            else
                peer_rank = self_rank - 1
            end if
        else
            if (mod(self_rank, 2) == 0) then
                peer_rank = self_rank - 1
            else
                peer_rank = self_rank + 1
            end if
        end if

        last_id = self_id

        if (peer_rank >= 0 .and. peer_rank < processes) then
            if (self_rank > peer_rank) then
                call rex_exchange_as_leader(energy, peer_rank)
            else
                call rex_exchange_as_follower(energy, peer_rank)
            end if
        end if

        exchange_attempts = exchange_attempts + 1
        ok = self_id /= last_id

        return
    end subroutine

    subroutine rex_exchange_as_leader(self_energy, peer_rank)
        real, intent(in) :: self_energy
        integer, intent(in) :: peer_rank
        real self_temp, peer_temp, peer_energy
        integer peer_id, tmp_id
        real u_rand
        real metropolis
        real message(2)
        integer status_(MPI_STATUS_SIZE)
        integer err

        call mpi_recv(message, 2, MPI_REAL, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, status_, err)
        peer_id = int(message(1))
        peer_energy = message(2)

        self_temp = replica_temps(self_id)
        peer_temp = replica_temps(peer_id)

        call rex_rand_uniform(u_rand)
        metropolis = exp((1 / self_temp - 1 / peer_temp) * (self_energy - peer_energy))

        if (u_rand < metropolis) then
            tmp_id = self_id
            self_id = peer_id
            peer_id = tmp_id
        end if

        message(1) = real(peer_id)
        message(2) = peer_energy
        call mpi_send(message, 2, MPI_REAL, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, err)

        return
    end subroutine

    subroutine rex_exchange_as_follower(self_energy, peer_rank)
        real, intent(in) :: self_energy
        integer, intent(in) :: peer_rank
        real message(2)
        integer status_(MPI_STATUS_SIZE)
        integer err

        message(1) = real(self_id)
        message(2) = self_energy
        call mpi_send(message, 2, MPI_REAL, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, err)
        call mpi_recv(message, 2, MPI_REAL, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, status_, err)
        self_id = int(message(1))

        return
    end subroutine

    subroutine rex_finalize()
        integer err
        call mpi_finalize(err)
        return
    end subroutine

    subroutine rex_seed(seed)
        integer(8), intent(in) :: seed
        integer(8) :: discard
        integer :: i

        sfc_state(1:3) = seed
        sfc_state(4) = 0

        do i = 1, 12
            call rex_rand_bits(discard)
        end do

        return
    end subroutine

    subroutine rex_rand_bits(val)
        integer(8), intent(out) :: val
        integer(8) :: a, b, c, x

        a = sfc_state(1)
        b = sfc_state(2)
        c = sfc_state(3)
        x = sfc_state(4)

        val = a + b + x
        x = x + 1
        a = xor(b, ishft(b, 11))
        b = c + ishft(c, 3)
        c = val + or(ishft(c, 24), ishft(c, 24 - 64))

        sfc_state = (/ a, b, c, x /)

        return
    end subroutine

    subroutine rex_rand_uniform(val)
        real, intent(out) :: val
        integer(8) :: bits

        call rex_rand_bits(bits)
        val = 0.5 + real(bits) / (2.0 ** 64)

        return
    end subroutine
end module
