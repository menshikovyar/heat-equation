module mpi_sol
    use mpi
    implicit none
    
    contains

    subroutine find_d_c_f_r(A, b, n, D, C, F, r) 
        integer, intent(in) :: n
        real, intent(in) :: b(n), A(n,n)
        real, intent(out) :: F(n,n),r(n)
        real ::  D(n,n), C(n,n)
        integer :: i, j, process_rank, ierr
        
        call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierr)
        
        if (process_rank == 0) then   
        do j = 1, n
            do i = 1, n
                if(i == j) then
                    D(i, j) = 1/A(i, j)
                    C(i, j) = 0
                else
                    D(i, j) = 0
                    C(i, j) = -A(i, j) 
                end if
            end do
        end do
        F = matmul(D, C) 
        r = matmul(D, b)
        end if
    end subroutine find_d_c_f_r


    SUBROUTINE find_x(x0, x1, n,local_n, eps, F, r) 
    	integer, parameter :: k = 100000	
        integer :: n, ierr, iter, process_rank, i, j, local_n, size_cluster
        real, intent(inout) :: x1(n)
        real ::  x0(n), F(n,n), r(n)
        real :: norm, eps, time_start, time_end 
        real, allocatable :: local_F(:, :), local_x1(:), local_r(:)
        
        
        call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierr) 
        
        allocate(local_F(n, local_n), local_x1(local_n), local_r(local_n)) 
    	iter = 0 !задаем счетчики
    	norm = 1 
    	call MPI_SCATTER(transpose(F), local_n * n, MPI_REAL,local_F, local_n * n, MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
    	call MPI_SCATTER(r, local_n, MPI_REAL,local_r, local_n, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

    	call cpu_time(time_start)
    	do while(iter < k) 
    	    iter = iter + 1
        	do j = 1, local_n
        	    local_x1(j) = dot_product(local_F(:, j), x0) + local_r(j)
        	end do
    	    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  	    call MPI_ALLGATHER(local_x1, local_n, MPI_REAL, x1, local_n, MPI_REAL, MPI_COMM_WORLD, ierr) 
            call calc_norm(x1 - x0, n, norm) 
            if (norm < eps) then
                exit
            end if
            x0 = x1
    	end do

    	call cpu_time(time_end)
    	
    	deallocate(local_F, local_x1, local_r)
    	
    	
        
    end SUBROUTINE find_x

    subroutine calc_norm(vector, n, norm) 
        real, dimension(:), intent(in) :: vector
        integer, intent(in) :: n
        real, intent(out) :: norm
        integer :: i
        real :: sum

        sum = 0.0
        do i = 1, n
        sum = sum + vector(i)**2
        end do
        norm = sqrt(sum)
    end subroutine calc_norm
    
end module


