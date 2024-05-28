program main
    use mpi
    implicit none

    integer, parameter :: n = 100
    real, parameter :: eps = 0.001
    real, parameter :: tau = 1.0/n, h = 1.0/n ! Шаги по времени и пространству 
    real :: U_et(n, n), F(n, n), A(n,n), B(n,n), Ux(n, n), Uxt(n), Bt(n)
    real :: D(n,n), C(n,n)
    real :: x_curr, t_curr, start, finish
    integer :: i, j, ierr, size_cluster, process_Rank, local_n

    interface
        subroutine solve(A, b, x0, n,local_n, eps) !подпрограмма, которая вызывается из другого файла для решения СЛАУ методом Якоби
            integer :: n, ierror,local_n
            real :: A(n, n), b(n), x1(n), D(n,n), C(n,n), F(n,n), eps, r(n)
            real, intent(inout) :: x0(n)
        end subroutine solve

        subroutine save_matrix_to_csv(matrix, n, filename) !подпрограмма, которая вызывается позже для записи полученной матрицы U в файл
                real, dimension(:,:), intent(in) :: matrix
                character(len=*), intent(in) :: filename
                integer :: i, j, n
        end subroutine save_matrix_to_csv
            
    end interface 
    
    
    call MPI_INIT(ierr) !инициализация MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierr) !содержит ранг текущего процесса
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size_cluster, ierr) !содержит общее кол-во процессов
    local_n = n / size_cluster !получение кол-ва элементов, обрабатываемых каждым процессом

    if (process_Rank == 0) then   !на процессе с рангом 0 инициализируем значений

        print *, "U(x,t) = 5*exp(x)*cos(t)"
        print *, "f(x,t) = -5*exp(x)*(sin(t) + cos(t))"
        print *, "U(x,0) = 5*exp(x)"
        print *, "U(0,t) = 5*cos(t)"
        print *, "U(1,t) = 5*exp(1)*cos(t)"
        print *, "tau = ", tau
        print *, "h = ", h


        do i = 1, n !начальное приближение для U и вычисление эталонного решения
            x_curr = h*(i-1)
            do j = 1, n
            	t_curr = tau*(j-1)
                U_et(i,j) = 5*exp(x_curr)*cos(t_curr)
                f(i, j) = -5*exp(x_curr)*(sin(t_curr) + cos(t_curr))
            end do
            t_curr = tau*(i-1)
            Ux(i, 1) = 5*exp(x_curr) ! Задаем НУ, для t=0
            Ux(1, i) = 5*cos(t_curr)! Задаем ГУ, для x=0
            Ux(n, i) = 5*exp(1.0)*cos(t_curr) ! Задаем ГУ, для x=1 (последней точки стержня)
        end do
	A = 0
        A(1,1) = 1
        A(n,n) = 1
        do i = 2, n-1
            A(i, i-1) = -tau/(h**2)
            A(i, i) = 1 + 2*tau/(h**2)
            A(i, i+1) = -tau/(h**2)
        end do
        
       
    end if

    call cpu_time(start)
    do j = 2, n
      Uxt = 1
      Bt = Ux(:, j - 1) + tau * f(:, j - 1) !заполняем вектор В на основе данных с прошлого слоя
      call solve(A, Bt, Uxt, n,local_n, eps) !решается СЛАУ 
      Ux(2:n-1, j) = Uxt(2:n-1) !решение сохраняется в массив
    end do
    call cpu_time(finish)

    if (process_Rank == 0) then   !выводим значения эталонного решения и обычного, и время выполнения программы
        print *, "time to find u:", finish - start
        print *, "U etalon"
        do i = 1, 10
            print *, U_et(i, 90)
        end do

        print *, "U iter"
        do i = 1, 10
            print *, Ux(i, 90)
        end do

        call save_matrix_to_csv(U_et(:, :), n - 1, "U_etalon.csv")
        call save_matrix_to_csv(Ux(:, :), n - 1, "U_x.csv")
    endif
    call MPI_FINALIZE(ierr) !завершение работы MPI, освобождение ресурсов


end program main

subroutine solve(A, b, x0, n,local_n, eps)
    use mpi_sol
    integer :: n, ierror
    real :: A(n, n), b(n), x1(n), D(n,n), C(n,n), F(n,n), eps, r(n)
    real, intent(inout) :: x0(n)

    
    call find_d_c_f_r(A, b, n, D, C, F, r)
    call find_x(x0, x1, n,local_n, eps, F, r)

end subroutine solve

subroutine save_matrix_to_csv(matrix, n, filename)
    real, dimension(:,:), intent(in) :: matrix
    character(len=*), intent(in) :: filename
  
    integer :: i, j, n

    open(10, file=filename, status='replace', action='write') ! Открываем файл для записи

    do i = 1, n
        do j = 1, n
            write(10, '(A, ES14.7, A)', advance='no') ' ', matrix(i,j), ' '
        end do
        write(10, *) ! Переходим на новую строку после каждой строки матрицы
    end do

    close(10) ! Закрываем файл
    print*, 'Матрица успешно сохранена в файл ', trim(filename)
end subroutine save_matrix_to_csv


