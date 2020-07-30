! start the module containing the matmul kernel
module mmul_mod
 use cudafor
contains
 ! mmul_kernel computes A*B into C where
 ! A is NxM, B is MxL, C is then NxL
 attributes(global) subroutine mmul_kernel( A, B, C, N, M, L )
  real :: A(N,M), B(M,L), C(N,L)
  integer, value :: N, M, L
  integer :: i, j, kb, k, tx, ty
  ! submatrices stored in shared memory
  real, shared :: Asub(16,16), Bsub(16,16)
  ! the value of C(i,j) being computed
  real :: Cij
  ! Get the thread indices
  tx = threadidx%x
  ty = threadidx%y

  ! This thread computes C(i,j) = sum(A(i,:) * B(:,j))
  i = (blockidx%x-1) * 16 + tx
  j = (blockidx%y-1) * 16 + ty
  Cij = 0.0
  ! Do the k loop in chunks of 16, the block size
  do kb = 1, M, 16
    ! Fill the submatrices
    ! Each of the 16x16 threads in the thread block
    ! loads one element of Asub and Bsub
    Asub(tx,ty) = A(i,kb+ty-1)
    Bsub(tx,ty) = B(kb+tx-1,j)
    ! Wait until all elements are filled
    call syncthreads()
    ! Multiply the two submatrices
    ! Each of the 16x16 threads accumulates the
    ! dot product for its element of C(i,j)
    do k = 1,16
		Cij = Cij + Asub(tx,k) * Bsub(k,ty)
    enddo
    ! Synchronize to make sure all threads are done
    ! reading the submatrices before overwriting them
    ! in the next iteration of the kb loop
    call syncthreads()
  enddo
  ! Each of the 16x16 threads stores its element
  ! to the global C array
  C(i,j) = Cij
  end subroutine mmul_kernel

! The host routine to drive the matrix multiplication
 subroutine mmul( A, B, C )
  real, dimension(:,:) :: A, B, C
  ! allocatable device arrays
  real, device, allocatable, dimension(:,:) :: Adev,Bdev,Cdev
  ! dim3 variables to define the grid and block shapes
  type(dim3) :: dimGrid, dimBlock

  ! Get the array sizes
  N = size( A, 1 )
  M = size( A, 2 )
  L = size( B, 2 )
  ! Allocate the device arrays
  allocate( Adev(N,M), Bdev(M,L), Cdev(N,L) )

  ! Copy A and B to the device
  Adev = A(1:N,1:M)
  Bdev(:,:) = B(1:M,1:L)

  ! Create the grid and block dimensions
  dimGrid = dim3( N/16, L/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Cdev, N, M, L)

  ! Copy the results back and free up memory
  C(1:N,1:L) = Cdev
  deallocate( Adev, Bdev, Cdev )
 end subroutine mmul
 
 ! The host routine to drive the double matrix multiplication
 subroutine mmul2( A, B, C, D)
  real, dimension(:,:) :: A, B, C, D
  ! allocatable device arrays
  real, device, allocatable, dimension(:,:) :: Adev,Bdev,Cdev,Ddev,Tempdev
  ! dim3 variables to define the grid and block shapes
  type(dim3) :: dimGrid, dimBlock

  ! Get the array sizes
  i = size( A, 1 )
  m = size( A, 2 )
  n = size( B, 2 )
  j = size( C, 2 )
  ! Allocate the device arrays
  allocate( Adev(i,m), Bdev(m,n), Tempdev(i,n), Cdev(n,j), Ddev(i,j) )

  ! Copy A and B to the device
  Adev = A(1:i,1:m)
  Bdev(:,:) = B(1:m,1:n)
  Cdev(:,:) = C(1:n,1:j)

  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D(1:i,1:j) = Ddev
  deallocate( Adev, Bdev, Cdev, Ddev, Tempdev)
 end subroutine mmul2
  
  ! The host routine to drive the double matrix multiplication twice
 subroutine mmul22( A, B, B2, C, D, D2)
  real, dimension(:,:) :: A, B, B2, C, D, D2
  ! allocatable device arrays
  real, device, allocatable, dimension(:,:) :: Adev,Bdev,Cdev,Ddev,Tempdev
  ! dim3 variables to define the grid and block shapes
  type(dim3) :: dimGrid, dimBlock

  ! Get the array sizes
  i = size( A, 1 )
  m = size( A, 2 )
  n = size( B, 2 )
  j = size( C, 2 )
  ! Allocate the device arrays
  allocate( Adev(i,m), Bdev(m,n), Tempdev(i,n), Cdev(n,j), Ddev(i,j))

  ! Copy A and B to the device
  Adev = A(1:i,1:m)
  Bdev(:,:) = B(1:m,1:n)
  Cdev(:,:) = C(1:n,1:j)
!first transform
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)
  ! Copy the results back
  D(1:i,1:j) = Ddev
!second transform

  Bdev(:,:) = B2(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D2(1:i,1:j) = Ddev
  deallocate( Adev, Bdev, Cdev, Ddev, Tempdev)
 end subroutine mmul22
 
   ! The host routine to drive the double matrix multiplication four times
 subroutine mmul24( A, B, B2, B3, B4, C, D, D2, D3, D4)
  real, dimension(:,:) :: A, B, B2, B3, B4, C, D, D2, D3, D4
  ! allocatable device arrays
  real, device, allocatable, dimension(:,:) :: Adev,Bdev,Cdev,Ddev,Tempdev
  ! dim3 variables to define the grid and block shapes
  type(dim3) :: dimGrid, dimBlock

  ! Get the array sizes
  i = size( A, 1 )
  m = size( A, 2 )
  n = size( B, 2 )
  j = size( C, 2 )
  ! Allocate the device arrays
  allocate( Adev(i,m), Bdev(m,n), Tempdev(i,n), Cdev(n,j), Ddev(i,j))

  ! Copy A and B to the device
  Adev = A(1:i,1:m)
  Bdev(:,:) = B(1:m,1:n)
  Cdev(:,:) = C(1:n,1:j)
!first transform
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)
  ! Copy the results back
  D(1:i,1:j) = Ddev
!second transform

  Bdev(:,:) = B2(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D2(1:i,1:j) = Ddev
  
  !third transform

  Bdev(:,:) = B3(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D3(1:i,1:j) = Ddev
  
  !fourth transform

  Bdev(:,:) = B4(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D4(1:i,1:j) = Ddev
  
  
  deallocate( Adev, Bdev, Cdev, Ddev, Tempdev)
 end subroutine mmul24
 
 subroutine mmul26( A, B, B2, B3, B4, B5, B6, C, D, D2, D3, D4, D5, D6)
  real, dimension(:,:) ::  A, B, B2, B3, B4, B5, B6, C, D, D2, D3, D4, D5, D6
  ! allocatable device arrays
  real, device, allocatable, dimension(:,:) :: Adev,Bdev,Cdev,Ddev,Tempdev
  ! dim3 variables to define the grid and block shapes
  type(dim3) :: dimGrid, dimBlock

  ! Get the array sizes
  i = size( A, 1 )
  m = size( A, 2 )
  n = size( B, 2 )
  j = size( C, 2 )
  ! Allocate the device arrays
  allocate( Adev(i,m), Bdev(m,n), Tempdev(i,n), Cdev(n,j), Ddev(i,j))

  ! Copy A and B to the device
  Adev = A(1:i,1:m)
  Bdev(:,:) = B(1:m,1:n)
  Cdev(:,:) = C(1:n,1:j)
!first transform
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)
  ! Copy the results back
  D(1:i,1:j) = Ddev
!second transform

  Bdev(:,:) = B2(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D2(1:i,1:j) = Ddev
  
  !third transform

  Bdev(:,:) = B3(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D3(1:i,1:j) = Ddev
  
  !fourth transform

  Bdev(:,:) = B4(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D4(1:i,1:j) = Ddev
  
  !fifth transform

  Bdev(:,:) = B5(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D5(1:i,1:j) = Ddev
  
  !sixth transform

  Bdev(:,:) = B6(1:m,1:n)
  ! Create the grid and block dimensions - first matmul
  dimGrid = dim3( i/16, n/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Adev, Bdev, Tempdev, i, m, n)
  
  ! Create the grid and block dimensions - second matmul
  dimGrid = dim3( i/16, j/16, 1 )
  dimBlock = dim3( 16, 16, 1 )
  call mmul_kernel<<<dimGrid,dimBlock>>>( Tempdev, Cdev, Ddev, i, n, j)

  ! Copy the results back and free up memory
  D6(1:i,1:j) = Ddev
  
  
  deallocate( Adev, Bdev, Cdev, Ddev, Tempdev)
 end subroutine mmul26
 
 
 end module mmul_mod

