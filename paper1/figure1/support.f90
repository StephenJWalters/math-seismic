module support
! pgfortran -c support.f90 -r8
! pgfortran -o sm3 sm3.f90 support.o  -O3 -acc -ta=nvidia:cuda9.2 -r8
contains
  SUBROUTINE lgwt(n1,xa,xb,x,w,pi)

! lgwt.m
!
! This script is for computing definite integrals using Legendre-Gauss 
! Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
! [a,b] with truncation order N
!
! Suppose you have a continuous function f(x) which is defined on [a,b]
! which you can evaluate at any x in [a,b]. Simply evaluate it at all of
! the values contained in the x vector to obtain a vector f. Then compute
! the definite integral using sum(f.*w)
!
! Written by Greg von Winckel - 02/25/2004
implicit none
integer, intent(in)::n1
real, intent(in)::xa,xb,pi
real,dimension(n1), intent(out)::x,w

integer n,n2,i,k,p(n1)
real eps,LG(n1,n1+1),inc
real,dimension(n1)::xu,y,y0,LP

N=N1-1
N2=N+2
do i=0,n
	p(i+1)=i
end do

!xu=linspace(-1,1,N1)'
inc=2./n
xu=0.
do i=0,n
	xu(i+1)=-1+i*inc
end do
!print*, xu
! Initial guess
y=cos((2.*p+1.)*pi/(2.*N+2.))+(0.27/N1)*sin(N*pi*xu/N2)

! Legendre-Gauss Vandermonde Matrix
LG=0

! Derivative of LGVM
Lp=0

! Compute the zeros of the N+1 Legendre Polynomial
! using the recursion relation and the Newton-Raphson method
y0=2
eps=epsilon(y0)
! Iterate until new points are uniformly within epsilon of old points
do while (maxval(abs(y-y0))>eps)
    LG(:,1)=1
    LG(:,2)=y
    do k=2,N1
        LG(:,k+1)=( (2*k-1)*y*LG(:,k)-(k-1)*LG(:,k-1) )/k
		!print *,'l',L(:,k+1)
    end do
    Lp=N2*( LG(:,N1)-y*LG(:,N2) )/(1-y**2)   
    y0=y
    y=y0-LG(:,N2)/Lp
end do

! Linear map from[-1,1] to [a,b]
x=(xa*(1-y)+xb*(1+y))/2      

! Compute the weights
w=(xb-xa)/((1-y**2)*Lp**2)*(real(N2)/real(N1))**2
end subroutine lgwt

Subroutine printelapsedtime(day,hours,minutes,seconds,milliseconds)
integer, dimension(8) :: time_array
integer hours, minutes,seconds,milliseconds,day,elapsedHours,elapsedMinutes,elapsedSeconds,elapsedMs

	CALL date_and_time(VALUES=time_array)
	elapsedHours=24*(time_array(3)-day)+(time_array(5)-hours)
	elapsedMinutes=time_array(6)-minutes
	elapsedSeconds=time_array(7)-seconds
	elapsedMs=time_array(8)-milliseconds
	
	if (elapsedMs<0) then
		elapsedMs=elapsedMs+1000
		elapsedSeconds=elapsedSeconds-1
	end if
	if (elapsedSeconds<0) then
		elapsedSeconds=elapsedSeconds+60
		elapsedMinutes=elapsedMinutes-1
	end if
	if (elapsedMinutes<0) then
		elapsedMinutes=elapsedMinutes+60
		elapsedHours=elapsedHours-1
	end if
	
	print *, 'time:', elapsedHours,'hrs,',elapsedMinutes,'mins',elapsedSeconds,'sec',elapsedMs,'ms'	
end Subroutine printelapsedtime


subroutine inv(a,c,n)!============================================================ !
!Inverse matrix ! Method: Based LU decomposition for Ax=b ! Marouf Abderahmane - Strasbourg University 2016
!----------------------------------------------------------- 
	implicit none 
	integer,intent(in)::n 
	real,dimension(n,n),intent(in)::a
	real,intent(out)::c(n,n) 
	real,dimension(n,n)::L, U,tempa
	real,dimension(n)::b, d, x
	real coeff
	integer i, j, k
	L=0.0
	U=0.0
	b=0.0 
	tempa=a
	do k=1, n-1 
		do i=k+1,n 
			coeff=tempa(i,k)/tempa(k,k) 
			L(i,k) = coeff 
			do j=k+1,n 
				tempa(i,j) = tempa(i,j)-coeff*tempa(k,j) 
			end do 
		end do 
	end do 
	do i=1,n 
		L(i,i) = 1.0 
	end do
	do j=1,n
		do i=1,j 
			U(i,j) = tempa(i,j)
		end do 
	end do 
	do k=1,n 
		b(k)=1.0 
		d(1) = b(1)
		do i=2,n
			d(i)=b(i)
			do j=1,i-1
				d(i) = d(i) - L(i,j)*d(j)
			end do
		end do 
		x(n)=d(n)/U(n,n)
		do i = n-1,1,-1
			x(i) = d(i)
			do j=n,i+1,-1
				x(i)=x(i)-U(i,j)*x(j)
			enddo
			x(i) = x(i)/u(i,i)
		end do
		do i=1,n
			c(i,k) = x(i)
		end do
		b(k)=0.0
	end do
end subroutine inv

Subroutine leftInverse(A,Linv,m,n)
integer, intent(in)::m,n
real,dimension(m,n),intent(in)::A
real,dimension(n,m),intent(out)::Linv
real,dimension(n,m)::At
real,dimension(n,n) :: Ata,AtaInv

At=Transpose(A)
Ata=matmul(at,a)
call inv(Ata,atainv,N)
Linv=matmul(atainv,at)

end subroutine

Subroutine rightInverse(A,Rinv,m,n)
integer, intent(in)::m,n
real,dimension(m,n),intent(in)::A
real,dimension(n,m),intent(out)::Rinv
real,dimension(n,m)::At
real,dimension(m,m) :: Aat,AatInv

At=Transpose(A)
Aat=matmul(a,at)
call inv(Aat,aatinv,m)
Rinv=matmul(at,aatinv)

end subroutine


Subroutine fofalpha(k,alpha,pi,f,m)
implicit none
	real,intent(in)::k,alpha,pi
	integer,intent(in)::m
	real,intent(out)::f
	f=k*cos(m/2.*(alpha*pi+pi))-m*alpha/2.*sin(m/2.*(alpha*pi+pi))
	!f=1.
End Subroutine fofalpha

Subroutine bisection(left,right,root,pi,k,m)
implicit none
	real,intent(in)::left,right,pi,k
	real,intent(out)::root
	real::lt,rt,mid,fl,fr,fm
	integer,intent(in)::m
	
	call fofalpha(k,left,pi,fl,m)
	call fofalpha(k,right,pi,fr,m)
	if (fl*fr>0.) then
		print*,'brackets are wrong for m=: ', m,'l:',left,'r',right
		return
	endif
	
	mid=(left+right)*0.5
	call fofalpha(k,mid,pi,fm,m)
	lt=left
	rt=right
	do while(abs(fm)>0.000002)
		if (fm*fl<0.) then
			rt=mid
		else
			lt=mid
		end if
		mid=(lt+rt)*0.5
		call fofalpha(k,lt,pi,fl,m)
		call fofalpha(k,rt,pi,fr,m)
		call fofalpha(k,mid,pi,fm,m)
		!print*,fl,fm,fr,lt,mid,rt,fm*fl
	end do 
	root=mid
end subroutine bisection
Subroutine savearray(filename,filelen,array,imax,jmax)
implicit none
integer,intent(in)::imax,jmax,filelen
real,dimension(imax,jmax),intent(in)::array
character(filelen),intent(in)::filename
integer::i,j
open(unit=2,file=filename)
do i=1,imax
	do j=1,jmax
		write(2,*, ADVANCE = "yes") array(i,j)
	end do
end do
close(2)
end subroutine

Subroutine savevector(filename,filelen,vector,imax)
implicit none
integer,intent(in)::imax,filelen
real,dimension(imax),intent(in)::vector
character(filelen),intent(in)::filename
integer::i
open(unit=2,file=filename)
do i=1,imax
	write(2,*, ADVANCE = "yes") vector(i)
end do
close(2)
end subroutine

Subroutine saveinteger(filename,filelen,param)
implicit none
integer,intent(in)::filelen
integer,intent(in)::param
character(filelen),intent(in)::filename
open(unit=2,file=filename)
write(2,*) param
close(2)
end subroutine

Subroutine savereal(filename,filelen,param)
implicit none
integer,intent(in)::filelen
real,intent(in)::param
character(filelen),intent(in)::filename
open(unit=2,file=filename)
write(2,*) param
close(2)
end subroutine

Subroutine printenddatetime()
integer, dimension(8) :: time_array
integer hours, minutes,seconds,milliseconds,day,elapsedHours,elapsedMinutes,elapsedSeconds,elapsedMs
CALL date_and_time(VALUES=time_array)
day=time_array(3)
hours=time_array(5)
minutes=time_array(6)
seconds=time_array(7)
milliseconds=time_array(8)
!-----Print the end time to screen-----!
print *, 'ended: day', day,'time:', hours,':', minutes,':', seconds
end Subroutine printenddatetime

subroutine linspace(x, a, b, n)
implicit none
  integer, intent(in)   :: n
  real,dimension(n),intent(out) :: x
  real,intent(in)  :: a
  real,intent(in)  :: b
  real :: dx
  integer   :: i
	dx = (b-a) / (n-1.)
	do i=1,n
	  x(i) = (i-1)*dx+a
	enddo
end subroutine linspace 
end module