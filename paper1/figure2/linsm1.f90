program LinearSeismic
!----compile this code using "pgfortran linsm1.f90 support.o -o lsm6 -O3 -r8"
!----could also use gfortran, but some modifications required
use support	!---support module is compiled with "pgfortran -c support.f90 -r8"
implicit none
real, parameter::pi=3.14159265358979323846

Integer,parameter::nt=5,nk=201,nx=560,ny=nx
real,Dimension(nt)::t
real,Dimension(nx)::x
real,Dimension(ny)::y
real,Dimension(nk)::k,integrandA2,integrandB2,integrandAB1,integrandAB3,bj1kr,kbj0kr
real,dimension(nx,ny,nt)::A2,B2,AB1,AB3,Ux,Uy,Umag
real,dimension(nk,nt)::A2kt,B2kt,AB1kt,AB3kt
real,dimension(nk,nx,ny)::A2b1kxy
real,dimension(nk,nx,ny,nt)::intA2kxyt
real,parameter::vp1=1.5,vs1=0.5,mu=vs1**2,lambda=vp1**2-2*vs1**2		!----solution is for scaled rho=1
real,parameter::a=0.1,F0=1.,G0=0.,alpha=sqrt(mu),beta=sqrt((lambda+2*mu)),maxK=90.
integer::ik,ix,iy,it,i
character(16)::filename
real::r,dk

!-----timing  variables
integer hours, minutes,seconds,milliseconds,day,elapsedHours,elapsedMinutes,elapsedSeconds,elapsedMs,dw
integer, dimension(8) :: time_array

!-----Set the start time variables-----!
CALL date_and_time(VALUES=time_array)
day=time_array(3)
hours=time_array(5)
minutes=time_array(6)
seconds=time_array(7)
milliseconds=time_array(8)
!-----Print the start time to screen-----!
print *, 'started: day', day,': time', hours,':', minutes,':', seconds

call linspace(x, -pi-1.,pi-1., nx)	!-----offset the initial pulse from (0,0) to (1,1) by moving the coordinate system
call linspace(y, -pi-1.,pi-1., ny)
call linspace(t,0.,4., nt)
call linspace(k,0.,maxK,nk)
dk=maxK/(nk-1.)

call savevector('xlin.out',len('xlin.out'),x,nx)
call savevector('ylin.out',len('ylin.out'),y,ny)
call savevector('tlin.out',len('tlin.out'),t,nt)

do ik=1,nk
	do it=1,nt
		A2kt(ik,it)=exp(-a**2*k(ik)**2/4.)*cos(k(ik)*alpha*t(it))*k(ik)**2
		B2kt(ik,it)=exp(-a**2*k(ik)**2/4.)*cos(k(ik)*beta*t(it))*k(ik)**2
		AB1kt(ik,it)=exp(-a**2*k(ik)**2/4.)*(cos(k(ik)*alpha*t(it))-cos(k(ik)*beta*t(it)))
	enddo
enddo
do ik=1,nk
	do ix=1,nx
		do iy=1,ny
			A2b1kxy(ik,ix,iy)=bessel_j1(k(ik)*sqrt(x(ix)**2+y(iy)**2))
		enddo
	enddo
enddo
do ix=1,nx
    do iy=1,ny
		r=sqrt(x(ix)**2+y(iy)**2)
		bj1kr=bessel_j1(k*r)
		kbj0kr=k*bessel_j0(k*r)
        do it=1,nt
		    integrandA2=A2kt(:,it)*bj1kr
            integrandB2=B2kt(:,it)*bj1kr
		    integrandAB1=AB1kt(:,it)*bj1kr
            integrandAB3=AB1kt(:,it)*kbj0kr
            
            call trapz(dk,integrandB2,nk,B2(ix,iy,it))
		    call trapz(dk,integrandA2,nk,A2(ix,iy,it))
		    call trapz(dk,integrandAB1,nk,AB1(ix,iy,it))
		    call trapz(dk,integrandAB3,nk,AB3(ix,iy,it))
		    
            Ux(ix,iy,it)=-0.5*a**2*x(ix)/r**3*(x(ix)**2*F0+y(iy)**2*G0)*B2(ix,iy,it)&
			-0.5*a**2*(F0-G0)*(2*x(ix)*(x(ix)**2-3*y(iy)**2)/r**5*AB1(ix,iy,it)&
                +x(ix)*y(iy)**2/r**3*A2(ix,iy,it)&
                -x(ix)*(x(ix)**2-3*y(iy)**2)/r**4*AB3(ix,iy,it) )
        
            Uy(ix,iy,it)=-0.5*a**2*y(iy)/r**3*(x(ix)**2*F0+y(iy)**2*G0)*B2(ix,iy,it)&
                +0.5*a**2*(F0-G0)*(-2*y(iy)*(3*x(ix)**2-y(iy)**2)/r**5*AB1(ix,iy,it)&
                +x(ix)**2*y(iy)/r**3*A2(ix,iy,it)&
                +y(iy)*(3*x(ix)**2-y(iy)**2)/r**4*AB3(ix,iy,it) )
        
            Umag(ix,iy,it)=sqrt(Ux(ix,iy,it)**2+Uy(ix,iy,it)**2)!*sign(1.,Ux(ix,iy,it)*x(ix))
        enddo
    enddo
enddo
do it=1,nt
	filename='data/umag'//char(48+it)//'.out'
	call savearray(filename,len(filename),Umag(:,:,it),nx,ny)
enddo

call printelapsedtime(day,hours,minutes,seconds,milliseconds)
end
subroutine trapz(dx, y, ny, integral)
    ! Calculates the integral of a function y with respect to x using the trapezoidal rule.
	implicit none
	integer,intent(in)	::ny
    real, dimension(ny),intent(in)  ::y  !---function y(x)
	real, intent(in)  	:: dx
    real, intent(out)	:: integral      !--- result
	integer::i
    ! Integrate using the trapezoidal rule
	integral=(y(1)+y(ny))*dx*0.5	!---half weighting for endpoints
	do i=2,ny-1
		integral=integral+y(i)*dx
	enddo
end subroutine