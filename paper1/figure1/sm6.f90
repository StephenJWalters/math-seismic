program SeismicIntegrator
!This routine integrates the systems of ODEs for the 2-d elastic equations
!
!pgfortran sm4.f90 support.o -o sm4 -O3 -acc -ta=nvidia:cuda9.2 -r8

!granite - Vp~5, Vs~3, rho~2600
!clay - Vp~1.5, Vs~0.5, rho~2200

use support	! pgfortran -c support.f90 -r8
implicit none
real, parameter::pi=3.14159265358979323846

Integer,parameter::nimages=4				!-----number of times to output the density for an image---!
Real,parameter::t=4./nimages,h=0.005		!-----time units per image-----time for initial pass through the RK4 routine-----!
real,parameter,Dimension(4)::deltvec=(/0.,h/2,h/2,h/),rk4wt=(/h/6.,h/3.,h/3.,h/6./)	!rk4 spacing and weighting

Integer,parameter::loops=NInt(t/h),totalsteps=loops*nimages

integer, parameter::nxx=7*80,nyy=nxx,dampwidth=56
integer, parameter::Mmax=nxx/5-1,Nmax=nyy/5-1

Integer,parameter::numsensors=nxx/10+1,sensorj=nyy-dampwidth-10		!-----location of sensor-----!
integer,Dimension(numsensors)::sensorsi
real,Dimension(numsensors)::sensorsx
integer,parameter::ds=(nxx-2*dampwidth)/(numsensors-1)

real,Dimension(totalsteps,numsensors)::uxout,uyout

real,Dimension(nxx)::x
real,Dimension(nyy)::y
real,Dimension(nxx,nyy)::u,ux,uy,damping,rhoinv,mu,lambda,trapwt
real,Dimension(0:Mmax,nxx)::Cx,Sx,gmCx
real,Dimension(0:Nmax,nyy)::Cy,Sy,gnCy
real,Dimension(0:Mmax)::gm
real,Dimension(0:Nmax)::gn

!arrays for the Fourier coefficients.
real,Dimension(0:Mmax,0:Nmax)::AAmn,Atemp,AI,Aprime,BBmn,Btemp,BI,Bprime,CCmn,Ctemp,CI,Cprime,DDmn,Dtemp,DI,Dprime

real::argmpx,argnpy,wt,sumux,sumuy,looptime,dampdepth,r	!---time through the RKF routine for each image---!
real,parameter::a=0.1,F0=1.,G0=1.,df=0.9,df2=19,x0=1.,y0=1.
real,parameter::vp1=1.5,vs1=0.5,rho1=2200.,vp2=5.,vs2=3.,rho2=2600./rho1,mu1=vs1**2,mu2=vs2**2,lambda1=vp1**2-2*vs1**2,lambda2=vp2**2-2*vs2**2
integer::m,n,i,j,rk4step,image,m1,n1,iLoop,stepnum=0
character(16)::filename
real::dx,dy

!-----timing  variables
integer hours, minutes,seconds,milliseconds,day,elapsedHours,elapsedMinutes,elapsedSeconds,elapsedMs,dw
integer, dimension(8) :: time_array

!-----initial velocity field
CCmn=0.
DDmn=0.
!-----Set the start time variables-----!
CALL date_and_time(VALUES=time_array)
day=time_array(3)
hours=time_array(5)
minutes=time_array(6)
seconds=time_array(7)
milliseconds=time_array(8)
!-----Print the start time to screen-----!
print *, 'started: day', day,': time', hours,':', minutes,':', seconds
print*, 'num images: ',nimages
print*,'loops per image: ',loops
!print*,'ds: ',ds
call linspace(x, -pi,pi, nxx)
dx=2*pi/(nxx-1)
!print*,x(1),x(nxx)
call linspace(y, -pi,pi, nyy)
dy=2*pi/(nyy-1)
call savevector('x.out',len('x.out'),x,nxx)
call savevector('y.out',len('y.out'),y,nyy)

call saveinteger('dampwidth.out',len('dampwidth.out'),dampwidth)
call saveinteger('nimages.out',len('nimages.out'),nimages)

call savereal('t.out',len('t.out'),t)
call saveinteger('totalsteps.out',len('totalsteps.out'),totalsteps)
call saveinteger('numsensors.out',len('numsensors.out'),numsensors)

call savevector('sensorx.out',len('sensorx.out'),x,nxx)
call savereal('sensory.out',len('sensory.out'),y(sensorj))


trapwt=1.
do i=1,nxx
	trapwt(i,1)=0.5
	trapwt(i,nyy)=0.5
enddo
do j=1,nyy
	trapwt(1,j)=0.5
	trapwt(nxx,j)=0.5
enddo
trapwt(1,1)=0.25
trapwt(nxx,nyy)=0.25
trapwt(1,nyy)=0.25
trapwt(nxx,1)=0.25
dw=real(dampwidth)
open(unit=2,file='details.txt')
	write(2,*, ADVANCE = "YES") 'F0 = ',F0
	write(2,*, ADVANCE = "YES") 'G0 = ',G0
	write(2,*, ADVANCE = "YES") 'IC a = ',a
	write(2,*, ADVANCE = "YES") 'damping df = ',df
	write(2,*, ADVANCE = "YES") 'damping width dw = ',dw
	write(2,*, ADVANCE = "YES") 'init pulse x0 = ',x0
	write(2,*, ADVANCE = "YES") 'init pulse y0 = ',y0
	write(2,*, ADVANCE = "YES") 'Mmodes = ',Mmax
	write(2,*, ADVANCE = "YES") 'Nmodes = ',Nmax
	write(2,*, ADVANCE = "YES") 'x points nxx= ',nxx
	write(2,*, ADVANCE = "YES") 'y points nyy= ',nyy	
	write(2,*, ADVANCE = "YES") 'num images= ',nimages	
	write(2,*, ADVANCE = "YES") 'time per image= ',t	
	write(2,*, ADVANCE = "YES") 'timestep= ',h	
close(2)

!locate the sensors
do i=1,numsensors
	
	sensorsi(i)=dampwidth+ds*(i-1)
	sensorsx(i)=x(sensorsi(i))
enddo
call savevector('sensorsx.out',len('sensorsx.out'),sensorsx,numsensors)
!create the damping fields for absorbing boundaries
damping=0.

do i=1,dampwidth
	dampdepth=1.-real(i)/dw
	do j=1,nyy
		damping(i,j)=damping(i,j)+exp(df*dampdepth)-1.
	enddo
enddo
do i=nxx-dampwidth+1,nxx
	dampdepth=1.-real(nxx-i+1)/dw
	do j=1,nyy
		damping(i,j)=damping(i,j)+exp(df*dampdepth)-1.
	enddo
enddo
do j=1,dampwidth
	dampdepth=1.-real(j)/dw
	do i=1,nxx
		damping(i,j)=damping(i,j)+exp(df*dampdepth)-1.
	enddo
enddo
do j=nyy-dampwidth+1,nyy
	dampdepth=1.-real(nyy-j+1)/dw
	do i=1,nxx
		damping(i,j)=damping(i,j)+exp(df*dampdepth)-1.
	enddo
enddo
damping=df2*damping
!Create and store the trig functions  cos(m/2(pi+x))  and  sin(m/2(pi+x)) .
do m=0,Mmax
    do i=1,nxx
        argmpx=(x(i)+pi)*m/2.
        Cx(m,i)=cos(argmpx)
        Sx(m,i)=sin(argmpx)
    end do
end do

do n=0,Nmax
    do j=1,nyy
        argnpy=(y(j)+pi)*n/2.
        Cy(n,j)=cos(argnpy)
        Sy(n,j)=sin(argnpy)
    end do
end do

!$acc kernels
do i=1,nxx
	do j=1,nyy
		ux(i,j)=F0*2.*(x(i)-x0)/a**2*exp(-((x(i)-x0)**2+(y(j)-y0)**2)/a**2)
		uy(i,j)=G0*2.*(y(j)-y0)/a**2*exp(-((x(i)-x0)**2+(y(j)-y0)**2)/a**2)
   end do
end do
!$acc end kernels
u=sqrt(ux**2+uy**2)
call savearray('data/disp00.out',len('data/disp00.out'),u,nxx,nyy)

!-----Integration weightings
gm=1.
gm(0)=0.5
gn=1.
gn(0)=0.5

AAmn=0.
BBmn=0.
!$acc kernels
do m=0,Mmax
    do n=0,Nmax!
	    do i=1,nxx
            do j=1,nyy
				wt=gm(m)*gn(n)*dx*dy*Cx(m,i)*Cy(n,j)/pi**2
                AAmn(m,n)=AAmn(m,n)+wt*ux(i,j)*trapwt(i,j)
                BBmn(m,n)=BBmn(m,n)+wt*uy(i,j)*trapwt(i,j)
            end do
        end do
    end do
end do
!$acc end kernels

ux=matmul(matmul(transpose(cx),aamn),cy)
uy=matmul(matmul(transpose(cx),bbmn),cy)

u=sqrt(ux**2+uy**2)
call savearray('data/disp0.out',len('data/disp0.out'),u,nxx,nyy)

!set up the density distribution
rhoinv=1.
mu=mu1
lambda=lambda1

do m=0,Mmax
	do i=1,nxx
		gmcx(m,i)=gm(m)*Cx(m,i)
	enddo
enddo

do n=0,Nmax
	do j=1,nyy
		gncy(n,j)=gn(n)*Cy(n,j)
	enddo
enddo

!Do the integration process.
do image=1,nimages
	print*,'interval ',image,' started'
	looptime=0.
	do iLoop=1,loops
		Atemp=AAmn
		Btemp=BBmn
		Ctemp=CCmn
		Dtemp=DDmn
	
		do rk4step=1,4			!call the 4 steps of the rk4 routine
			AI=Atemp+deltvec(rk4step)*Aprime
			BI=Btemp+deltvec(rk4step)*Bprime
			CI=Ctemp+deltvec(rk4step)*Cprime
			DI=Dtemp+deltvec(rk4step)*Dprime
						
			call SeismicRHS(AI,BI,CI,DI,Cprime,Dprime,Cx,Sx,Cy,Sy,GmCx,GnCy,nxx,nyy,trapwt,Mmax,Nmax,damping,rhoinv,mu,lambda,gm,gn,dx,dy,pi)
			
			Aprime=CI
			Bprime=DI
			
			AAmn=AAmn+rk4wt(rk4step)*Aprime
			BBmn=BBmn+rk4wt(rk4step)*Bprime
			CCmn=CCmn+rk4wt(rk4step)*Cprime
			DDmn=DDmn+rk4wt(rk4step)*Dprime
			
		end do
		looptime=looptime+h
		stepnum=stepnum+1
		sumux=0.
		sumuy=0.
		
		!!$acc kernels
		do n1=1,numsensors
			sumux=0.
			sumuy=0.
			m1=sensorsi(n1)
			do n=0,Nmax
				do m=0,Mmax
					sumux=sumux+AAmn(m,n)*Cx(m,m1)*Cy(n,sensorj)
					sumuy=sumuy+BBmn(m,n)*Cx(m,m1)*Cy(n,sensorj)
				enddo
			end do
			uxout(stepnum,n1)=sumux
			uyout(stepnum,n1)=sumuy
		end do
		!!$acc end kernels
	end do
	print*,'looptime: ',looptime

!Compute the Displacement function and save it.
ux=matmul(matmul(transpose(cx),aamn),cy)
uy=matmul(matmul(transpose(cx),bbmn),cy)
u=sqrt(ux**2+uy**2)
print*,'stepnum: ',stepnum
print *,'max displacement=  ',maxval(u)

	if (image<10) then
		filename='data/disp'//char(48+image)//'.out'		!--------disp0 to disp9------!
	elseif (image<100) then
		filename='data/disp'//char(48+image/10)//char(48+image-10*(image/10))//'.out'		!--------disp10 to disp99----!
	elseif (image<1000) then
		filename='data/disp'//char(48+image/100)//char(48+image/10-10*(image/100))//char(48+image-10*(image/10))//'.out'	!--------disp100 to disp999----!
	endif
call savearray(filename,len(filename),u,nxx,nyy)	
call savearray('data/ux.out',len('data/ux.out'),uxout,totalsteps,numsensors)
call savearray('data/uy.out',len('data/uy.out'),uyout,totalsteps,numsensors)
call printelapsedtime(day,hours,minutes,seconds,milliseconds)
	
end do
call printenddatetime()


end program SeismicIntegrator

Subroutine SeismicRHS(AI,BI,CI,DI,Cprime,Dprime,Cx,Sx,Cy,Sy,GmCx,GnCy,nxx,nyy,trapwt,Mmax,Nmax,damping,rhoinv,mu,lambda,gm,gn,dx,dy,pi)
!This routine evaluates the right-hand sides of the four sets of
!ordinary differential equations for the Fourier coefficients
use mmul_mod
implicit none
Integer,intent(in)::nxx,nyy,Mmax,Nmax
real,Dimension(nxx,nyy),intent(in)::damping,rhoinv,mu,lambda,trapwt
!trig basis functions
real,Dimension(0:Mmax,nxx),intent(in)::Cx,Sx,GmCx
real,Dimension(0:Nmax,nyy),intent(in)::Cy,Sy,GnCy
!Arrays for the Fourier coefficients.
real,Dimension(0:Mmax,0:Nmax),intent(in)::AI,BI,CI,DI
real,Dimension(0:Mmax,0:Nmax),intent(out)::Cprime,Dprime
real,Dimension(0:Mmax)::gm
real,Dimension(0:Nmax)::gn
real,Dimension(nxx,0:Nmax)::C
real,Dimension(0:Mmax,0:Nmax)::mmA25,mnA25,nnA25,mmB25,mnB25,nnB25
!dummy spatial variables and counters
real,Dimension(nxx,nyy)::uxxx,uxxy,uxyy,uyyy,uyxy,uyxx,ux,uy,vx,vy
integer::i,j,m,n
real,intent(in)::dx,dy,pi
real::wt

!  calculate the spatial variables from the fourier coefficients
!$acc kernels
do m=0,Mmax
    do n=0,Nmax
		mmA25(m,n)=m*m*AI(m,n)*0.25
		mnA25(m,n)=m*n*AI(m,n)*0.25
		nnA25(m,n)=n*n*AI(m,n)*0.25
		mmB25(m,n)=m*m*BI(m,n)*0.25
		mnB25(m,n)=m*n*BI(m,n)*0.25
		nnB25(m,n)=n*n*BI(m,n)*0.25
	end do
end do
!$acc end kernels
call mmul( transpose(cx), -mmA25, C )
call mmul(C,cy,uxxx)
call mmul( transpose(sx), mnA25, C )
call mmul(C,sy,uxxy)
call mmul( transpose(cx), -nnA25, C )
call mmul(C,cy,uxyy)
call mmul( transpose(cx), -mmB25, C )
call mmul(C,cy,uyxx)
call mmul( transpose(sx), mnB25, C )
call mmul(C,sy,uyxy)
call mmul( transpose(cx), -nnB25, C )
call mmul(C,cy,uyyy)

vx=matmul(matmul(transpose(cx),CI),cy)
vy=matmul(matmul(transpose(cx),DI),cy)

! evaluate Cprime and Dprime (i.e. u'') from the spatial variables.
ux=rhoinv*((lambda+2.*mu)*uxxx+(lambda+mu)*uyxy+mu*uxyy)-vx*damping
uy=rhoinv*((lambda+2.*mu)*uyyy+(lambda+mu)*uxxy+mu*uyxx)-vy*damping

call mmul(  trapwt*ux*dx*dy/pi**2,transpose(gncy), C )
call mmul(gmcx,C,Cprime)
call mmul(  trapwt*uy*dx*dy/pi**2,transpose(gncy), C )
call mmul(gmcx,C,Dprime)

end Subroutine SeismicRHS
