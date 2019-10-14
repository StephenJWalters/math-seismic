program SeismicIntegrator
!This routine integrates the systems of ODEs for the 2-d elastic equations
!
!pgfortran sm7.f90 support.o mmul_mod.o -o sm7 -O3 -acc -ta=nvidia:cuda9.2 -r8 -Mcuda

!granite - Vp~5, Vs~3, rho~2600
!clay - Vp~1.5, Vs~0.5, rho~2200
!use mmul_mod
use support	! pgfortran -c support.f90 -r8
implicit none
real, parameter::pi=3.14159265358979323846

Integer,parameter::nimages=8				!-----number of times to output the density for an image---!
Real,parameter::t=0.6/nimages,h=0.0005		!-----time units per image-----timestep for each pass through the RK4 routine-----!

real,parameter,Dimension(4)::odestep=(/0.,h/2,h/2,h/),odewt=(/h/6.,h/3.,h/3.,h/6./)	!rk4 spacing and weighting

Integer,parameter::loops=NInt(t/h),totalsteps=loops*nimages

integer, parameter::nxx=20*80,nyy=nxx/5,dampwidth=nxx/5,dampwidthy=nyy/5
integer, parameter::Mmax=nxx/5-1,Nmax=nyy/5-1

Integer,parameter::numsensors=nxx/20+1,sensorj=nyy-0*dampwidth-0		!-----location of sensor-----!
integer,Dimension(numsensors)::sensorsi
real,Dimension(numsensors)::sensorsx
integer,parameter::ds=(nxx-2*dampwidth)/(numsensors-1)

real,Dimension(totalsteps,numsensors)::vxout,vyout

real,Dimension(nxx)::x
real,Dimension(nyy)::y,vpy,vsy,rhoy
real,Dimension(nxx,nyy)::u,ux,uy,damping,rhoinv,mu,lambda,mux,muy,lambdax,lambday,trapwt
real,Dimension(0:Mmax,nxx)::Cx,Sx,CxinvT
real,Dimension(0:Nmax,nyy)::Cy,Sy
real,Dimension(nxx,0:Mmax)::Cxinv
real,Dimension(nyy,0:Nmax)::Cyinv
real,Dimension(0:Mmax,nxx)::gmCx
real,Dimension(0:Nmax,nyy)::gnCy

real,Dimension(0:Mmax)::gm
real,Dimension(0:Nmax)::gn

!arrays for the Fourier coefficients.
real,Dimension(0:Mmax,0:Nmax)::AAmn,Atemp,AI,Aprime,BBmn,Btemp,BI,Bprime,CCmn,Ctemp,CI,Cprime,DDmn,Dtemp,DI,Dprime,Lmn,Mmn

real::argmpx,argnpy,wt,sumvx,sumvy,looptime,dampdepth,r	!---time through the RKF routine for each image---!
real,parameter::a=0.05,F0=0.01,G0=0.01,df=0.9,df2=19.,x0=0.5,y0=-0.1,x1=-0.5,x2=2.5,y1=-0.4,y2=0.,v=(x2-x1)*(y2-y1)*0.25
real,parameter::alphax=pi/(x2-x1),alphay=pi/(y2-y1)
real,parameter::vp1=0.4,vs1=0.3,rho1=300.,vp2=3.5,vs2=1.8,rho2=900.
real::mu1,mu2,lambda1,lambda2
integer::m,n,i,j,step,image,m1,n1,iLoop,stepnum=0
character(16)::filename
real::dx,dy,time1=0.,time1a=0.,time2=0.,time2a=0.,time3=0.,time4=0.,time4a=0.,time4b=0.
real,dimension(4)::time=0.
!-----timing  variables
integer hours, minutes,seconds,milliseconds,day,elapsedHours,elapsedMinutes,elapsedSeconds,elapsedMs,dw
integer, dimension(8) :: time_array
print*,'nxx: ',nxx
!-----initial velocity field
do j=1,nyy
	
enddo

CCmn=0.
DDmn=0.
mu1=rho1*vs1**2
mu2=rho2*vs2**2
lambda1=rho1*(vp1**2-2.*vs1**2)
lambda2=rho2*(vp2**2-2.*vs2**2)
!-----Set the start time variables-----!
CALL date_and_time(VALUES=time_array)
day=time_array(3)
hours=time_array(5)
minutes=time_array(6)
seconds=time_array(7)
milliseconds=time_array(8)
!-----Print the start time to screen-----!
print *, 'started: day', day,': time', hours,':', minutes,':', seconds
print*,'loops: ',loops
print*,'ds: ',ds
call linspace(x, x1,x2, nxx)
dx=(x2-x1)/(nxx-1)
print*,x(1),x(nxx)
call linspace(y, y1,y2, nyy)
dy=(y2-y1)/(nyy-1)
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
trapwt=trapwt*dx*dy/v
call savevector('x.out',len('x.out'),x,nxx)
call savevector('y.out',len('y.out'),y,nyy)

call saveinteger('dampwidth.out',len('dampwidth.out'),dampwidth)
call saveinteger('dampwidthy.out',len('dampwidthy.out'),dampwidthy)
call saveinteger('nimages.out',len('nimages.out'),nimages)

call savereal('t.out',len('t.out'),t)
call saveinteger('totalsteps.out',len('totalsteps.out'),totalsteps)
call saveinteger('numsensors.out',len('numsensors.out'),numsensors)

call savevector('sensorx.out',len('sensorx.out'),x,nxx)
!call savereal('sensorx.out',len('sensorx.out'),x(sensori))
call savereal('sensory.out',len('sensory.out'),y(sensorj))
call savereal('x1.out',len('x1.out'),x1)
call savereal('x2.out',len('x2.out'),x2)
call savereal('y1.out',len('y1.out'),y1)
call savereal('y2.out',len('y2.out'),y2)

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
print*,'y1: ',y(1)
do j=1,dampwidthy
	dampdepth=1.-real(j)/dw
	do i=1,nxx
		damping(i,j)=damping(i,j)+exp(df*dampdepth)-1.
	enddo
enddo
do j=nyy-dampwidthy+1,nyy
	dampdepth=1.-real(nyy-j+1)/dw
	do i=1,nxx
		!damping(i,j)=damping(i,j)+exp(df*dampdepth)-1.
	enddo
enddo
print*,'dampdepthy: ',y(dampwidthy)
damping=df2*damping
!Create and store the trig functions  cos(m/2(pi+x))  and  sin(m/2(pi+x)) .
do m=0,Mmax
    do i=1,nxx
        argmpx=(x(i)-x1)*m*alphax
        Cx(m,i)=cos(argmpx)
        Sx(m,i)=sin(argmpx)
    end do
end do

do n=0,Nmax
    do j=1,nyy
        argnpy=(y(j)-y1)*n*alphay
        Cy(n,j)=cos(argnpy)
        Sy(n,j)=sin(argnpy)
    end do
end do
print*,'v: ',v
print*,'cx: ',maxval(abs(cx))
call rightinverse(Cx,Cxinv,Mmax+1,nxx)
call rightinverse(Cy,Cyinv,Nmax+1,nyy)
CxinvT=transpose(cxinv)
print*,'CxinvT: ',maxval(abs(CxinvT))

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

AAmn=0.
BBmn=0.
!$acc kernels
do m=0,Mmax
    do n=0,Nmax!
	    do i=1,nxx
            do j=1,nyy
				wt=gm(m)*gn(n)*Cx(m,i)*Cy(n,j)
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

!set up the density and wavespeed distribution

do j=1,nyy
	vpy(j)=vp2 - (vp2-vp1) *Exp((y(j)/0.05))
	vsy(j)=vs2 - (vs2-vs1) *Exp((y(j)/0.05))
	rhoy(j)=rho2 - (rho2-rho1) *Exp((y(j)/0.05))
enddo

do i=1,nxx
	do j=1,nyy
		mu(i,j)=vsy(j)**2*rhoy(j)
		lambda(i,j)=(vpy(j)**2-2.*vsy(j)**2)*rhoy(j)
		rhoinv(i,j)=1./rhoy(j)
   end do
end do

Lmn=matmul(matmul(CxinvT,lambda),cyinv)
Mmn=matmul(matmul(CxinvT,mu),cyinv)

!calculate and store mu, lambda and their first derivatives
mu=0.
lambda=0.
mux=0.
muy=0.
lambdax=0.
lambday=0.
!$acc kernels
do i=1,nxx	
    do j=1,nyy
        do m=0,Mmax
            do n=0,Nmax
                mu(i,j)=mu(i,j)+Mmn(m,n)*Cx(m,i)*Cy(n,j)
                lambda(i,j)=lambda(i,j)+Lmn(m,n)*Cx(m,i)*Cy(n,j)
                mux(i,j)=mux(i,j)-Mmn(m,n)*m*alphax*Sx(m,i)*Cy(n,j)
                muy(i,j)=muy(i,j)-Mmn(m,n)*n*alphay*Cx(m,i)*Sy(n,j)
                lambdax(i,j)=lambdax(i,j)-Lmn(m,n)*m*alphax*Sx(m,i)*Cy(n,j)
                lambday(i,j)=lambday(i,j)-Lmn(m,n)*n*alphay*Cx(m,i)*Sy(n,j)
            end do
        end do
    end do
end do
!$acc end kernels

mu=matmul(matmul(transpose(cx),Mmn),cy)
lambda=matmul(matmul(transpose(cx),Lmn),cy)

!Do the integration process.
do image=1,nimages
	print*,'interval ',image,' started'
	looptime=0.
	do iLoop=1,loops
		Atemp=AAmn
		Btemp=BBmn
		Ctemp=CCmn
		Dtemp=DDmn
	
		do step=1,4			!call the 4 steps of the rk4 routine
			AI=Atemp+odestep(step)*Aprime
			BI=Btemp+odestep(step)*Bprime
			CI=Ctemp+odestep(step)*Cprime
			DI=Dtemp+odestep(step)*Dprime
						
			call SeismicRHS(AI,BI,CI,DI,Cprime,Dprime,Cx,Sx,Cy,Sy,GmCx,GnCy,alphax,alphay,nxx,nyy,trapwt,Mmax,Nmax,damping,rhoinv,mu,&
					lambda,gm,gn,mux,muy,lambdax,lambday,time)
			
			Aprime=CI
			Bprime=DI
			
			AAmn=AAmn+odewt(step)*Aprime
			BBmn=BBmn+odewt(step)*Bprime
			CCmn=CCmn+odewt(step)*Cprime
			DDmn=DDmn+odewt(step)*Dprime
			
		end do
		looptime=looptime+h
		stepnum=stepnum+1
		
		!!$acc kernels
		do n1=1,numsensors
			sumvx=0.
			sumvy=0.
			m1=sensorsi(n1)
			do n=0,Nmax
				do m=0,Mmax
					sumvx=sumvx+CCmn(m,n)*Cx(m,m1)*Cy(n,sensorj)
					sumvy=sumvy+DDmn(m,n)*Cx(m,m1)*Cy(n,sensorj)
				enddo
			end do
			vxout(stepnum,n1)=sumvx
			vyout(stepnum,n1)=sumvy
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
call printelapsedtime(day,hours,minutes,seconds,milliseconds)
		
call savearray('data/vx.out',len('data/vx.out'),vxout,totalsteps,numsensors)
call savearray('data/vy.out',len('data/vy.out'),vyout,totalsteps,numsensors)

end do
call printenddatetime()

print*,'time1:=',time(1)*1000.
print*,'time2:=',time(2)*1000.
print*,'time3:=',time(3)*1000.
print*,'time4:=',time(4)*1000.

print*,'nxx: ',nxx
end program SeismicIntegrator

Subroutine SeismicRHS(AI,BI,CI,DI,Cprime,Dprime,Cx,Sx,Cy,Sy,GmCx,GnCy,alphax,alphay,nxx,nyy,trapwt,Mmax,Nmax,damping,rhoinv,mu,&
			lambda,gm,gn,mux,muy,lambdax,lambday,time)
!This routine evaluates the right-hand sides of the four sets of
!ordinary differential equations for the Fourier coefficients
use mmul_mod
implicit none
Integer,intent(in)::nxx,nyy,Mmax,Nmax
real,Dimension(nxx,nyy),intent(in)::damping,rhoinv,mu,lambda,mux,muy,lambdax,lambday,trapwt
!trig basis functions
real,Dimension(0:Mmax,nxx),intent(in)::Cx,Sx
real,Dimension(0:Nmax,nyy),intent(in)::Cy,Sy
real,Dimension(0:Mmax,nxx)::gmCx
real,Dimension(0:Nmax,nyy)::gnCy
real,Dimension(nxx,0:Nmax)::C
!Arrays for the Fourier coefficients.
real,Dimension(0:Mmax,0:Nmax),intent(in)::AI,BI,CI,DI
real,Dimension(0:Mmax,0:Nmax),intent(out)::Cprime,Dprime
real,Dimension(0:Mmax,0:Nmax)::mAI,nAI,mBI,nBI,mmA25,mnA25,nnA25,mmB25,mnB25,nnB25
real,Dimension(0:Mmax)::gm
real,Dimension(0:Nmax)::gn
!dummy spatial variables and counters
real,Dimension(nxx,nyy)::uxxx,uxxy,uxyy,uyyy,uyxy,uyxx,ux,uy,vx,vy,uxx,uxy,uyx,uyy

integer::i,j,m,n
real,intent(in)::alphax,alphay
real,dimension(4),intent(inout)::time
real::wt
real :: start, finish
!  calculate the spatial variables from the fourier coefficients
!	first, vectorise all the m,n terms into matrices
call cpu_time(start)
do m=0,Mmax
    do n=0,Nmax
		mmA25(m,n)=m*m*AI(m,n)*alphax*alphax
		mnA25(m,n)=m*n*AI(m,n)*alphax*alphay
		nnA25(m,n)=n*n*AI(m,n)*alphay*alphay
		mmB25(m,n)=m*m*BI(m,n)*alphax*alphax
		mnB25(m,n)=m*n*BI(m,n)*alphax*alphay
		nnB25(m,n)=n*n*BI(m,n)*alphay*alphay
		mAI(m,n)=AI(m,n)*m*alphax
		nAI(m,n)=AI(m,n)*n*alphay
		mBI(m,n)=BI(m,n)*m*alphax
		nBI(m,n)=BI(m,n)*n*alphay
	end do
end do
call cpu_time(finish)
time(1)=time(1)+finish-start
!	now use two matrix multiply operations to determine each spatial variable
call cpu_time(start)
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
call mmul( transpose(sx), -mAI, C )
call mmul(C,cy,uxx)
call mmul( transpose(cx), -nAI, C )
call mmul(C,sy,uxy)
call mmul( transpose(sx), -mBI, C )
call mmul(C,cy,uyx)
call mmul( transpose(cx), -nBI, C )
call mmul(C,sy,uyy)
call mmul( transpose(cx), CI, C )
call mmul(C,cy,vx)
call mmul( transpose(cx), DI, C )
call mmul(C,cy,vy)
call cpu_time(finish)
time(2)=time(2)+finish-start
!	calculate the displacement from the elastic wave equation, with damping
call cpu_time(start)
ux=rhoinv*((lambda+2.*mu)*uxxx+(lambda+mu)*uyxy+mu*uxyy+lambdax*(uxx+uyy)+2.*mux*uxx+muy*(uyx+uxy))-vx*damping
uy=rhoinv*((lambda+2.*mu)*uyyy+(lambda+mu)*uxxy+mu*uyxx+lambday*(uxx+uyy)+2.*muy*uyy+mux*(uyx+uxy))-vy*damping
call cpu_time(finish)
time(3)=time(3)+finish-start
!	finally, calculate the time derivatives of the coefficients from ux'' and uy''
call cpu_time(start)
call mmul(  trapwt*ux,transpose(gncy), C )
call mmul(gmcx,C,Cprime)
call mmul(  trapwt*uy,transpose(gncy), C )
call mmul(gmcx,C,Dprime)
call cpu_time(finish)
time(4)=time(4)+finish-start
end Subroutine SeismicRHS
