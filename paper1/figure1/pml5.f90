program SeismicIntegrator
!This routine integrates the systems of ODEs for the 2-d elastic equations

!pgfortran pml5.f90 support.o mmul_mod.o -o pml1 -O3 -acc -gpu=cc60 -r8 -cuda

use support	! pgfortran -c support.f90 -r8
implicit none
real, parameter::pi=3.14159265358979323846
integer, parameter::nxout=7*80,nyout=nxout,outratio=1,nx=outratio*nxout,ny=outratio*nyout,dampwidth=nx/20

Integer,parameter::nimages=4				!-----number of times to output the density for an image---!
Real,parameter::t=4./nimages,h=2./nxout		!-----time units per image-----time for initial pass through the RK4 routine-----!
real,parameter,Dimension(4)::rkpts=(/0.,h/2,h/2,h/),rkwts=(/h/6.,h/3.,h/3.,h/6./)	!rk4 spacing and weighting

Integer,parameter::loops=NInt(t/h),totalsteps=loops*nimages

integer, parameter::Mmax=nx/5-1,Nmax=ny/5-1

Integer,parameter::numsensors=nx/10+1,sensorj=ny-dampwidth-10		!-----location of sensor-----!
integer,Dimension(numsensors)::sensorsi
real,Dimension(numsensors)::sensorsx
integer,parameter::ds=(nx-2*dampwidth)/(numsensors-1)

real,Dimension(totalsteps,numsensors)::uxout,uyout

real,Dimension(nx)::x
real,Dimension(ny)::y
real,Dimension(nxout)::xout
real,Dimension(nyout)::yout
real,Dimension(nx,ny)::u,ux,uy,damping,rhoinv,mu,lambda,trapwt,alphax,alphaxx,alphay,alphayy,betax,betay,betaxx,betayy
real,Dimension(nxout,nyout)::uout,betaxout,betaxxout
real,Dimension(0:Mmax,nx)::Cx,Sx,gmCx
real,Dimension(nx,0:Mmax)::CxT,SxT
real,Dimension(0:Nmax,ny)::Cy,Sy,gnCy
real,Dimension(ny,0:Nmax)::gnCyTrans
real,Dimension(0:Mmax)::gm
real,Dimension(0:Nmax)::gn

!arrays for the Fourier coefficients.
real,Dimension(0:Mmax,0:Nmax)::Amn,Atemp,AI,Aprime,Bmn,Btemp,BI,Bprime,Cmn,Ctemp,CI,Cprime,Dmn,Dtemp,DI,Dprime
real,Dimension(0:Mmax,0:Nmax)::E11,E12,E21,E22,E11i,E12i,E21i,E22i,E11p,E12p,E21p,E22p,E11temp,E12temp,E21temp,E22temp

real::argmpx,argnpy,wt,sumux,sumuy,looptime,dampdepth,r	!---time through the RKF routine for each image---!
real,parameter::a=0.1,F0=1,G0=1,df=5.4,df2=3.42,x0=1.,y0=1.,betamax=60,betaindex=1.6
real,parameter::x1=-pi,x2=-x1,y1=x1,y2=-y1

real,parameter::vp1=1.5,vs1=0.5,mu1=vs1**2,lambda1=vp1**2-2*vs1**2
integer::m,n,i,j,rk4step,image,m1,n1,iLoop,stepnum=0
character(16)::filename
real::dx,dy,dw,alphaindex
real::wtsq=(0.5/x2*pi)**2
!-----timing  variables
integer, dimension(8) :: time_array

!-----initial velocity field
Cmn=0.
Dmn=0.

E11=0.
E12=0.
E21=0.
E22=0.
!-----Print the start time to screen-----!
call tic(time_array)
call printdatetime('begun')
print*, 'dampwidth: ',dampwidth
print*, 'zeta: ',betaindex
print*, 'betamax: ',betamax

print*, 'num images: ',nimages
print*,'loops per image: ',loops
print*,'resolution: ',nx,' by ',ny
call linspace(x, dx,x1,x2, nx)
call linspace(y, dy,y1,y2, ny)

call savevector('x.out',x,nx)
call savevector('y.out',y,ny)

do i=1,nxout
	xout(i)=x(i*outratio)
enddo 
do j=1,nyout
	yout(j)=y(j*outratio)
enddo 

call savevector('xout.out',xout,nxout)
call savevector('yout.out',yout,nyout)

call saveinteger('dampwidth.out',dampwidth)
call saveinteger('nimages.out',nimages)

call savereal('t.out',t)
call saveinteger('totalsteps.out',totalsteps)
call saveinteger('numsensors.out',numsensors)

call savevector('sensorx.out',x,nx)
call savereal('sensory.out',y(sensorj))

trapwt=1.
do i=1,nx
	trapwt(i,1)=0.5
	trapwt(i,ny)=0.5
enddo
do j=1,ny
	trapwt(1,j)=0.5
	trapwt(nx,j)=0.5
enddo
trapwt(1,1)=0.25
trapwt(nx,ny)=0.25
trapwt(1,ny)=0.25
trapwt(nx,1)=0.25

trapwt=trapwt*dx*dy/x1**2

dw=real(dampwidth)
open(unit=2,file='damping.txt')
	write(2,*, ADVANCE = "YES") df
	write(2,*, ADVANCE = "YES") df2
close(2)

open(unit=2,file='details.txt')
	write(2,*, ADVANCE = "YES") 'damping width dw = ',dw
	write(2,*, ADVANCE = "YES") 'Mmodes = ',Mmax
	write(2,*, ADVANCE = "YES") 'Nmodes = ',Nmax
	write(2,*, ADVANCE = "YES") 'x points nx= ',nx
	write(2,*, ADVANCE = "YES") 'y points ny= ',ny	
	write(2,*, ADVANCE = "YES") 'num images= ',nimages	
	write(2,*, ADVANCE = "YES") 'time per image= ',t	
	write(2,*, ADVANCE = "YES") 'timestep= ',h	
close(2)

!locate the sensors
do i=1,numsensors
	
	sensorsi(i)=dampwidth+ds*(i-1)
	sensorsx(i)=x(sensorsi(i))
enddo
call savevector('sensorsx.out',sensorsx,numsensors)
!create the damping fields for absorbing boundaries
damping=0.

!acc kernels
do i=1,dampwidth
	dampdepth=(1-real(i)/dampwidth)
	do j=1,ny
		damping(i,j)=exp(df*dampdepth)-1.
		damping(nx+1-i,j)=exp(df*dampdepth)-1.
	enddo
enddo
!acc end kernels

do j=1,dampwidth
	dampdepth=(1.-real(j)/dampwidth)
	do i=1,nx
		damping(i,j)=max(damping(i,j),exp(df*dampdepth)-1.)
		damping(i,ny+1-j)=max(damping(i,ny+1-j),exp(df*dampdepth)-1.)
	enddo
enddo
damping=df2*damping/df
call savearray('damping.out',damping,nx,ny)
!Create and store the trig functions  cos(m/2(pi+x))  and  sin(m/2(pi+x)) .
do m=0,Mmax
    do i=1,nx
        argmpx=(x(i)-x1)*m/2./x2*pi
        Cx(m,i)=cos(argmpx)
        Sx(m,i)=sin(argmpx)
    end do
end do
CxT=Transpose(Cx)
SxT=Transpose(Sx)
do n=0,Nmax
    do j=1,ny
        argnpy=(y(j)-y1)*n/2./y2*pi
        Cy(n,j)=cos(argnpy)
        Sy(n,j)=sin(argnpy)
    end do
end do

!$acc kernels
do i=1,nx
	do j=1,ny
		ux(i,j)=F0*2.*(x(i)-x0)/a**2*exp(-((x(i)-x0)**2+(y(j)-y0)**2)/a**2)
		uy(i,j)=G0*2.*(y(j)-y0)/a**2*exp(-((x(i)-x0)**2+(y(j)-y0)**2)/a**2)
   end do
end do
!$acc end kernels
u=sqrt(ux**2+uy**2)
call savearray('data/disp00.out',u,nx,ny)

!-----Integration weightings
gm=1.
gm(0)=0.5
gn=1.
gn(0)=0.5

Amn=0.
Bmn=0.
!$acc kernels
do m=0,Mmax
    do n=0,Nmax!
	    do i=1,nx
            do j=1,ny
				wt=gm(m)*gn(n)*Cx(m,i)*Cy(n,j)
                Amn(m,n)=Amn(m,n)+wt*ux(i,j)*trapwt(i,j)
                Bmn(m,n)=Bmn(m,n)+wt*uy(i,j)*trapwt(i,j)
            end do
        end do
    end do
end do
!$acc end kernels

ux=matmul(matmul(CxT,Amn),cy)
uy=matmul(matmul(CxT,Bmn),cy)

u=sqrt(ux**2+uy**2)
do i=1,nxout
	do j=1,nyout
		uout(i,j)=u(i*outratio,j*outratio)
	enddo
enddo
call savearray('data/disp0.out',uout,nxout,nyout)

!set up the density distribution
rhoinv=1.
mu=mu1
lambda=lambda1

do m=0,Mmax
	do i=1,nx
		gmcx(m,i)=gm(m)*Cx(m,i)
	enddo
enddo

do n=0,Nmax
	do j=1,ny
		gncytrans(j,n)=gn(n)*Cy(n,j)
	enddo
enddo

betax=0.
betay=0.
betaxx=0.
betayy=0.

do i=1,dampwidth
	do j=1,ny
		betax(i,j)=((x(dampwidth)-x(i))/(x(dampwidth)-x(1)))**betaindex
		betax(nx+1-i,j)=betax(i,j)
		betaxx(i,j)=-betaindex*((x(dampwidth)-x(i))**(betaindex-1.))/(x(dampwidth)-x(1))**betaindex
		betaxx(nx+1-i,j)=-betaxx(i,j)
	enddo
enddo
call savearray('data/betax.out',betaxx,nx,ny)
do j=1,dampwidth
	do i=1,nx
		betay(i,j)=((y(dampwidth)-y(j))/(y(dampwidth)-y(1)))**betaindex
		betay(i,ny+1-j)=betay(i,j)
		betayy(i,j)=-betaindex*((y(dampwidth)-y(j))**(betaindex-1))/(y(dampwidth)-y(1))**betaindex
		betayy(i,ny+1-j)=-betayy(i,j)
	enddo
enddo

betax=betamax*betax
betaxx=betamax*betaxx
betay=betamax*betay
betayy=betamax*betayy

do i=1,nxout
	do j=1,nyout
		betaxout(i,j)=betax(i*outratio,j*outratio)
		betaxxout(i,j)=betaxx(i*outratio,j*outratio)
	enddo
enddo
call savearray('data/betaxout.out',betaxout,nxout,nyout)
call savearray('data/betaxxout.out',betaxxout,nxout,nyout)

rhoinv=rhoinv*trapwt
damping=damping*trapwt
!Do the integration process.
do image=1,nimages
	print*,'interval ',image,' started'
	looptime=0.
	do iLoop=1,loops
		Atemp=Amn
		Btemp=Bmn
		Ctemp=Cmn
		Dtemp=Dmn
		E11temp=E11
		E12temp=E12
		E21temp=E21
		E22temp=E22
	
		do rk4step=1,4			!call the 4 steps of the rk4 routine
			AI=Atemp+rkpts(rk4step)*Aprime
			BI=Btemp+rkpts(rk4step)*Bprime
			CI=Ctemp+rkpts(rk4step)*Cprime
			DI=Dtemp+rkpts(rk4step)*Dprime
			E11I=E11temp+rkpts(rk4step)*E11p
			E12I=E12temp+rkpts(rk4step)*E12p
			E21I=E21temp+rkpts(rk4step)*E21p
			E22I=E22temp+rkpts(rk4step)*E22p
						
			call SeismicRHS(betax,betay,betaxx,betayy,E11,E12,E21,E22,E11p,E12p,E21p,E22p,AI,BI,CI,DI,Cprime,Dprime,CxT,SxT,Cy,Sy,GmCx,GnCyTrans,nx,ny,Mmax,Nmax,damping,rhoinv,mu,lambda,wtsq)
			
			Aprime=CI
			Bprime=DI
			
			Amn=Amn+rkwts(rk4step)*Aprime
			Bmn=Bmn+rkwts(rk4step)*Bprime
			Cmn=Cmn+rkwts(rk4step)*Cprime
			Dmn=Dmn+rkwts(rk4step)*Dprime
			E11=E11+rkwts(rk4step)*E11p
			E12=E12+rkwts(rk4step)*E12p
			E21=E21+rkwts(rk4step)*E21p
			E22=E22+rkwts(rk4step)*E22p
			
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
					sumux=sumux+Amn(m,n)*Cx(m,m1)*Cy(n,sensorj)
					sumuy=sumuy+Bmn(m,n)*Cx(m,m1)*Cy(n,sensorj)
				enddo
			end do
			uxout(stepnum,n1)=sumux
			uyout(stepnum,n1)=sumuy
		end do
		!!$acc end kernels
	end do
	print*,'looptime: ',looptime

!Compute the Displacement function and save it.
ux=matmul(matmul(transpose(cx),Amn),cy)
uy=matmul(matmul(transpose(cx),Bmn),cy)
u=sqrt(ux**2+uy**2)
do i=1,nxout
	do j=1,nyout
		uout(i,j)=u(i*outratio,j*outratio)
	enddo
enddo
print*,'stepnum: ',stepnum
print *,'max displacement=  ',maxval(uout(dampwidth:nxout+1-dampwidth,dampwidth:nyout+1-dampwidth))

	if (image<10) then
		filename='data/disp'//char(48+image)//'.out'		!--------disp0 to disp9------!
	elseif (image<100) then
		filename='data/disp'//char(48+image/10)//char(48+image-10*(image/10))//'.out'		!--------disp10 to disp99----!
	elseif (image<1000) then
		filename='data/disp'//char(48+image/100)//char(48+image/10-10*(image/100))//char(48+image-10*(image/10))//'.out'	!--------disp100 to disp999----!
	endif
!call savearray(filename,u,nx,ny)
call savearray(filename,uout,nxout,nyout)	
call savearray('data/ux.out',uxout,totalsteps,numsensors)
call savearray('data/uy.out',uyout,totalsteps,numsensors)
print*,filename
call toc(time_array)
	
end do

call printdatetime('finito')


end program SeismicIntegrator

Subroutine SeismicRHS(bx,by,bxx,byy,E11,E12,E21,E22,E11p,E12p,E21p,E22p,AI,BI,CI,DI,Cprime,Dprime,CxT,SxT,Cy,Sy,GmCx,GnCyTrans,nx,ny,Mmax,Nmax,damping,rhoinv,mu,lambda,wtsq)
!This routine evaluates the right-hand sides of the four sets of
!ordinary differential equations for the Fourier coefficients
use mmul_mod
implicit none
Integer,intent(in)::nx,ny,Mmax,Nmax
real,Dimension(nx,ny),intent(in)::damping,rhoinv,mu,lambda,bx,bxx,by,byy
!trig basis functions
real,Dimension(0:Mmax,nx),intent(in)::GmCx
real,Dimension(nx,0:Mmax),intent(in)::CxT,SxT
real,Dimension(0:Nmax,ny),intent(in)::Cy,Sy
real,Dimension(ny,0:Nmax),intent(in)::GnCyTrans
!Arrays for the Fourier coefficients.
real,Dimension(0:Mmax,0:Nmax),intent(in)::AI,BI,CI,DI,e11,e12,e21,e22

real,Dimension(0:Mmax,0:Nmax),intent(out)::Cprime,Dprime,E11p,E12p,E21p,E22p
real,Dimension(0:Mmax,0:Nmax)::mmA25,mnA25,nnA25,mmB25,mnB25,nnB25
real,Dimension(0:Mmax,0:Nmax)::mE11,nE11,mE12,nE12,mE21,nE21,mE22,nE22,ma5,na5,mb5,nb5

!dummy spatial variables and counters
real,Dimension(nx,ny)::uxxx,uxxy,uxyy,uyyy,uyxy,uyxx,ux,uy,vx,vy,uxx,uyy,uxy,uyx
real,Dimension(nx,ny)::w11,w12,w21,w22,w11t,w12t,w21t,w22t,w11x,w12y,w21x,w22y
integer::m,n,i,j
real,intent(in)::wtsq
!  calculate the spatial variables from the fourier coefficients
do m=0,Mmax
    do n=0,Nmax
		mmA25(m,n)=m*m*AI(m,n)*wtsq
		mnA25(m,n)=m*n*AI(m,n)*wtsq
		nnA25(m,n)=n*n*AI(m,n)*wtsq
		mmB25(m,n)=m*m*BI(m,n)*wtsq
		mnB25(m,n)=m*n*BI(m,n)*wtsq
		nnB25(m,n)=n*n*BI(m,n)*wtsq
		mE11(m,n)=-m*E11(m,n)*0.5
		nE11(m,n)=-n*E11(m,n)*0.5
		mE12(m,n)=-m*E12(m,n)*0.5
		nE12(m,n)=-n*E12(m,n)*0.5
		mE21(m,n)=-m*E21(m,n)*0.5
		nE21(m,n)=-n*E21(m,n)*0.5
		mE22(m,n)=-m*E22(m,n)*0.5
		nE22(m,n)=-n*E22(m,n)*0.5
		mA5(m,n)=-m*AI(m,n)*0.5
		nA5(m,n)=-n*AI(m,n)*0.5
		mB5(m,n)=-m*BI(m,n)*0.5
		nB5(m,n)=-n*BI(m,n)*0.5
	end do
end do

call mmul26(CxT, -mmA25,-nnA25,-mmB25,-nnB25,CI,DI,cy,uxxx,uxyy,uyxx,uyyy,vx,vy)
call mmul26(CxT, AI,BI,E11,E12,E21,E22,cy,ux,uy,w11,w12,w21,w22)
call mmul22(SxT, mnA25,mnB25,sy,uxxy,uyxy)
call mmul24(SxT,mE11,mE21,mA5,mB5,cy,w11x,w21x,uxx,uyx)
call mmul24(CxT,nE12,nE22,nA5,nB5,sy,w12y,w22y,uxy,uyy)
w11t=rhoinv*((lambda+2.*mu)*uxx-bx*w11)
w12t=rhoinv*(mu*uxy-by*w12)
w21t=rhoinv*(mu*uyx-bx*w21)
w22t=rhoinv*((lambda+2.*mu)*uyy-by*w22)

ux=rhoinv*((lambda+2.*mu)*(uxxx)+(lambda+mu)*uyxy+mu*(uxyy)+(by-bx)*w11x-bxx*w11+(bx-by)*w12y-byy*w12-(bx+by)*vx-bx*by*ux)
uy=rhoinv*((lambda+2.*mu)*(uyyy)+(lambda+mu)*uxxy+mu*(uyxx)+(bx-by)*w22y-byy*w22+(by-bx)*w21x-bxx*w21-(bx+by)*vy-bx*by*uy)

call mmul26(gmcx,ux,uy,w11t,w12t,w21t,w22t,GnCyTrans,Cprime,Dprime,E11p,E12p,E21p,E22p)
end Subroutine SeismicRHS
