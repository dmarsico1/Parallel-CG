program mpi_cart_poisson_3d

use mpi

implicit none


integer :: M=300
integer :: My=300
integer :: i,j,l
integer :: N=300
integer :: ierr,total_sum,myid,k,z_dim,x_dim,y_dim
integer :: status(mpi_status_size)
real*8, dimension(:), allocatable :: x,y,z
real*8, dimension(:,:,:), allocatable :: f,r,p,q,z1
real*8 :: hx,hz,hy,rold,rnew,rnew_tot,rold_tot,alpha_tot,gam,gam_tot,tol
character(len=100) :: str_form,temp_str
real*8 :: start,finish
integer :: comm2d
integer, dimension(3) :: outdims,outcoords,coords
integer, dimension(3) ::  dims
logical :: reorder 
logical, dimension(3) :: isperiodic,outperiods
integer :: ndim

dims(1)=0
dims(2)=0
dims(3)=0

isperiodic(1)=.false.
isperiodic(2)=.false.
isperiodic(3)=.false.
reorder=.true.
ndim=3

hx=8.0*atan(1.0)/(N+1.0)
hy=8.0*atan(1.0)/(My+1.0)
hz=8.0*atan(1.0)/(M+1.0)

tol=1.0/1000.0


call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,k,ierr)

call mpi_dims_create(k,ndim,dims,ierr)

call mpi_cart_create(mpi_comm_world,ndim,dims,isperiodic,reorder,comm2d,ierr)

call mpi_cart_get(comm2d,ndim,outdims,outperiods,outcoords,ierr)

call mpi_comm_rank(comm2d,myid,ierr)


if(outcoords(1)==0)then
	z_dim=M-(outdims(1)-1)*(M/outdims(1))
else
	z_dim=M/outdims(1)
endif

if(outcoords(2)==0)then
	x_dim=N-(outdims(2)-1)*(N/outdims(2))
else
	x_dim=N/outdims(2)
endif

if(outcoords(3)==0)then
	y_dim=My-(outdims(3)-1)*(My/outdims(3))
else
	y_dim=My/outdims(3)
endif

allocate(f(z_dim,x_dim,y_dim),x(x_dim),y(y_dim),z(z_dim),z1(z_dim,x_dim,y_dim),q(z_dim,x_dim,y_dim),&
		r(z_dim,x_dim,y_dim),p(z_dim+2,x_dim+2,y_dim+2))


if(outcoords(1)==0)then
	do i=1,z_dim
		z(i)=dble(i)*hz
	enddo
else
	do i=1,z_dim
		z(i)=(M-(outdims(1)-1)*(M/outdims(1)))*hz+z_dim*(outcoords(1)-1)*hz+dble(i)*hz
	enddo
endif

if(outcoords(2)==0)then
	do i=1,x_dim
		x(i)=dble(i)*hx
	enddo
else
	do i=1,x_dim
		x(i)=(N-(outdims(2)-1)*(N/outdims(2)))*hx+x_dim*(outcoords(2)-1)*hx+dble(i)*hx
	enddo
endif

if(outcoords(3)==0)then
	do i=1,y_dim
		y(i)=dble(i)*hy
	enddo
else
	do i=1,y_dim
		y(i)=(My-(outdims(3)-1)*(My/outdims(3)))*hy+y_dim*(outcoords(3)-1)*hy+dble(i)*hy
	enddo
endif


call forcing(f,x,y,z,z_dim,x_dim,y_dim)

p=0.0

r=f
p(2:z_dim+1,2:x_dim+1,2:y_dim+1)=r

call cpu_time(start)

rold=sum(r*r)
call mpi_allreduce(rold,rold_tot,1,mpi_real8,mpi_sum,comm2d,ierr)
do l=1,N
	call neighbors(myid,k)
	call matmult(q,p,z_dim,x_dim,y_dim,hx,hz,hy)
	gam=sum(q*p(2:z_dim+1,2:x_dim+1,2:y_dim+1))
	call mpi_allreduce(gam,gam_tot,1,mpi_real8,mpi_sum,comm2d,ierr)
	alpha_tot=rold_tot/gam_tot
	z1=z1+alpha_tot*p(2:z_dim+1,2:x_dim+1,2:y_dim+1)
	r=r-alpha_tot*q
	rnew=sum(r*r)
	call mpi_allreduce(rnew,rnew_tot,1,mpi_real8,mpi_sum,comm2d,ierr)
	if(sqrt(rnew_tot)<tol)then
		exit
	endif
	p(2:z_dim+1,2:x_dim+1,2:y_dim+1) = r + (rnew_tot/rold_tot)*p(2:z_dim+1,2:x_dim+1,2:y_dim+1)
	rold_tot=rnew_tot
	call mpi_barrier(comm2d,ierr)
enddo

call cpu_time(finish)

if(myid==0)then
	write(*,*) sqrt(rnew_tot),l,finish-start
endif


if(myid<10)then
	str_form='(I1)'
elseif(myid>=10)then
	str_form='(I2)'
endif

write(temp_str,str_form) myid

if(myid==4)then
open(unit=12,status='replace',file='test'//'1'//'.txt',form='formatted')
do i=1,x_dim
	write(12,*) (z1(12,i,j), j=1,y_dim)
enddo
close(12)
endif

call mpi_finalize(ierr)

contains

subroutine forcing(f,x,y,z,M,N,My)

implicit none

integer :: M,N,My
real*8, intent(in), dimension(N) :: x
real*8, intent(in), dimension(M) :: z
real*8, intent(in), dimension(My) :: y
real*8, intent(out), dimension(M,N,My) :: f
integer :: i,j,k

do i=1,M
	do j=1,N
		do k=1,My
			f(i,j,k)=4.0*sin(x(j))*sin(z(i))*sin(y(k))-sin(2.0*x(j))*sin(2.0*z(i))*sin(2.0*y(k))&
					+3.0*sin(4.0*x(j))*sin(4.0*z(i))*sin(4.0*y(k))
		enddo
	enddo
enddo

end subroutine


subroutine neighbors(myid,k)

implicit none

integer :: right,left,above,below,front,back,myid,k

call mpi_cart_shift(comm2d,0,1,below,above,ierr)

call mpi_cart_shift(comm2d,1,1,left,right,ierr)

call mpi_cart_shift(comm2d,2,1,back,front,ierr)


call mpi_sendrecv(p(z_dim+1,2:x_dim+1,2:y_dim+1),x_dim*y_dim,mpi_double,above,0,&
				p(1,2:x_dim+1,2:y_dim+1),x_dim*y_dim,mpi_double,below,0,&
				comm2d,status,ierr)

call mpi_sendrecv(p(2,2:x_dim+1,2:y_dim+1),x_dim*y_dim,mpi_double,below,1,&
				p(z_dim+2,2:x_dim+1,2:y_dim+1),x_dim*y_dim,mpi_double,above,1,&
				comm2d,status,ierr)


call mpi_sendrecv(p(2:z_dim+1,x_dim+1,2:y_dim+1),z_dim*y_dim,mpi_double,right,2,&
				p(2:z_dim+1,1,2:y_dim+1),z_dim*y_dim,mpi_double,left,2,&
				comm2d,status,ierr)
			
call mpi_sendrecv(p(2:z_dim+1,2,2:y_dim+1),z_dim*y_dim,mpi_double,left,3,&
				p(2:z_dim+1,x_dim+2,2:y_dim+1),z_dim*y_dim,mpi_double,right,3,&
				comm2d,status,ierr)

call mpi_sendrecv(p(2:z_dim+1,2:x_dim+1,y_dim+1),z_dim*x_dim,mpi_double,front,4,&
				p(2:z_dim+1,2:x_dim+1,1),z_dim*x_dim,mpi_double,back,4,&
				comm2d,status,ierr)
			
call mpi_sendrecv(p(2:z_dim+1,2:x_dim+1,2),z_dim*x_dim,mpi_double,back,5,&
				p(2:z_dim+1,2:x_dim+1,y_dim+2),z_dim*x_dim,mpi_double,front,5,&
				comm2d,status,ierr)

end subroutine neighbors


subroutine MatMult(p,q,M,N,My,hx,hz,hy)

implicit none
real*8, dimension(M+2,N+2,My+2), intent(in) :: q !this is what gets INPUT.  Layer of ghost cells around q.
real*8, dimension(M,N,My), intent(out) :: p !this is what gets OUTPUT.  No ghost cells around p.
real*8, intent(in) :: hx,hz,hy
integer, intent(in) :: M,N,My
integer :: i,j,k

do i=2,M+1
	do j=2,N+1
		do k=2,My+1
			p(i-1,j-1,k-1)=(1.0/hx**2)*(q(i,j+1,k)-2.0*q(i,j,k)+q(i,j-1,k))+(1.0/hz**2)*(q(i+1,j,k)-2.0*q(i,j,k)+q(i-1,j,k))&
						+(1.0/hy**2)*(q(i,j,k+1)-2.0*q(i,j,k)+q(i,j,k-1))
		enddo
	enddo
enddo


end subroutine MatMult


end program mpi_cart_poisson_3d

