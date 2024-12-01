module qfactor_common
implicit none
	integer:: qx, qy, qz, q1, q2, q1m1, q2m1, Normal_index
	integer(1), allocatable:: rsboundary(:,:), reboundary(:,:), rboundary_tmp(:, :), sign2d(:,:)
	real:: xreg(0:1), yreg(0:1), zreg(0:1), cut_coordinate, delta, &
	point0(0:2), ev1(0:2), ev2(0:2), ev3(0:2)
	real, allocatable:: rsF(:, :, :), reF(:, :, :), bnr(:, :), &
	q_perp(:,:), q(:, :), length(:, :), twist(:, :)
	logical:: twistFlag, vflag, q0flag, cflag, csflag, scottFlag
	logical, allocatable:: tangent_Flag(:, :)
end module qfactor_common


module trace_common
implicit none
	integer:: nx, ny, nz, nxm1, nym1, nzm1, nxm2, nym2, nzm2, maxsteps, maxsteps_foot, r0max(0:2), &
	binary_index_xs, binary_index_ys, binary_index_zs, index_try_xs, index_try_ys, index_try_zs
	integer(2), allocatable:: binary_values(:)
	real:: xmax, ymax, zmax, xmin, ymin, zmin, pmin(0:2), pmax(0:2), &
	step, min_step, min_step_foot, tol, min_incline, NaN, dxa_uni, dya_uni, dza_uni
	real(8), parameter:: pi=3.141592653589793D0
	logical:: RK4flag, grad3DFlag, stretchFlag, uni_stretch_Flag
	real, allocatable:: xa(:), ya(:), za(:), dxa(:), dya(:), dza(:)
!----------------------------------------------------------------------------
! for RKF45
	real:: c2,c3,c4,c5,c6, &
	a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65, &
	b1,b3,b4,b5,b6, &
	ce1,ce3,ce4,ce5,ce6;
!----------------------------------------------------------------------------
INTERFACE
	subroutine indexes(vp, vpBound, index_i, index_j, index_k)
	implicit none
	real:: vp(0:2), vpBound(0:2)
	integer:: index_i, index_j, index_k
	END  subroutine indexes
	
	subroutine weights(vp, round, weight)
	implicit none
	real:: vp(0:2), weight(0:1,0:1,0:1)
	integer:: round(0:1,0:2)
	END  subroutine weights
END INTERFACE
!----------------------------------------------------------------------------
procedure (indexes),  pointer :: vp_index     => null ()
procedure (weights),  pointer :: round_weight => null ()
end module trace_common


module field_common
implicit none
real, allocatable:: Bfield(:, :, :, :), CurlB(:,:, :, :), grad_unit_vec_Bfield(:, :, :, :, :)
!----------------------------------------------------------------------------
INTERFACE
	subroutine calculate_curlB_grid(i, j, k, CurlBp)
	implicit none
	integer:: i, j, k
	real:: CurlBp(0:2)
	END  subroutine calculate_curlB_grid

	subroutine calculate_grad_unit_vec_B_grid(i, j, k, grad_unit_vec_B)
	implicit none
	integer:: i, j, k
	real:: grad_unit_vec_B(0:2,0:2)
	END  subroutine calculate_grad_unit_vec_B_grid
END INTERFACE
procedure (calculate_curlB_grid),           pointer :: curlB_grid           => null ()
procedure (calculate_grad_unit_vec_B_grid), pointer :: grad_unit_vec_B_grid => null ()
end module field_common


subroutine vp_index0(vp, vpBound, index_i, index_j, index_k)
use trace_common
implicit none
real:: vp(0:2), vpBound(0:2), xp, yp, zp
integer:: i, index_i, index_j, index_k, binary_index, index_try
!----------------------------------------------------------------------------
vpBound=vp
do i=0,2
	if ( .not. (vpBound(i) .ge. pmin(i))) then 
	! this way can avoid the crash of vp(i) .eq. NaN (caused by B=0), compared with vp(i) .lt. 0.0
		vpBound(i)=pmin(i)
	else if (vpBound(i) .ge. pmax(i)) then
		vpBound(i)=pmax(i)
	endif
enddo
!----------------------------------------------------------------------------
! this way is slower than the later way
!index_i=count(vpBound(0) .ge. xa(0:nxm2))-1
!index_j=count(vpBound(1) .ge. ya(0:nym2))-1
!index_k=count(vpBound(2) .ge. za(0:nzm2))-1
!----------------------------------------------------------------------------
!binary (tree) search
xp=vpBound(0)
yp=vpBound(1)
zp=vpBound(2)
!----------------------------------------------------------------------------
binary_index=binary_index_xs
index_try=index_try_xs  !binary_values(binary_index_xs)
do while(binary_index .ge. 1)
	binary_index=binary_index-1
	if (xp .ge. xa(index_try))  then
		if (index_try+binary_values(binary_index) .le. nxm2) &
		index_try=index_try+binary_values(binary_index)		
	else		
		index_try=index_try-binary_values(binary_index)		
	endif
enddo

if ((xp .ge. xa(index_try)) .and. (index_try .le. nxm2)) then
	index_i=index_try
else
	index_i=index_try-1
endif
!----------------------------------------------------------------------------
binary_index=binary_index_ys
index_try=index_try_ys
do while(binary_index .ge. 1)
	binary_index=binary_index-1
	if (yp .ge. ya(index_try))  then
		if (index_try+binary_values(binary_index) .le. nym2) &
		index_try=index_try+binary_values(binary_index)
	else
		index_try=index_try-binary_values(binary_index)		
	endif
	
enddo

if (yp .ge. ya(index_try) .and. (index_try .le. nym2)) then
	index_j=index_try
else
	index_j=index_try-1
endif
!!----------------------------------------------------------------------------
binary_index=binary_index_zs
index_try=index_try_zs
do while(binary_index .ge. 1)
	binary_index=binary_index-1
	if (zp .ge. za(index_try))  then
		if (index_try+binary_values(binary_index) .le. nzm2) &
		index_try=index_try+binary_values(binary_index)
	else
		index_try=index_try-binary_values(binary_index)		
	endif
enddo

if (zp .ge. za(index_try) .and. (index_try .le. nzm2)) then
	index_k=index_try
else
	index_k=index_try-1
endif

end subroutine vp_index0


subroutine vp_index_uni_stretch(vp, vpBound, index_i, index_j, index_k)
use trace_common
implicit none
real:: vp(0:2), vpBound(0:2), xp, yp, zp
integer:: index_i, index_j, index_k
!----------------------------------------------------------------------------
xp=vp(0)
yp=vp(1)
zp=vp(2)
!----------------------------------------------------------------------------
if ( .not. (xp .gt. xmin)) then
	xp=xmin; index_i=0
else if (xp .ge. xmax) then
	xp=xmax; index_i=nxm2
else 
	index_i=floor((xp-xmin)/dxa_uni)
	index_i=minval([index_i, nxm2])
endif
!----------------------------------------------------------------------------
if ( .not. (yp .gt. ymin)) then
	yp=ymin; index_j=0
else if (yp .ge. ymax) then
	yp=ymax; index_j=nym2
else 
	index_j=floor((yp-ymin)/dya_uni)
	index_j=minval([index_j, nym2])
endif
!----------------------------------------------------------------------------
if ( .not. (zp .gt. zmin)) then
	zp=zmin; index_k=0
else if (zp .ge. zmax) then
	zp=zmax; index_k=nzm2
else 
	index_k=floor((zp-zmin)/dza_uni)
	index_k=minval([index_k, nzm2])
endif
!----------------------------------------------------------------------------
vpBound=[xp, yp, zp]

end subroutine vp_index_uni_stretch


subroutine round_weight0(vp, round, weight)
use trace_common
implicit none
real:: w(0:1,0:2), weight(0:1,0:1,0:1), vp(0:2)
integer:: round(0:1,0:2), i, j, k
!----------------------------------------------------------------------------
do i=0,2
	if ( .not. (vp(i) .gt. 0.0)) then
! compared with vp(i) .le. 0.0, this way can avoid the crash from vp(i) .eq. NaN (by B=0)
		round(0,i)=0
		w(1,i)=0.0
	else if (vp(i) .ge. pmax(i)) then
		round(0,i)=r0max(i)
		w(1,i)=1.0
	else
		round(0,i)=floor(vp(i))
		w(1,i)=vp(i)-round(0,i)
	endif
enddo
 
round(1,:)=round(0,:)+1
w(0,:)=1.0-w(1,:)
forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)

end subroutine round_weight0


subroutine round_weight_stretch(vp, round, weight)
use trace_common
implicit none
real:: w(0:1,0:2), weight(0:1,0:1,0:1), vp(0:2), vpBound(0:2)
integer:: round(0:1,0:2), i, j, k, index_i, index_j, index_k
!----------------------------------------------------------------------------
call vp_index(vp, vpBound, index_i, index_j, index_k)

round(0,:)=[index_i, index_j, index_k]
round(1,:)=round(0,:)+1

w(0,0)=(xa(index_i+1)-vpBound(0))/dxa(index_i)
w(0,1)=(ya(index_j+1)-vpBound(1))/dya(index_j)
w(0,2)=(za(index_k+1)-vpBound(2))/dza(index_k)

w(1,:)=1.0-w(0,:)

forall(i=0:1,j=0:1,k=0:1) weight(i,j,k)=w(i,0)*w(j,1)*w(k,2)

end subroutine round_weight_stretch


subroutine curlB_grid0(i, j, k, CurlBp)
use trace_common
use field_common
implicit none
integer:: i, j, k
real:: gradBp(0:2,0:2), CurlBp(0:2)
!----------------------------------------------------------------------------
if (i .eq. 0) then
	gradBp(0,1:2) =-1.5*Bfield(1:2,0,j,k)   +2.0*Bfield(1:2,1,j,k)     -0.5*Bfield(1:2,2,j,k)
else  if (i .eq. nxm1) then
	gradBp(0,1:2) = 1.5*Bfield(1:2,nxm1,j,k)-2.0*Bfield(1:2,nxm1-1,j,k)+0.5*Bfield(1:2,nxm1-2,j,k)
else
	gradBp(0,1:2) = (Bfield(1:2,i+1,j,k)-Bfield(1:2,i-1,j,k))*0.5
endif
!----------------------------------------------------------------------------
if (j .eq. 0) then
	gradBp(1,0:2:2) =-1.5*Bfield(0:2:2,i,0,k)   +2.0*Bfield(0:2:2,i,1,k)     -0.5*Bfield(0:2:2,i,2,k)
else if (j .eq. nym1) then
	gradBp(1,0:2:2) = 1.5*Bfield(0:2:2,i,nym1,k)-2.0*Bfield(0:2:2,i,nym1-1,k)+0.5*Bfield(0:2:2,i,nym1-2,k)
else
	gradBp(1,0:2:2) = (Bfield(0:2:2,i,j+1,k)-Bfield(0:2:2,i,j-1,k))*0.5
endif
!----------------------------------------------------------------------------
if (k .eq. 0) then
	gradBp(2,0:1) =-1.5*Bfield(0:1,i,j,0)   +2.0*Bfield(0:1,i,j,1)     -0.5*Bfield(0:1,i,j,2)
else if (k .eq. nzm1) then
	gradBp(2,0:1) = 1.5*Bfield(0:1,i,j,nzm1)-2.0*Bfield(0:1,i,j,nzm1-1)+0.5*Bfield(0:1,i,j,nzm1-2)
else
	gradBp(2,0:1) = (Bfield(0:1,i,j,k+1)-Bfield(0:1,i,j,k-1))*0.5
endif
!----------------------------------------------------------------------------
curlBp=[gradBp(1,2)-gradBp(2,1), gradBp(2,0)-gradBp(0,2), gradBp(0,1)-gradBp(1,0)]

END subroutine curlB_grid0


subroutine curlB_grid_stretch(i, j, k, CurlBp)
use field_common
implicit none
integer:: i, j, k, ci0, cj0, ck0, s
real:: gradBp(0:2,0:2), coef(0:2, 0:2), CurlBp(0:2)
!----------------------------------------------------------------------------
call diff_coefficent(i, j, k, ci0, cj0, ck0, coef)
forall(s=1:2)   gradBp(0,s)=sum(coef(:, 0)*Bfield(s, i+ci0:i+ci0+2, j, k))
forall(s=0:2:2) gradBp(1,s)=sum(coef(:, 1)*Bfield(s, i, j+cj0:j+cj0+2, k))
forall(s=0:1)   gradBp(2,s)=sum(coef(:, 2)*Bfield(s, i, j, k+ck0:k+ck0+2))
 curlBp(0)=gradBp(1,2)-gradBp(2,1)
 curlBp(1)=gradBp(2,0)-gradBp(0,2)
 curlBp(2)=gradBp(0,1)-gradBp(1,0)
END subroutine curlB_grid_stretch


subroutine diff_coefficent(i, j, k, ci0, cj0, ck0, coef)
use trace_common
implicit none
integer:: i, j, k, ci0, cj0, ck0
real:: coef(0:2, 0:2)
!----------------------------------------------------------------------------
if (i .eq. 0) then
 	coef(0,0)=-(2.0*dxa(0)+dxa(1))/(dxa(0)*(dxa(0)+dxa(1)))
 	coef(1,0)=(dxa(0)+dxa(1))/(dxa(0)*dxa(1))
 	coef(2,0)=-dxa(0)/(dxa(1)*(dxa(0)+dxa(1)))
 	ci0=0
else if (i .eq. nxm1) then 	
 	coef(0,0)=  dxa(nxm1-1)/(dxa(nxm1-2)*(dxa(nxm1-1)+dxa(nxm1-2)))
 	coef(1,0)=-(dxa(nxm1-1)+dxa(nxm1-2))/(dxa(nxm1-1)*dxa(nxm1-2))
 	coef(2,0)= (2.0*dxa(nxm1-1)+dxa(nxm1-2))/(dxa(nxm1-1)*(dxa(nxm1-1)+dxa(nxm1-2)))
 	ci0=-2
else
	coef(0,0)=-dxa(i)/(dxa(i-1)*(dxa(i)+dxa(i-1)))
	coef(1,0)=(dxa(i)-dxa(i-1))/(dxa(i)*dxa(i-1))
	coef(2,0)= dxa(i-1)/(dxa(i)*(dxa(i)+dxa(i-1)))
	ci0=-1
endif
!----------------------------------------------------------------------------
if (j .eq. 0) then
 	coef(0,1)=-(2.0*dya(0)+dya(1))/(dya(0)*(dya(0)+dya(1)))
 	coef(1,1)=(dya(0)+dya(1))/(dya(0)*dya(1))
 	coef(2,1)=-dya(0)/(dya(1)*(dya(0)+dya(1)))
 	cj0=0
else if (j .eq. nym1) then
 	coef(0,1)=  dya(nym1-1)/(dya(nym1-2)*(dya(nym1-1)+dya(nym1-2)))
 	coef(1,1)=-(dya(nym1-1)+dya(nym1-2))/(dya(nym1-1)*dya(nym1-2))
 	coef(2,1)= (2.0*dya(nym1-1)+dya(nym1-2))/(dya(nym1-1)*(dya(nym1-1)+dya(nym1-2)))
 	cj0=-2
else
	coef(0,1)=-dya(j)/ (dya(j-1)*(dya(j)+dya(j-1)))
	coef(1,1)=(dya(j)-dya(j-1))/(dya(j)*dya(j-1))
	coef(2,1)= dya(j-1)/ (dya(j)*(dya(j)+dya(j-1)))
	cj0=-1
endif
!----------------------------------------------------------------------------
if (k .eq. 0) then
 	coef(0,2)=-(2*dza(0)+dza(1))/(dza(0)*(dza(0)+dza(1)))
 	coef(1,2)=(dza(0)+dza(1))/(dza(0)*dza(1))
 	coef(2,2)=-dza(0)/(dza(1)*(dza(0)+dza(1)))
 	ck0=0
else if (k .eq. nzm1) then 	
 	coef(0,2)=  dza(nzm1-1)/(dza(nzm1-2)*(dza(nzm1-1)+dza(nzm1-2))) 
 	coef(1,2)=-(dza(nzm1-1)+dza(nzm1-2))/(dza(nzm1-1)*dza(nzm1-2))
 	coef(2,2)= (2*dza(nzm1-1)+dza(nzm1-2))/(dza(nzm1-1)*(dza(nzm1-1)+dza(nzm1-2)))
 	ck0=-2
else
	coef(0,2)=-dza(k)/ (dza(k-1)*(dza(k)+dza(k-1)))
	coef(1,2)=(dza(k)-dza(k-1))/(dza(k)*dza(k-1))
	coef(2,2)= dza(k-1)/ (dza(k)*(dza(k)+dza(k-1)))
	ck0=-1
endif

end subroutine diff_coefficent


! trilinear interpolation
subroutine interpolateB(vp, bp)
use field_common
use trace_common
implicit none
real:: vp(0:2), bp(0:2), weight(0:1,0:1,0:1)
integer:: round(0:1,0:2), i
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
end subroutine interpolateB


subroutine interpolate_unit_vec_B(vp, unit_vec_bp)
!Unit vector B
use field_common
use trace_common
implicit none
real:: vp(0:2), bp(0:2), unit_vec_bp(0:2), weight(0:1,0:1,0:1)
integer:: round(0:1,0:2), i
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
unit_vec_bp= bp/norm2(bp)
end subroutine interpolate_unit_vec_B


subroutine interpolateAlpha(vp, alpha)
use field_common
use trace_common
implicit none
real:: weight(0:1,0:1,0:1), vp(0:2), bp(0:2), CurlBp(0:2), alpha
integer:: round(0:1,0:2), i
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
forall(i=0:2) CurlBp(i)=sum(weight*curlB(i, round(:,0), round(:,1), round(:,2)))
alpha=dot_product(curlbp, bp)/dot_product(bp, bp)
end subroutine interpolateAlpha


subroutine interpolate_k1(vp, unit_vec_bp, alpha, alphaFlag, ds_factor)
use trace_common
use field_common
implicit none
real:: weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), CurlBp(0:2), alpha, ds_factor
integer:: round(0:1,0:2), i
logical:: alphaFlag
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
if(alphaFlag) then
	forall(i=0:2) CurlBp(i)=sum(weight*curlB(i, round(:,0), round(:,1), round(:,2)))
	alpha=dot_product(curlbp, bp)/dot_product(bp, bp)
endif
unit_vec_bp= bp/norm2(bp)

if (stretchFlag) then
	ds_factor=norm2(unit_vec_bp/[dxa(round(0,0)),dya(round(0,1)),dza(round(0,2))])
else
	ds_factor=1.0
endif
end subroutine interpolate_k1


!d vp/dt=Bp/Bt
subroutine RK4(dt, vp0, vp1, alpha, alphaFlag)
implicit none
real:: dt, ds, alpha, vp0(0:2), vp1(0:2), k1(0:2), k2(0:2), k3(0:2), k4(0:2), ds_factor
logical:: alphaFlag
!----------------------------------------------------------------------------
! the unit of ds is same as the physical unit (e.g. the unit of xreg, yreg, zreg)
! the unit of dt is the scale of the cell (a self-adaptive fashion that varying from cell to cell)
call interpolate_k1(vp0, k1, alpha, alphaFlag, ds_factor)
ds=dt/ds_factor
call interpolate_unit_vec_B(vp0+ds*1./3.*k1,       k2)
call interpolate_unit_vec_B(vp0+ds*(-1./3.*k1+k2), k3)
call interpolate_unit_vec_B(vp0+ds*(k1-k2+k3),     k4)

vp1=vp0+ds/8.0*(k1+3.*(k2+k3)+k4)

end subroutine RK4


!d vp/d r_i=Bp/B_i 
subroutine RK4_Boundary(ds, vp0, vp1, b_dim)
implicit none
real ::  ds, vp0(0:2), vp1(0:2), bp(0:2), k1(0:2), k2(0:2), k3(0:2), k4(0:2)
integer:: b_dim
!----------------------------------------------------------------------------
call interpolateB(vp0, bp)
k1=bp/bp(b_dim)
k1(b_dim)=1.0
call interpolateB(vp0+ds*1./3.*k1, bp)
k2=bp/bp(b_dim)
k2(b_dim)=1.0
call interpolateB(vp0+ds*(-1./3.*k1+k2), bp)
k3=bp/bp(b_dim)
k3(b_dim)=1.0
call interpolateB(vp0+ds*(k1-k2+k3), bp)
k4=bp/bp(b_dim)
k4(b_dim)=1.0

vp1=vp0+ds/8.0*(k1+3.*(k2+k3)+k4)

end subroutine RK4_Boundary


subroutine RKF45(dt, vp0, vp1, alpha, alphaFlag, tol_this)
use trace_common
implicit none
real:: k1(0:2), k2(0:2), k3(0:2), k4(0:2), k5(0:2), k6(0:2), vp0(0:2), vp1(0:2), dvp(0:2), &
dt, ds0, ds1, ds, alpha, error, min_error, tol_this, tol_this_1, scale_dt, ds_factor
logical:: continue_flag, alphaFlag
integer:: rb, rb_index
!----------------------------------------------------------------------------
continue_flag=.true.
call interpolate_k1(vp0, k1, alpha, alphaFlag, ds_factor)
tol_this_1=tol_this/ds_factor

do while ( continue_flag ) 
	ds=dt/ds_factor
	call interpolate_unit_vec_B(vp0+ds*a21*k1,                                   k2)
	call interpolate_unit_vec_B(vp0+ds*(a31*k1+ a32*k2),                         k3)   
	call interpolate_unit_vec_B(vp0+ds*(a41*k1+ a42*k2+ a43*k3),                 k4)   
	call interpolate_unit_vec_B(vp0+ds*(a51*k1+ a52*k2+ a53*k3+ a54*k4),         k5)   
	call interpolate_unit_vec_B(vp0+ds*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5), k6)

	vp1 = vp0 + ds*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
	dvp = ds*(ce1*k1+ce3*k3+ce4*k4+ce5*k5+ce6*k6)
	error = norm2(dvp)

!----------------------------------------------------------------------------
	call vp_rboundary(vp1, rb, rb_index)
	if (rb .eq. 0) then
 		continue_flag = error .gt. tol_this_1 .and. (abs(dt) .gt. min_step)
		if (error .gt. 0.) then
			dt=dt* ((tol_this_1/error)**0.2)*0.9
		else
			dt=sign(100., dt)
		endif
		if (abs(dt) .gt. 100.) dt=sign(100., dt)
		if (abs(dt) .lt. min_step) dt=sign(min_step, dt)
!----------------------------------------------------------------------------
	else
		continue_flag =abs(dt) .gt. min_step
		if (continue_flag) then

			select case (rb)
			case(1:6)
				if(mod(rb, 2) .eq. 1) then
					ds0=pmin(rb_index)- vp0(rb_index)
					ds1= vp1(rb_index)-pmin(rb_index)
				else
					ds0=pmax(rb_index)- vp0(rb_index)
					ds1= vp1(rb_index)-pmax(rb_index)
				endif
				if (abs(ds0+ds1) .le. 0.05*norm2(vp1-vp0)) then
					 dt=sign(min_step, dt)
				else if (abs(dt*ds0/(ds0+ds1)) .lt. min_step) then
					continue_flag = .false.
				else
					! then if a next do loop exist, continue_flag will be .false. in that loop
					! because dt will \approx 0.5*sign(min_step, dt) in that loop
					dt=sign(abs(dt*ds0/(ds0+ds1))-0.5*min_step, dt)
				endif
			case(7)
				dt=dt*0.5
			end select
			! once dt is scaled to sign(min_step, dt), only once left to be processed in the do loop
			if (abs(dt) .lt. min_step) dt=sign(min_step, dt)
		endif
!----------------------------------------------------------------------------
	endif
enddo

end subroutine RKF45


subroutine vp_rboundary(vp, rb, rb_index)
use trace_common
implicit none
real:: vp(0:2)
integer:: k, rb, rb_index
logical:: boundary_mark(1:6)
!----------------------------------------------------------------------------
forall(k=0:2)
	boundary_mark(5-2*k)=vp(k)<pmin(k)
	boundary_mark(6-2*k)=vp(k)>pmax(k)
endforall

select case(count(boundary_mark))
	case(0)
		rb=0  !inside
	case(1)
		do k=1,6
			if (boundary_mark(k)) then
				rb=k
				exit
			endif
		enddo
		! Findloc is introducted from Fortran 2008 and later
		! rb=Findloc(boundary_mark, .true., 1)
		
		rb_index=(6-rb)/2
	case default
		rb=7
end select

end subroutine vp_rboundary


subroutine correct_foot(vp, vp1, sign_dt, rb)
use trace_common
implicit none
real:: dt, ds0, ds1, alpha, &
vp(0:2), vp0(0:2), vp1(0:2), vp_orig(0:2), vp1_orig(0:2), vp_tmp(0:2)
integer:: sign_dt, rb, rb_index, it
!---------------------------------------------------------------------------
call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 0) return

if (any((vp .eq. pmin) .or. (vp .eq. pmax))) then
	vp1=vp
	return
endif

vp0=vp
vp1_orig=vp1
vp_orig =vp
if (rb .ne. 7) then 
	if( mod(rb, 2) .eq. 1) then
		ds0=pmin(rb_index)- vp0(rb_index)
		ds1= vp1(rb_index)-pmin(rb_index)
	else
		ds0=pmax(rb_index)- vp0(rb_index)
		ds1= vp1(rb_index)-pmax(rb_index)
	endif

	if  (RK4flag) then
		dt=    step*abs(ds0/(ds0+ds1))*sign_dt*0.95
	else
		dt=min_step*abs(ds0/(ds0+ds1))*sign_dt*0.95
	endif
	
	if (abs(dt) .ge. min_step_foot) then
		call RK4(dt, vp0, vp, alpha, .false.)
		do while( .not.( all(pmin<=vp .and. vp<=pmax)) .and. (abs(dt) .ge. min_step_foot) )
			dt= dt*0.9
			if (abs(dt) .ge. min_step_foot) call RK4(dt, vp0, vp, alpha, .false.)
		enddo
	endif
endif

dt=min_step_foot*sign_dt
vp_tmp=vp
it=0

do while( all(pmin<=vp_tmp .and. vp_tmp<=pmax))
	call RK4(dt, vp_tmp, vp1, alpha, .false.)
	it=it+1	
	if (it .ge. maxsteps_foot) then
		vp0=vp_orig
		vp1=vp1_orig
		exit
	endif
	vp0=vp_tmp
	vp_tmp=vp1
end do

call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 7) return
	
if( mod(rb, 2) .eq. 1) then
	ds0=pmin(rb_index)- vp0(rb_index)
	ds1= vp1(rb_index)-pmin(rb_index)
else
	ds0=pmax(rb_index)- vp0(rb_index)
 	ds1= vp1(rb_index)-pmax(rb_index)
endif

if (abs(ds0+ds1) .le. 0.05*norm2(vp0-vp1)) then 	
	if((ds0+ds1) .ne. 0.0) then 
		vp1=(vp0*ds1+vp1*ds0)/(ds0+ds1)
	else
		vp1=vp0
	endif
else
	call RK4_Boundary(ds0, vp0, vp1, rb_index)
endif

end subroutine correct_foot


subroutine trace_bline(vp0, rs, re, rbs, rbe, length0, twist0, twistFlag, incline)
!+
! NAME :
!   trace_bline
!   
!PURPOSE:
!     Uses RK4 or RKF45 solver to integrate a line of force for a 3D magnetic field

!INPUTS:
!     vp0        - the position to start the line at
!     twistFlag  - compute the twist or not
!     incline    - abs(bn/norm2(bp)) at the cross section, for adjusting step or tol

!OUTPUTS:
!     rs,  re    - pixel coordinates of the start/end point of field lines
!     rbs, rbe   - categorize field lines based on where they thread the boundary of the field cube
!     length0    - length of the field line
!     twist0     - twist number of the field line
!-
use trace_common
implicit none
real:: vp0(0:2), vp1(0:2), vp(0:2), vp_tmp(0:2), rs(0:2), re(0:2), bp(0:2), dt, dL, dL0, length0, &
dtwist, twist0, alpha, alpha0, tol_this, tol_this1, step_this, incline, incline_this
integer:: it, sign_dt, rb, rbs, rbe
logical:: twistFlag, z0Flag
!----------------------------------------------------------------------------
twist0 =0.
length0=0.
!----------------------------------------------------------------------------
if (.not. all(pmin<=vp0 .and. vp0<=pmax)) then
	rbs=7
 	rbe=7
  	return
endif
!----------------------------------------------------------------------------

z0flag= vp0(2) .eq. zmin

if (incline .le. min_incline)  then
	incline_this = min_incline
else
	incline_this = incline
endif

if (RK4flag) then 
	step_this=step*incline_this
	if (step_this .le. min_step) step_this=min_step
else
	step_this=min_step
	tol_this=tol*(incline_this**1.5)
endif

do sign_dt=-1,1,2

	if (z0flag) then		
		call interpolateB(vp0, bp)		
		if (bp(2)*sign_dt .le. 0.) then
			if (sign_dt .eq. -1) then
				rs=vp0; rbs=1
			else
				re=vp0; rbe=1
			endif		
			cycle
		endif
	endif
	
	vp=vp0
	it=0
	dt=step_this*sign_dt
	dL=0.
	do while (all(pmin<=vp .and. vp<=pmax) .and. it < maxsteps)
		length0=length0+dL
		
		if (RK4flag) then  
		 	call RK4  (dt, vp, vp1, alpha, twistflag) 
		 	! alpha @ vp
		else		
			call RKF45(dt, vp, vp1, alpha, twistflag, tol_this)
		endif
		
		dL0=dL		
		dL=norm2(vp1-vp)
				
		if (twistflag) then
			if (it .ne. 0) then
				dtwist=(alpha0+alpha)/2.*dL0
				twist0=twist0+dtwist
			endif			
			alpha0=alpha
		endif
		
		it=it+1	
		vp_tmp=vp
		vp=vp1		
	end do
	
	call correct_foot(vp_tmp, vp1, sign_dt, rb)
	
	if (rb .eq. 0 .or. rb .eq. 7) then 
  		rbs=rb; rbe=rb
		return
	endif

	dL=norm2(vp1-vp_tmp)
	
	length0=length0+dL
	if (twistflag) then 
		call interpolateAlpha(vp1, alpha)
		dtwist=(alpha0+alpha)/2.*dL
		twist0=twist0+dtwist
	endif
	
	if (sign_dt .eq. -1) then
		rs=vp1; rbs=rb
	else
		re=vp1; rbe=rb
	endif
enddo

if (twistflag) twist0=twist0/(4.0*pi)

END subroutine trace_bline


subroutine read_ax()
use trace_common
implicit none
integer:: i, j, k, binary_index_top
!----------------------------------------------------------------------------
allocate(xa(0:nxm1))
open(unit=8,file='xa.bin', access='stream', status='old')
read(8) xa
 close(8)
	 
allocate(ya(0:nym1))
OPEN(unit=8,file='ya.bin', access='stream', status='old')
read(8) ya
 close(8)
	 
allocate(za(0:nzm1))
OPEN(unit=8,file='za.bin', access='stream', status='old')
read(8) za
 close(8)

allocate(dxa(0:nxm2))
allocate(dya(0:nym2))
allocate(dza(0:nzm2))

forall(i=0:nxm2) dxa(i)=xa(i+1)-xa(i)
forall(j=0:nym2) dya(j)=ya(j+1)-ya(j)
forall(k=0:nzm2) dza(k)=za(k+1)-za(k)

xmin=xa(0); xmax=xa(nxm1)
ymin=ya(0); ymax=ya(nym1)
zmin=za(0); zmax=za(nzm1)
!----------------------------------------------------------------------------
uni_stretch_Flag=(minval(abs(dxa))/maxval(abs(dxa)) .gt. 0.99) .and. &
                 (minval(abs(dya))/maxval(abs(dya)) .gt. 0.99) .and. &
                 (minval(abs(dza))/maxval(abs(dza)) .gt. 0.99)
                                  
if (uni_stretch_Flag) then
	dxa_uni=(xmax-xmin)/nxm1
	dya_uni=(ymax-ymin)/nym1
	dza_uni=(zmax-zmin)/nzm1
else
!for binary search
	binary_index_xs=floor(dlog10(dble(nxm2))/dlog10(2.0D0))
	binary_index_ys=floor(dlog10(dble(nym2))/dlog10(2.0D0))
	binary_index_zs=floor(dlog10(dble(nzm2))/dlog10(2.0D0))
	binary_index_top=maxval([binary_index_xs, binary_index_ys, binary_index_zs])
	allocate(binary_values(0:binary_index_top))
	forall(i=0:binary_index_top) binary_values(i)=2**i
	index_try_xs=binary_values(binary_index_xs)
	index_try_ys=binary_values(binary_index_ys)
	index_try_zs=binary_values(binary_index_zs)
endif

end subroutine read_ax


subroutine logical2int(logical_in, interger_out)
implicit none
logical:: logical_in
integer:: interger_out
if (logical_in) then
	interger_out=1
else
	interger_out=0
endif
end subroutine logical2int


subroutine cross_product(v1, v2, v3)
implicit none
real:: v1(0:2), v2(0:2), v3(0:2)
v3(0)=dble(v1(1))*v2(2)-dble(v1(2))*v2(1)
v3(1)=dble(v1(2))*v2(0)-dble(v1(0))*v2(2)
v3(2)=dble(v1(0))*v2(1)-dble(v1(1))*v2(0)
end subroutine cross_product


subroutine initialize()
use qfactor_common
use trace_common
use field_common
implicit none
external vp_index0, vp_index_uni_stretch, round_weight0, round_weight_stretch, &
grad_unit_vec_B_grid0, grad_unit_vec_B_grid_stretch, curlB_grid0, curlB_grid_stretch
!----------------------------------------------------------------------------
integer:: i, j, k, s, &
twistFlag_int, RK4flag_int, scottFlag_int, csFlag_int, curlB_out_int, &
q0flag_int, vflag_int, cflag_int, reclen, nbridges, nx_mag, ny_mag
real, allocatable:: Bfield_tmp(:, :, :, :), magnetogram(:, :)
real:: point1(0:2), point2(0:2), vp(0:2), bp(0:2), delta_mag
logical:: xa_exist, ya_exist, za_exist, ifort_flag, curlB_out
!----------------------------------------------------------------------------
open(unit=8, file='head.txt', status='old')
read(8, *) nx, ny, nz, nbridges, delta, maxsteps, &
           xreg, yreg, zreg, step, tol, &
           twistFlag_int, RK4flag_int, scottFlag_int, csFlag_int, curlB_out_int
 close(8)
 
 twistFlag=twistFlag_int .eq. 1
   RK4Flag=  RK4flag_int .eq. 1
 scottFlag=scottFlag_int .eq. 1
    csFlag=   csFlag_int .eq. 1
 curlB_out=curlB_out_int .eq. 1
!----------------------------------------------------------------------------
! for RKF45
if (.not. RK4flag) then
	a21=   1./4.
	a31=   3./32.;   a32=    9./32.
	a41=1932./2197.; a42=-7200./2197.; a43=  7296./2197.
	a51= 439./216.;  a52=-8.;          a53=  3680./513.;   a54= -845./4104.
	a61=  -8./27.;   a62= 2.;          a63= -3544./2565.;  a64= 1859./4104.; a65=-11./40.
	b1 =  16./135.;   b3= 6656./12825.; b4= 28561./56430.;  b5=   -9./50.;    b6=  2./55. 
	ce1=   1./360.;  ce3= -128./4275.; ce4= -2197./75240.; ce5=    1./50.;   ce6=  2./55.
endif
!----------------------------------------------------------------------------
nxm1=nx-1; nym1=ny-1; nzm1=nz-1
nxm2=nx-2; nym2=ny-2; nzm2=nz-2
NaN =transfer(2143289344, 1.0)
!----------------------------------------------------------------------------
! read Bx, By, Bz
allocate(Bfield(0:2, 0:nxm1, 0:nym1, 0:nzm1))
allocate(Bfield_tmp( 0:nxm1, 0:nym1, 0:nzm1, 0:2))
open(unit=8, file='b3d.bin', access='stream', status='old')
read(8) Bfield_tmp
 close(8)
 
! if reclen .eq. 4, compiled by gfortran; if reclen .eq. 1, compiled by ifort
inquire(iolength=reclen) 1.0 
ifort_flag=(reclen .eq. 1)

if (ifort_flag) then 
	CALL OMP_set_num_threads(nbridges)
else
	! the setting of 4 threads achieves the best performance here for gfortran
	CALL OMP_set_num_threads(minval([4,nbridges])) 
endif

!switch the order of indexes for a better efficiency
!$OMP PARALLEL DO  PRIVATE(k), schedule(DYNAMIC)
do k=0, nzm1
	forall(s=0:2) Bfield(s,:,:,k)=Bfield_tmp(:,:,k,s)
enddo
!$OMP END PARALLEL DO

deallocate(Bfield_tmp)

if (reclen .eq. 4) CALL OMP_set_num_threads(nbridges)
!----------------------------------------------------------------------------
inquire(file='xa.bin', exist=xa_exist)
inquire(file='ya.bin', exist=ya_exist)
inquire(file='za.bin', exist=za_exist)
stretchFlag= xa_exist .and. ya_exist .and. za_exist

if (stretchFlag) then
	call read_ax()
	if (uni_stretch_Flag) then
		vp_index=>vp_index_uni_stretch
	else	
		vp_index=>vp_index0
	endif
	round_weight        =>round_weight_stretch
	curlB_grid          =>curlB_grid_stretch
	grad_unit_vec_B_grid=>grad_unit_vec_B_grid_stretch
else
	xmin=0.0; xmax=nxm1
	ymin=0.0; ymax=nym1
	zmin=0.0; zmax=nzm1
	
	round_weight        =>round_weight0
	curlB_grid          =>curlB_grid0
	grad_unit_vec_B_grid=>grad_unit_vec_B_grid0
endif
r0max=[nxm2, nym2, nzm2]
pmin=[xmin,ymin,zmin]
pmax=[xmax,ymax,zmax]
!----------------------------------------------------------------------------
if (twistFlag .or. curlB_out) then
	allocate(curlB(0:2, 0:nxm1, 0:nym1, 0:nzm1))

	!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC) 
	do k=0, nzm1
	do j=0, nym1
	do i=0, nxm1
		call curlB_grid(i, j, k, curlB(:,i,j,k))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
	
	if (curlB_out) then
		open(unit=8, file='curlB.bin', access='stream', status='replace')
		write(8) curlB
		close(8)
		deallocate(curlB, Bfield)
		if (stretchFlag) then
			deallocate(xa, ya, za, dxa, dya, dza)
			if (.not. uni_stretch_Flag) deallocate(binary_values)
		endif
		call exit
	endif
endif
!----------------------------------------------------------------------------
if (stretchFlag) then
	delta_mag=minval([dxa,dya])

	nx_mag=(xmax-xmin)/delta_mag+1
	ny_mag=(ymax-ymin)/delta_mag+1

	allocate(magnetogram(0:nx_mag-1, 0:ny_mag-1))
	do j=0, ny_mag-1
	do i=0, nx_mag-1
		vp=[xmin+delta_mag*i, ymin+delta_mag*j, zmin]
		call interpolateB(vp, bp)
		magnetogram(i, j)=bp(2)		
	enddo
	enddo

	open(unit=8, file='magnetogram.bin', access='stream', status='replace')
	write(8) magnetogram
	close(8)
	 
	deallocate(magnetogram)

	open(8,file='mag_info.txt', status='replace')	
	write(8,*) nx_mag, ny_mag, delta_mag
	close(8)
endif
!----------------------------------------------------------------------------
q0Flag=(zreg(0) .eq. zmin) .and. (zreg(1) .eq. zmin) .and. (.not. csflag)
 vFlag=(xreg(1) .ne. xreg(0)) .and. (yreg(1) .ne. yreg(0)) .and. (zreg(1) .ne. zreg(0)) .and. (.not. csflag)
 cFlag=(.not. vflag) .and. (.not. q0flag)

if (csflag) then
	Normal_index =-1
	point0=[xreg(0),yreg(0),zreg(0)]
	point1=[xreg(1),yreg(1),zreg(0)]
	point2=[xreg(0),yreg(0),zreg(1)]
	
	q1=norm2(point1-point0)/delta+1
	q2=norm2(point2-point0)/delta+1
	qx=0; qy=0; qz=0
	!ev*: elementary vector of the cut plane
	ev1=(point1-point0)/norm2(point1-point0)
	ev2=(point2-point0)/norm2(point2-point0)
	call cross_product(ev1, ev2, ev3)
	ev3=ev3/norm2(ev3)
else
	qx=(xreg(1)-xreg(0))/delta+1
	qy=(yreg(1)-yreg(0))/delta+1
	qz=(zreg(1)-zreg(0))/delta+1

	if (vflag) then
		q1=qx; q2=qy; Normal_index=2
	endif
	if (qx .eq. 1) then 
		q1=qy; q2=qz; Normal_index=0; cut_coordinate=xreg(0)
	endif	
	if (qy .eq. 1) then 
		q1=qx; q2=qz; Normal_index=1; cut_coordinate=yreg(0)
	endif
	if (qz .eq. 1) then 
		q1=qx; q2=qy; Normal_index=2; cut_coordinate=zreg(0)
	endif
endif

q2m1=q2-1; q1m1=q1-1

call logical2int(q0flag, q0flag_int)
call logical2int( cflag,  cflag_int)
call logical2int( vflag,  vflag_int)

! tell IDL these informations
open(unit=8, file='tail.txt', status='replace')
write(8, *) qx, qy, qz, q1, q2
write(8, *) q0flag_int, cflag_int, vflag_int
 close(8)
!----------------------------------------------------------------------------
 !maxsteps    =nint(4*(nx+ny+nz)/step)*10
min_incline  =0.05
if (stretchFlag) then
	min_step=minval([step,delta/maxval([dxa,dya,dza])])
else
	min_step=minval([step,delta])
endif

min_step_foot=min_step/2.0
if (RK4flag) then
	maxsteps_foot=    step/min_step_foot*4
else
	maxsteps_foot=min_step/min_step_foot*4
endif
!----------------------------------------------------------------------------
grad3DFlag= scottFlag .or. vflag .or. stretchflag
if (grad3DFlag) then
	allocate(grad_unit_vec_Bfield(0:2, 0:2, 0:nxm1, 0:nym1, 0:nzm1))
	!$OMP PARALLEL DO  PRIVATE(i, j, k), schedule(DYNAMIC) 
	do k=0, nzm1
	do j=0, nym1
	do i=0, nxm1
		call grad_unit_vec_B_grid(i, j, k, grad_unit_vec_Bfield(:,:,i,j,k))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
endif
!----------------------------------------------------------------------------

allocate(     q(0:q1m1, 0:q2m1))
allocate(   bnr(0:q1m1, 0:q2m1))
allocate(length(0:q1m1, 0:q2m1))
allocate(reF(0:2, 0:q1m1, 0:q2m1))
allocate(reboundary(-2:q1+1, -2:q2+1))
allocate(rboundary_tmp(0:q1m1, 0:q2m1))

if (twistFlag) allocate( twist(0:q1m1, 0:q2m1))
if (scottFlag) allocate(q_perp(0:q1m1, 0:q2m1))

if (vflag .or. cflag) then
	allocate(rsboundary(-2:q1+1, -2:q2+1))
	allocate(rsF(0:2, 0:q1m1, 0:q2m1))
	allocate(tangent_Flag(-1:q1, -1:q2))
endif

end subroutine initialize
