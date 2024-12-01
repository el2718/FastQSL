subroutine grad_unit_vec_B_grid0(i, j, k, grad_unit_vec_B)
use trace_common
use field_common
implicit none
integer:: i, j, k
real:: grad_unit_vec_B(0:2, 0:2)
!----------------------------------------------------------------------------
if (i .eq. 0) then
	grad_unit_vec_B(0,:) = &
	-1.5*Bfield(:,0,j,k)/norm2(Bfield(:,0,j,k)) &
	+2.0*Bfield(:,1,j,k)/norm2(Bfield(:,1,j,k)) &
	-0.5*Bfield(:,2,j,k)/norm2(Bfield(:,2,j,k))
else if (i .eq. nxm1) then	
	grad_unit_vec_B(0,:) = &
	 1.5*Bfield(:,i  ,j,k)/norm2(Bfield(:,i  ,j,k)) &
	-2.0*Bfield(:,i-1,j,k)/norm2(Bfield(:,i-1,j,k)) &
	+0.5*Bfield(:,i-2,j,k)/norm2(Bfield(:,i-2,j,k))
else
	grad_unit_vec_B(0,:) = &
	(Bfield(:,i+1,j,k)/Norm2(Bfield(:,i+1,j,k))-Bfield(:,i-1,j,k)/Norm2(Bfield(:,i-1,j,k)))*0.5
endif
!----------------------------------------------------------------------------
if (j .eq. 0) then
	grad_unit_vec_B(1,:) = &
	-1.5*Bfield(:,i,0,k)/norm2(Bfield(:,i,0,k)) &
	+2.0*Bfield(:,i,1,k)/norm2(Bfield(:,i,1,k)) &
	-0.5*Bfield(:,i,2,k)/norm2(Bfield(:,i,2,k))
else if (j .eq. nym1) then	
	grad_unit_vec_B(1,:) = &
	 1.5*Bfield(:,i,j  ,k)/norm2(Bfield(:,i,j  ,k)) &
	-2.0*Bfield(:,i,j-1,k)/norm2(Bfield(:,i,j-1,k)) &
	+0.5*Bfield(:,i,j-2,k)/norm2(Bfield(:,i,j-2,k))	
else
	grad_unit_vec_B(1,:) = &
	(Bfield(:,i,j+1,k)/Norm2(Bfield(:,i,j+1,k))-Bfield(:,i,j-1,k)/Norm2(Bfield(:,i,j-1,k)))*0.5
endif
!----------------------------------------------------------------------------
if (k .eq. 0) then
	grad_unit_vec_B(2,:) = &
	-1.5*Bfield(:,i,j,0)/norm2(Bfield(:,i,j,0)) &
	+2.0*Bfield(:,i,j,1)/norm2(Bfield(:,i,j,1)) &
	-0.5*Bfield(:,i,j,2)/norm2(Bfield(:,i,j,2))
else if (k .eq. nzm1) then
	grad_unit_vec_B(2,:) = &
	 1.5*Bfield(:,i,j,k  )/norm2(Bfield(:,i,j,k  )) &
	-2.0*Bfield(:,i,j,k-1)/norm2(Bfield(:,i,j,k-1)) &
	+0.5*Bfield(:,i,j,k-2)/norm2(Bfield(:,i,j,k-2))
else
	grad_unit_vec_B(2,:) = &
	(Bfield(:,i,j,k+1)/Norm2(Bfield(:,i,j,k+1))-Bfield(:,i,j,k-1)/Norm2(Bfield(:,i,j,k-1)))*0.5
endif

END subroutine grad_unit_vec_B_grid0


subroutine grad_unit_vec_B_grid_stretch(i, j, k, grad_unit_vec_B)
use field_common
implicit none
integer:: i, j, k, ci0, cj0, ck0, t
real:: coef(0:2, 0:2), grad_unit_vec_B(0:2,0:2)
!----------------------------------------------------------------------------
call diff_coefficent(i, j, k, ci0, cj0, ck0, coef)
grad_unit_vec_B=0.0
do t=0,2
	grad_unit_vec_B(0,:)=grad_unit_vec_B(0,:)+&
		coef(t, 0)*Bfield(:, i+ci0+t, j, k)/norm2(Bfield(:, i+ci0+t, j, k))
enddo
do t=0,2
	grad_unit_vec_B(1,:)=grad_unit_vec_B(1,:)+&
		coef(t, 1)*Bfield(:, i, j+cj0+t, k)/norm2(Bfield(:, i, j+cj0+t, k))
enddo
do t=0,2
	grad_unit_vec_B(2,:)=grad_unit_vec_B(2,:)+&
		coef(t, 2)*Bfield(:, i, j, k+ck0+t)/norm2(Bfield(:, i, j, k+ck0+t))		
enddo

END subroutine grad_unit_vec_B_grid_stretch


subroutine interpolate_grad_unit_vec_B(vp, unit_vec_bp, grad_unit_vec_Bp)
use trace_common
use field_common
implicit none
real:: weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), &
grad_unit_vec_Bp(0:2,0:2), grad_unit_vec_B_cell(0:2,0:2,0:1,0:1,0:1)
integer:: round(0:1,0:2), i, j, k
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
unit_vec_bp= bp/norm2(bp)

if (grad3DFlag) then
	forall(i=0:2, j=0:2) grad_unit_vec_Bp(i,j)=sum(weight*grad_unit_vec_Bfield(i, j, round(:,0), round(:,1), round(:,2)))
else
!the output is identical as the upper, but don't require grad_unit_vec_Bfield
!the efficiency is 1/4.79(gfortran) or 1/2.21(ifort) times of the upper
	do k=0,1
	do j=0,1
	do i=0,1
		if (weight(i,j,k) .ne. 0.0) then
			call grad_unit_vec_B_grid(round(i,0), round(j,1), round(k,2), grad_unit_vec_B_cell(:,:,i,j,k))
		else
			!avoid NaN
			grad_unit_vec_B_cell(:,:,i,j,k)=0.0
		endif
	enddo
	enddo
	enddo
	forall(i=0:2, j=0:2) grad_unit_vec_Bp(i,j)=sum(weight*grad_unit_vec_B_cell(i, j, :, :, :))
endif

end subroutine interpolate_grad_unit_vec_B


subroutine interpolateAll(vp, unit_vec_bp, grad_unit_vec_Bp, alpha, alphaFlag, ds_factor)
use trace_common
use field_common
implicit none
real:: weight(0:1,0:1,0:1), vp(0:2), bp(0:2), unit_vec_bp(0:2), grad_unit_vec_Bp(0:2,0:2), &
grad_unit_vec_B_cell(0:2,0:2,0:1,0:1,0:1), alpha, ds_factor, CurlBp(0:2)
integer:: round(0:1,0:2), i, j, k
logical:: alphaFlag
!----------------------------------------------------------------------------
call round_weight(vp, round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
unit_vec_bp= bp/norm2(bp)

if (grad3DFlag) then
	forall(i=0:2, j=0:2) grad_unit_vec_Bp(i,j)=sum(weight*grad_unit_vec_Bfield(i, j, round(:,0), round(:,1), round(:,2)))
else
	do k=0,1
	do j=0,1
	do i=0,1
		if (weight(i,j,k) .ne. 0.0) then
			call grad_unit_vec_B_grid(round(i,0), round(j,1), round(k,2), grad_unit_vec_B_cell(:,:,i,j,k))
		else
			grad_unit_vec_B_cell(:,:,i,j,k)=0.0
		endif
	enddo
	enddo
	enddo
	forall(i=0:2, j=0:2) grad_unit_vec_Bp(i,j)=sum(weight*grad_unit_vec_B_cell(i, j, :, :, :))
endif

if(alphaFlag) then
	forall(i=0:2) CurlBp(i)=sum(weight*curlB(i, round(:,0), round(:,1), round(:,2)))
	alpha=dot_product(curlbp, bp)/dot_product(bp, bp)
endif

if (stretchFlag) then
	ds_factor=norm2(unit_vec_bp/[dxa(round(0,0)),dya(round(0,1)),dza(round(0,2))])
else
	ds_factor=1.0
endif

end subroutine interpolateAll


subroutine f_scott(vector9, vector9_k)
implicit none
real:: vector9(0:8), vector9_k(0:8), grad_unit_vec_Bp(0:2,0:2)
integer:: i
!----------------------------------------------------------------------------
call interpolate_grad_unit_vec_B(vector9(0:2), vector9_k(0:2), grad_unit_vec_Bp)
forall(i=0:2) vector9_k(3+i)= dot_product(vector9(3:5), grad_unit_vec_Bp(0:2,i))
forall(i=0:2) vector9_k(6+i)= dot_product(vector9(6:8), grad_unit_vec_Bp(0:2,i))
end subroutine f_scott


subroutine f_scott_k1(vector9, vector9_k, alpha, alphaFlag, ds_factor)
implicit none
real:: vector9(0:8), vector9_k(0:8), grad_unit_vec_Bp(0:2,0:2), alpha, ds_factor
integer:: i
logical:: alphaFlag
!----------------------------------------------------------------------------
call interpolateAll(vector9(0:2), vector9_k(0:2), grad_unit_vec_Bp, alpha, alphaFlag, ds_factor)
forall(i=0:2) vector9_k(3+i)= dot_product(vector9(3:5), grad_unit_vec_Bp(0:2,i))
forall(i=0:2) vector9_k(6+i)= dot_product(vector9(6:8), grad_unit_vec_Bp(0:2,i))
end subroutine f_scott_k1


subroutine f_scott_boundary(vector9, vector9_k, b_dim)
implicit none
real:: vector9(0:8), vector9_k(0:8)
integer:: b_dim
!----------------------------------------------------------------------------
call f_scott(vector9, vector9_k)
vector9_k=vector9_k/vector9_k(b_dim)
vector9_k(b_dim)=1.0
end subroutine f_scott_boundary


subroutine RK4_scott(dt, vector9, vector9_1, alpha, alphaFlag)
implicit none
real:: vector9(0:8), vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8), &
dt, ds, alpha, ds_factor
logical:: alphaFlag
!----------------------------------------------------------------------------
call f_scott_k1(vector9, k1, alpha, alphaFlag, ds_factor)
ds=dt/ds_factor
call f_scott(vector9+ds*1./3.*k1,       k2)
call f_scott(vector9+ds*(-1./3.*k1+k2), k3)
call f_scott(vector9+ds*(k1-k2+k3),     k4)

vector9_1=vector9+ds/8.0*(k1+3.*k2+3.*k3+k4)
end subroutine RK4_scott


subroutine RK4_scott_boundary(ds, vector9, vector9_1, b_dim)
implicit none
real:: ds, vector9(0:8),vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8)
integer:: b_dim
!----------------------------------------------------------------------------
call f_scott_boundary(vector9,                   k1, b_dim)
call f_scott_boundary(vector9+ds*1./3.*k1,       k2, b_dim)
call f_scott_boundary(vector9+ds*(-1./3.*k1+k2), k3, b_dim)
call f_scott_boundary(vector9+ds*(k1-k2+k3),     k4, b_dim)
vector9_1=vector9+ds/8.0*(k1+3.*k2+3.*k3+k4)
end subroutine RK4_scott_boundary


subroutine RKF45_scott(dt, vector9, vector9_1, alpha, alphaFlag)
use trace_common
implicit none
real:: vector9(0:8), vector9_1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8), k5(0:8), k6(0:8), &
dt, ds, ds0, ds1, error, min_error, ds_factor, vp0(0:2), vp1(0:2), dvp(0:2), scale_dt, alpha, tol_this_1
logical:: continue_flag, alphaFlag
integer:: rb, rb_index
!----------------------------------------------------------------------------
vp0=vector9(0:2)
continue_flag=.true.
call f_scott_k1(vector9, k1, alpha, alphaFlag, ds_factor)
tol_this_1=tol/ds_factor

do while ( continue_flag ) 
	ds=dt/ds_factor
	call f_scott(vector9+ds*a21*k1,                               k2)
	call f_scott(vector9+ds*(a31*k1+a32*k2),                      k3)
	call f_scott(vector9+ds*(a41*k1+a42*k2+a43*k3),               k4)   
	call f_scott(vector9+ds*(a51*k1+a52*k2+a53*k3+a54*k4),        k5)   
	call f_scott(vector9+ds*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5), k6)

	vector9_1 = vector9+ds*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
	
	dvp=ds*(ce1*k1(0:2)+ce3*k3(0:2)+ce4*k4(0:2)+ce5*k5(0:2)+ce6*k6(0:2))
	error = norm2(dvp)
	
!----------------------------------------------------------------------------
	vp1=vector9_1(0:2)
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

end subroutine RKF45_scott


subroutine correct_foot_scott(vector9, vector9_1, sign_dt, rb)
use trace_common
implicit none
real:: dt, dt0, ds0, ds1, vp(0:2), vp0(0:2), vp1(0:2), alpha, &
vector9(0:8), vector9_0(0:8), vector9_1(0:8), vector9_orig(0:8), vector9_1_orig(0:8)
integer:: sign_dt, rb, rb_index, it
!---------------------------------------------------------------------------
vp=vector9(0:2)
vp1=vector9_1(0:2)
dt0=dt

call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 0) return

if (any((vp .eq. pmin) .or. (vp .eq. pmax))) then
	vector9_1=vector9
	return
endif

vector9_orig  =vector9
vector9_1_orig=vector9_1
vector9_0=vector9

if (rb .ne. 7) then 
	vp0=vector9_0(0:2)
	vp1=vector9_1(0:2)
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
		call RK4_scott(dt, vector9_0, vector9, alpha, .false.)
		vp=vector9(0:2)
		do while(.not.(all(pmin<=vp .and. vp<=pmax)) .and. (abs(dt) .ge. min_step_foot))
			dt= dt*0.9
			if (abs(dt) .ge. min_step_foot) call RK4_scott(dt, vector9_0, vector9, alpha, .false.)
			vp=vector9(0:2)
		enddo
	endif
endif

dt=min_step_foot*sign_dt
vp=vector9(0:2)
it=0

do while(all(pmin<=vp .and. vp<=pmax))
	call RK4_scott(dt, vector9, vector9_1, alpha, .false.)
	it=it+1
	if (it .ge. maxsteps_foot) then
		vector9_0=vector9_orig
		vector9_1=vector9_1_orig
		dt=dt0
		exit
	endif	
	vector9_0=vector9
	vector9  =vector9_1
	vp       =vector9(0:2)
end do

vp =vector9_0(0:2)
vp1=vector9_1(0:2)

call vp_rboundary(vp1, rb, rb_index)
if (rb .eq. 7) return

if(mod(rb, 2) .eq. 1) then
	 ds0=pmin(rb_index)-  vp(rb_index)
	 ds1= vp1(rb_index)-pmin(rb_index)
else
	 ds0=pmax(rb_index)-  vp(rb_index)
	 ds1= vp1(rb_index)-pmax(rb_index)
endif

if (abs(ds0+ds1) .le. 0.05*norm2(vp0-vp1)) then
 	call RK4_scott(dt*ds0/(ds0+ds1), vector9_0, vector9_1, alpha, .false.)
else
	call RK4_scott_Boundary(ds0, vector9_0, vector9_1, rb_index)
endif

end subroutine correct_foot_scott


!Scott_2017_ApJ_848_117
subroutine trace_scott(vp0, q0, q_perp0, rs, re, rbs, rbe, length0, twist0, twistFlag)
use trace_common
implicit none
real:: vp(0:2), vp0(0:2), q0, q_perp0, dt, dL, dL0, length0, alpha, alpha0, twist0, dtwist, &
rs(0:2), re(0:2), Bn_s, Bn_e, b0_square, b0(0:2), bs(0:2), be(0:2), &
vector9(0:8), vector9_0(0:8), vector9_1(0:8), vector9_s(0:8), vector9_e(0:8), &
u0(0:2), us(0:2), ue(0:2), v0(0:2), vs(0:2), ve(0:2), vs1(0:2), ve1(0:2), us1(0:2), ue1(0:2)
integer:: it, sign_dt, rb, rbe, rbs, e_index, s_index, maxdim, index1, index2
logical:: z0flag, twistFlag
!----------------------------------------------------------------------------
twist0 =0.0
length0=0.0
!----------------------------------------------------------------------------
z0flag= vp0(2) .eq. zmin
call interpolateB(vp0, b0)

b0_square=dot_product(b0, b0)

maxdim=sum(maxloc(abs(b0)))-1
index1=mod(maxdim+1,3)
index2=mod(maxdim+2,3)

v0(maxdim)= b0(index1)
v0(index1)=-b0(maxdim)
v0(index2)=0.

v0=v0/norm2(v0)
call cross_product(b0, v0, u0)
u0=u0/norm2(u0)
!----------------------------------------------------------------------------
do sign_dt=-1, 1, 2	
	vector9(0:2)=vp0
	vector9(3:5)=u0
	vector9(6:8)=v0
	
	if (z0flag) then
		if ( b0(2)*sign_dt .le. 0.) then 
			if (sign_dt .eq. -1) then
				vector9_s=vector9; rbs=1
			else
				vector9_e=vector9; rbe=1
			endif		
			cycle
		endif
	endif	
	
	it=0
	dL=0.
	vp=vp0
	if (RK4flag) then
		dt=    step*sign_dt
	else
		dt=min_step*sign_dt
	endif	
	
	do while(all(pmin<=vp .and. vp<=pmax) .and. abs(it) < maxsteps)

		length0=length0+dL		
		
		if (RK4flag) then  
		 	call   RK4_scott(dt, vector9, vector9_1, alpha, twistflag)
		else			
			call RKF45_scott(dt, vector9, vector9_1, alpha, twistflag)
		endif
		
		dL0=dL
		dL =norm2(vector9_1(0:2)-vector9(0:2))
		
		if (twistflag) then
			if (it .ne. 0) then
				dtwist=(alpha0+alpha)/2.*dL0
				twist0=twist0+dtwist
			endif			
			alpha0=alpha
		endif
				
		it       =it+sign_dt	
		vector9_0=vector9
		vector9  =vector9_1
		vp       =vector9_1(0:2)
	end do
	
	call correct_foot_scott(vector9_0, vector9_1, sign_dt, rb)
	
	if (rb .eq. 0 .or. rb .eq. 7) then 
		rbs=rb; rbe=rb
		q0=NaN; q_perp0=NaN
		return
	endif

	dL=norm2(vector9_1(0:2)-vector9_0(0:2))
	length0=length0+dL
	
	if (twistflag) then
		vp=vector9_1(0:2)
		call interpolateAlpha(vp, alpha)
		dtwist=(alpha0+alpha)/2.*dL
		twist0=twist0+dtwist
	endif
	
	if (sign_dt .eq. -1) then
		vector9_s=vector9_1; rbs=rb
	else
		vector9_e=vector9_1; rbe=rb
	endif
enddo

if (twistflag) twist0=twist0/(4.0*pi)
!----------------------------------------------------------------------------
rs=vector9_s(0:2)
call interpolateB(rs, bs)
s_index=(6-rbs)/2
Bn_s=bs(s_index)
us=vector9_s(3:5)
vs=vector9_s(6:8)

re=vector9_e(0:2)
call interpolateB(re, be)
e_index=(6-rbe)/2
Bn_e=be(e_index)
ue=vector9_e(3:5)
ve=vector9_e(6:8)

us1=us-us(s_index)/Bn_s*bs
vs1=vs-vs(s_index)/Bn_s*bs
ue1=ue-ue(e_index)/Bn_e*be
ve1=ve-ve(e_index)/Bn_e*be

q0      =abs(dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
        +    dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
        -2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
      (b0_square / abs(Bn_s*Bn_e))

ue1=ue-dot_product(ue,be)/norm2(be)*(be/norm2(be))
ve1=ve-dot_product(ve,be)/norm2(be)*(be/norm2(be))
us1=us-dot_product(us,bs)/norm2(bs)*(bs/norm2(bs))
vs1=vs-dot_product(vs,bs)/norm2(bs)*(bs/norm2(bs))

q_perp0 =abs(dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
        +    dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
        -2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
            (b0_square / (norm2(bs)*norm2(be)))
     
end subroutine trace_scott
