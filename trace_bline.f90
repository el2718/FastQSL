module trace_common
use field_common
implicit none
	integer:: maxsteps, maxsteps_foot, normal_index
	real:: step, min_step, min_step_foot, tol, min_incline, NaN, ev3(0:2)
	logical:: RK4flag, twistFlag
!----------------------------------------------------------------------------
! for RKF45
	real:: b1, b3, b4, b5, b6, ce1, ce3, ce4, ce5, ce6, &
	a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65
contains

!d vp/dt=Bp/Bt
subroutine RK4(dt, vector0, vector1, scott_this, alpha, alphaFlag)
implicit none
real:: dt, ds, vector0(0:8), vector1(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8), ds_factor, alpha
logical:: scott_this, alphaFlag
!----------------------------------------------------------------------------
! the unit of ds is same as the physical unit (e.g. the unit of xreg, yreg, zreg)
! the unit of dt is the scale of the cell (a self-adaptive fashion that varying from cell to cell)
call interpolate_dvds(vector0, k1, scott_this, ds_factor, alpha, alphaFlag)
ds=dt/ds_factor
call interpolate_dvds(vector0+ds*1./3.*k1,       k2, scott_this)
call interpolate_dvds(vector0+ds*(-1./3.*k1+k2), k3, scott_this)
call interpolate_dvds(vector0+ds*(k1-k2+k3),     k4, scott_this)

vector1=vector0+ds/8.0*(k1+3.*(k2+k3)+k4)
end subroutine RK4


!d vp/d r_i=Bp/B_i
subroutine RK4_Boundary(ds, vector0, vector1, scott_this, b_dim)
implicit none
real:: ds, vector0(0:8), vector1(0:8), dvds(0:8), k1(0:8), k2(0:8), k3(0:8), k4(0:8)
integer:: b_dim
logical:: scott_this
!----------------------------------------------------------------------------
call interpolate_dvds(vector0,                   dvds, scott_this)
k1=dvds/dvds(b_dim)
call interpolate_dvds(vector0+ds*1./3.*k1,       dvds, scott_this)
k2=dvds/dvds(b_dim)
call interpolate_dvds(vector0+ds*(-1./3.*k1+k2), dvds, scott_this)
k3=dvds/dvds(b_dim)
call interpolate_dvds(vector0+ds*(k1-k2+k3),     dvds, scott_this)
k4=dvds/dvds(b_dim)

vector1=vector0+ds/8.0*(k1+3.*(k2+k3)+k4)
end subroutine RK4_Boundary


subroutine RKF45(dt, vector0, vector1, scott_this, alpha, alphaFlag, tol_this)
implicit none
real:: vector0(0:8), vector1(0:8), vp0(0:2), vp1(0:2), vp2(0:2), &
k1(0:8), k2(0:8), k3(0:8), k4(0:8), k5(0:8), k6(0:8), dvp(0:2), &
dt, dt_try, ds, ds_factor, ds0, ds1, error, min_error, tol_this, tol_this_1, scale_dt, alpha
logical:: continue_flag, alphaFlag, scott_this
integer:: rb, rb_index
!----------------------------------------------------------------------------
call interpolate_dvds(vector0, k1, scott_this, ds_factor, alpha, alphaFlag)

tol_this_1=tol_this/ds_factor
vp0=vector0(0:2)

continue_flag=.true.
do while (continue_flag)
	ds=dt/ds_factor
	call interpolate_dvds(vector0+ds* a21*k1,                                  k2, scott_this)
	call interpolate_dvds(vector0+ds*(a31*k1+ a32*k2),                         k3, scott_this)
	call interpolate_dvds(vector0+ds*(a41*k1+ a42*k2+ a43*k3),                 k4, scott_this)   
	call interpolate_dvds(vector0+ds*(a51*k1+ a52*k2+ a53*k3+ a54*k4),         k5, scott_this)   
	call interpolate_dvds(vector0+ds*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5), k6, scott_this)

	vector1 = vector0 + ds*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
	vp1     = vector1(0:2)
	dvp = ds*(ce1*k1(0:2)+ce3*k3(0:2)+ce4*k4(0:2)+ce5*k5(0:2)+ce6*k6(0:2))

	!dvp is very small, distance(vp1, vp1-dvp) not accurate enough
	if (sphericalFlag) then
		error= norm2([dvp(0)*vp1(2)*cos(vp1(1)),dvp(0)*vp1(2),dvp(2)])
	else
		error= norm2(dvp)
	endif

	if (abs(dt) .gt. min_step) then
		call trim_size(vp0, vp1, rb, rb_index, ds0, ds1)
		select case (rb)
		case(0)
			continue_flag = error .gt. tol_this_1
			if (continue_flag) dt=sign(maxval([abs(dt*0.618), min_step]), dt)
		case(1:6)

			if (abs(ds) .ge. 2.*distance(vp1, vp0)) then
				dt_try=dt/2.
			else if (abs(dt*ds0/(ds0+ds1)) .lt. min_step) then
				dt_try=0. ! to make continue_flag = .false.
			else
				! then if a next do loop exist, continue_flag will be .false. in that loop
				! because dt will \approx 0.5*sign(min_step, dt) in that loop
				dt_try=dt*abs(ds0/(ds0+ds1))-0.5*sign(min_step, dt)
			endif

			continue_flag = abs(dt_try) .ge. min_step
			if (continue_flag) dt = dt_try
		case(7)
			dt=sign(maxval([abs(dt*0.618), min_step]), dt)
			continue_flag = .true.
		end select
	else
	! once dt is scaled to sign(min_step, dt), only once left to be processed in the do loop
		continue_flag = .false.
	endif
enddo

! 0.618^5 \approx 0.090145133
min_error=(tol_this_1*0.090145133)/((100./abs(dt))**5.)
if (error .lt. min_error) then
	if (error .ne. 0.0) dt=sign(100., dt)
else
	scale_dt=((tol_this_1/error)**0.2)*0.618
	dt=dt*scale_dt
endif

dt=sign(maxval([abs(dt), min_step]), dt)

end subroutine RKF45


subroutine correct_foot(vector, vector_1, dt, rb, scott_this)
implicit none
real:: dt, ds0, ds1, vp(0:2), vp1(0:2), bp(0:2), alpha, &
vector(0:8), vector_0(0:8), vector_1(0:8), vector_orig(0:8), vector_1_orig(0:8)
integer:: rb, rb_index, it
logical:: scott_this
!---------------------------------------------------------------------------
vp=vector(0:2)
vp1=vector_1(0:2)

call trim_size(vp, vp1, rb, rb_index, ds0, ds1)

if (rb .eq. 0) then
	return
else
	if ((period0flag .and. any((vp(1:2) .eq. pmin(1:2)) .or. (vp(1:2) .eq. pmax(1:2)))) &
	     .or. ((.not. period0flag) .and. any((vp .eq. pmin) .or. (vp .eq. pmax)))) then
		vector_1=vector
		return
	endif
endif

vector_orig  =vector
vector_1_orig=vector_1

if (rb .ne. 7) then
	if (abs(dt*ds0/(ds0+ds1)) .ge. 1.5*min_step_foot) then
	! if the residual distance to the boundary is smaller 1.5*min_step_foot, 
	! then only not more than 2 steps is needed in the following do while loop

		dt=dt*abs(ds0/(ds0+ds1))-0.5*sign(min_step_foot, dt)
		call RK4(dt, vector_orig, vector, scott_this, alpha, .false.)
		! the output vector should be inside
		! the residual distance to the boundary \approx 0.5*min_step_foot
		if (.not. inside(vector(0:2))) vector=vector_orig
	endif
endif

dt=sign(min_step_foot, dt)
it=0
do while(inside(vector(0:2)) .and. (it .le. maxsteps_foot))
	call RK4(dt, vector, vector_1, scott_this, alpha, .false.)
	it=it+1
	vector_0=vector
	vector  =vector_1
end do

if (inside(vector_1(0:2))) then
	vector_0=vector_orig
	vector_1=vector_1_orig
endif

vp =vector_0(0:2)
vp1=vector_1(0:2)

call trim_size(vp, vp1, rb, rb_index, ds0, ds1)
if (rb .eq. 7) return

call interpolateB(vp, bp)
if (abs(bp(rb_index)) .le. 0.05*norm2(bp)) then
	vector_1 = (vector_0*ds1 + vector_1*ds0)/(ds0+ds1)
else
	call RK4_Boundary(ds0, vector_0, vector_1, scott_this, rb_index)
endif

end subroutine correct_foot


! vp is inside, vp1 is outside
subroutine trim_size(vp, vp1, rb, rb_index, ds0, ds1)
real:: ds0, ds1, vp(0:2), vp1(0:2), vp_tmp(0:2)
integer:: k, rb, rb_index
logical:: boundary_mark(1:6)
!----------------------------------------------------------------------------
vp_tmp=vp1
if (sphericalFlag) vp_tmp(0)=modulo(vp1(0), two_pi)

forall(k=0:2)
	boundary_mark(5-2*k) = .not. (vp_tmp(k) .ge. pmin(k))
	boundary_mark(6-2*k) = .not. (vp_tmp(k) .le. pmax(k))
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
		rb_index=(6-rb)/2
	case default
		rb=7
end select
!----------------------------------------------------------------------------
if ((rb .ne. 0) .and. (rb .ne. 7)) then
	if(mod(rb, 2) .eq. 1) then
		ds0=pmin(rb_index)-  vp(rb_index)
		ds1= vp1(rb_index)-pmin(rb_index)
	else
		ds0=pmax(rb_index)-  vp(rb_index)
		ds1= vp1(rb_index)-pmax(rb_index)
	endif
endif

end subroutine trim_size


subroutine calculate_incline(vp, bp, tangentFlag, bnp, incline)
implicit none
real:: vp(0:2), bp(0:2), bp_car(0:2), bnp, incline
logical:: tangentFlag
!----------------------------------------------------------------------------
if (Normal_index .ne. -1) then
	bnp=bp(normal_index)
else
	if (sphericalFlag) then
		bp_car(0)= bp(2)*cos(vp(1))*cos(vp(0))
		bp_car(1)= bp(2)*cos(vp(1))*sin(vp(0))
		bp_car(2)= bp(2)*sin(vp(1))
	else
		bp_car   = bp
	endif
	bnp=dot_product(bp_car, ev3)
endif

tangentFlag=incline .lt. min_incline

if (tangentFlag) then
	if ((vp(2) .eq. zmin) .and. normal_index .eq. 2) then
		incline=min_incline
	else
		bnp=norm2(bp)
		! set for the plane perpendicular to the field line
		incline=1.0
	endif
else
	incline=abs(bnp/norm2(bp))
endif

end subroutine calculate_incline


subroutine interpolate_dvds(vector, dvds, scott_this, ds_factor, alpha, alphaFlag)
implicit none
integer:: round(0:1,0:2), i, j, k
real:: vector(0:8), dvds(0:8), dbdcp(0:2,0:2), dbdc_cell(0:2,0:2,0:1,0:1,0:1), &
weight(0:1,0:1,0:1), bp(0:2), uni_vec_bp(0:2), r, sin_lat_rad, cos_lat_rad, CurlBp(0:2)
real, optional:: ds_factor, alpha
logical:: scott_this
logical, optional:: alphaFlag
!----------------------------------------------------------------------------
call round_weight(vector(0:2), round, weight)
forall(i=0:2) Bp(i)=sum(weight*Bfield(i, round(:,0), round(:,1), round(:,2)))
uni_vec_bp= bp/norm2(bp)

if (sphericalFlag) then
	r=vector(2)
	cos_lat_rad=cos(vector(1))
	dvds(0)=uni_vec_bp(0)/(r*cos_lat_rad)
	dvds(1)=uni_vec_bp(1)/ r
	dvds(2)=uni_vec_bp(2)
else
	dvds(0:2)=uni_vec_bp
endif

if(present(ds_factor))then
	if (stretchFlag) then
		ds_factor=norm2(dvds(0:2)/[dxa(round(0,0)), dya(round(0,1)), dza(round(0,2))])
	else
		ds_factor=1.0
	endif
	
	if(present(alphaFlag)) then
	if(alphaFlag) then
		forall(i=0:2) CurlBp(i)=sum(weight*curlB(i, round(:,0), round(:,1), round(:,2)))
		alpha=dot_product(CurlBp, bp)/dot_product(bp, bp)
	endif
	endif
endif

if (scott_this) then

	if (grad3DFlag) then
		forall(i=0:2, j=0:2) dbdcp(i,j)=sum(weight*dbdc_field(i, j, round(:,0), round(:,1), round(:,2)))
	else
	!this output is identical as the upper, it don't require dbdc_field while takes more time
		do k=0,1
		do j=0,1
		do i=0,1
			if (weight(i,j,k) .ne. 0.0) then
				call dbdc_grid(round(i,0), round(j,1), round(k,2), dbdc_cell(:,:,i,j,k))
			else
				!avoid NaN
				dbdc_cell(:,:,i,j,k)=0.0
			endif
		enddo
		enddo
		enddo
		forall(i=0:2, j=0:2) dbdcp(i,j)=sum(weight*dbdc_cell(i,j,:,:,:))
	endif
	!--------------------------------------------------------------------------
	if (sphericalFlag) then
		sin_lat_rad=sin(vector(1))

		dbdcp(:,0)=[dbdcp(0,0)-uni_vec_bp(1)*sin_lat_rad+uni_vec_bp(2)*cos_lat_rad, &
					sin_lat_rad*dvds(0)*r+dbdcp(1,0), -dvds(0)*cos_lat_rad+dbdcp(2,0)]
		dbdcp(:,1)=[dbdcp(0,1), dbdcp(1,1)+dvds(2), dbdcp(2,1)-dvds(1)]

		dbdcp(0,:)=dbdcp(0,:)/(r*cos_lat_rad)
		dbdcp(1,:)=dbdcp(1,:)/ r
	endif

	forall(i=0:2) dvds(3+i)= dot_product(vector(3:5), dbdcp(:,i))
	forall(i=0:2) dvds(6+i)= dot_product(vector(6:8), dbdcp(:,i))
endif
end subroutine interpolate_dvds


subroutine trace_bline(vp, rs, re, rbs, rbe, length0, twist0, bspe, inclineFlag, q0, q_perp0)
implicit none
integer:: it, sign_dt, rb, rbe, rbs, e_index, s_index, maxdim, index1, index2
real:: vp(0:2), dt, dL, dL0, alpha, alpha0, dtwist, rs(0:2), re(0:2), Bn_s, Bn_e, &
vector(0:8), vector_0(0:8), vector_1(0:8), vector_s(0:8), vector_e(0:8), &
u0(0:2), us(0:2), ue(0:2), v0(0:2), vs(0:2), ve(0:2), vs1(0:2), ve1(0:2), us1(0:2), ue1(0:2), &
bs(0:2), be(0:2), btmp(0:2), bp(0:2), bp_square, bnp, length0, incline, step_this, tol_this
logical:: alphaFlag, tangentFlag, scott_this
real, optional:: twist0, bspe(0:2, -1:1), q0, q_perp0
logical, optional:: inclineFlag
!----------------------------------------------------------------------------
if (.not. inside(vp)) then
	rbs=7; rs=NaN
	rbe=7; re=NaN
	length0=NaN
	if (present(twist0))  twist0 =NaN
	if (present(bspe))    bspe   =NaN
	return
endif
!----------------------------------------------------------------------------
if (present(twist0)) then
	alphaFlag=twistFlag
else
	alphaFlag=.false.
endif

if (alphaFlag) then
	twist0=0.
	alpha =0.
endif

length0=0.
!----------------------------------------------------------------------------
call interpolateB(vp, bp)

scott_this=present(q0)
if (scott_this) then
	bp_square=dot_product(bp, bp)

	maxdim=sum(maxloc(abs(bp)))-1
	index1=mod(maxdim+1,3)
	index2=mod(maxdim+2,3)

	v0(maxdim)= bp(index1)
	v0(index1)=-bp(maxdim)
	v0(index2)= 0.

	v0=v0/norm2(v0)
	call cross_product(bp, v0, u0)
	u0=u0/norm2(u0)
else
	u0=0.
	v0=0.
endif
!----------------------------------------------------------------------------
if (present(inclineFlag)) then
	call calculate_incline(vp, bp, tangentFlag, bnp, incline)

	if (RK4flag) then
		step_this=step*incline
		if (step_this .le. min_step) step_this=min_step
	else
		step_this=min_step
		tol_this=tol*(incline**1.5)
	endif
else
	tol_this=tol
	step_this=step
endif
!----------------------------------------------------------------------------
do sign_dt = -1, 1, 2

	vector=[vp, u0, v0]

	if ((vp(2) .eq. zmin) .and. (bp(2)*sign_dt .lt. 0.)) then
		if (sign_dt .eq. -1) then
			vector_s=vector; rbs=1; rs=vp; bs=bp
		else
			vector_e=vector; rbe=1; re=vp; be=bp
		endif
		cycle
	endif
	
	it=0
	dL=0.
	
	dt=step_this*sign_dt

	do while(inside(vector(0:2)) .and. abs(it) < maxsteps)
		length0=length0+dL
		
		if (RK4flag) then
		 	call RK4  (dt, vector, vector_1, scott_this, alpha, alphaFlag)
		else
		 	call RKF45(dt, vector, vector_1, scott_this, alpha, alphaFlag, tol_this)
		endif

		dL0=dL
		dL = distance(vector_1(0:2), vector(0:2))

		if (alphaFlag) then
			if (it .ne. 0) then
				dtwist=(alpha0+alpha)/2.*dL0
				twist0=twist0+dtwist
			endif
			alpha0=alpha
		endif

		it      =it+sign_dt
		vector_0=vector
		vector  =vector_1
	end do
	
	call correct_foot(vector_0, vector_1, dt, rb, scott_this)

	if (rb .eq. 0 .or. rb .eq. 7) then 
		if (sign_dt .eq. -1) then
			rs=vector_1(0:2); rbs=rb
		else
			re=vector_1(0:2); rbe=rb
		endif
		if (present(q0))      q0     =NaN
		if (present(q_perp0)) q_perp0=NaN
		return
	endif

	dL = distance(vector_1(0:2), vector_0(0:2))
	length0=length0+dL
	
	call interpolateB(vector_1(0:2), btmp, alpha, alphaFlag)
	
	if (alphaFlag) then
		dtwist=(alpha0+alpha)/2.*dL
		twist0=twist0+dtwist
	endif

	if (sign_dt .eq. -1) then
		vector_s=vector_1; rs=vector_s(0:2); rbs=rb; bs=btmp
	else
		vector_e=vector_1; re=vector_e(0:2); rbe=rb; be=btmp
	endif
enddo

if (present(bspe)) then
	bspe(:,-1)=bs
	bspe(:, 0)=bp
	bspe(:, 1)=be
endif

if (alphaFlag) twist0=twist0/(4.0*pi)
!----------------------------------------------------------------------------
if (.not. present(q0)) return

!Scott_2017_ApJ_848_117
s_index=(6-rbs)/2
Bn_s=bs(s_index)
us=vector_s(3:5)
vs=vector_s(6:8)

e_index=(6-rbe)/2
Bn_e=be(e_index)
ue=vector_e(3:5)
ve=vector_e(6:8)

us1=us-us(s_index)/Bn_s*bs
vs1=vs-vs(s_index)/Bn_s*bs
ue1=ue-ue(e_index)/Bn_e*be
ve1=ve-ve(e_index)/Bn_e*be

q0      =abs(dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
        +    dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
        -2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
      (bp_square / abs(Bn_s*Bn_e))

if (present(q_perp0)) then
ue1=ue-dot_product(ue,be)/norm2(be)*(be/norm2(be))
ve1=ve-dot_product(ve,be)/norm2(be)*(be/norm2(be))
us1=us-dot_product(us,bs)/norm2(bs)*(bs/norm2(bs))
vs1=vs-dot_product(vs,bs)/norm2(bs)*(bs/norm2(bs))

q_perp0 =abs(dot_product(Ue1,Ue1)*dot_product(Vs1,Vs1)  &
        +    dot_product(Us1,Us1)*dot_product(Ve1,Ve1)  &
        -2.0*dot_product(Ue1,Ve1)*dot_product(Us1,Vs1))/&
            (bp_square / (norm2(bs)*norm2(be)))
endif

end subroutine trace_bline

end module trace_common
