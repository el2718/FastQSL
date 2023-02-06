PRO qfactor, bx, by, bz, xa=xa, ya=ya, za=za, xreg=xreg, yreg=yreg, zreg=zreg, csFlag=csFlag, $
             factor=factor, delta=delta,  RK4Flag=RK4Flag, step=step, tol=tol, maxsteps=maxsteps, $
             scottFlag=scottFlag, twistFlag=twistFlag, curlB_out=curlB_out, odir=odir, fstr=fstr, $
             nbridges=nbridges, no_preview=no_preview, tmpB=tmpB, RAMtmp=RAMtmp
;+
; PURPOSE:
;   Calculate the squashing factor Q at the photosphere or a cross section or a box volume, 
;   given a 3D magnetic field with Cartesian, uniform or stretched grids
; 
;   For details see:
;   Zhang, P., Chen, J.*, Liu, R. and Wang, C., 2022, FastQSL: A Fast Computation Method for Quasi-separatrix Layers. The Astrophysical Journal, 937, 26
;
;   FORTRAN procedures used
;   qfactor.f90
;   trace_bline.f90
;   trace_scott.f90
;
; ------COMPILATION 
;   For Linux and MacOS (either by ifort or gfortran): 
;      ifort -o qfactor.x qfactor.f90 -fopenmp -O3 -xHost -ipo
;   gfortran -o qfactor.x qfactor.f90 -fopenmp -Ofast -march=native
;
;   please specify the path of qfactor.x at the line of "spawn, 'qfactor.x' " in this file, 
;   or move qfactor.x to the $PATH (e.g. /usr/local/bin/) of the system
;
;   For Windows: (executing "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" in cmd first would be necessary)
;      ifort /o qfactor.exe qfactor.f90 /Qopenmp /O3 /QxHost /Qipo
;   please replace qfactor.x by qfactor.exe for the path specifing
;
;   -O3, -xHost, -ipo, -Ofast, -march=native are for a better efficiency;
;   -Ofast would be problematic for MacOS, then please substitutue -O3 for it; 
;
; INPUTS
;   Bx, By, Bz: 3D magnetic field, will be forcibly converted to float arrays while writing 'b3d.bin'
;   
; OPTIONAL INPUTS  
;
;   xa, ya, za:  coordinates of grids in 1D arrays. The size of xa, ya, za must be consistant with the size of the 3D magnetic field, 
;                and the values in x{yz}a should increase by order
;   stretchFlag= keyword_set(xa) and keyword_set(ya) and keyword_set(za), which means bx, by, bz are on stretched grids, otherwise are on uniform grids; 
;                if invoked, the units of x{yz}reg, delta, length are same as x{yz}a, 
;                while the units of step and tol are the scale of a cell (a self-adaptive fashion that varying from cell to cell)
;
;   xreg, yreg, zreg: coordinates of the output region, in arrays of two elements, default is to include all available pixels at the photosphere.
;                 If stretchFlag=0B, the unit of coordinates is pixel, and the origin of the 3D field is [0, 0, 0].
;                 --- If (xreg[0] ne xreg[1]) AND (yreg[0] ne yreg[1]) AND (zreg[0] ne zreg[1]) AND NOT csFlag
;                     calculate Q in a box volume
;                     invoke vFlag
;                 --- If (xreg[0] eq xreg[1]) OR (yreg[0] eq yreg[1]) or (zreg[0] eq zreg[1])
;                     calculate Q in a cross section perpendicular to X or Y or Z axis
;                     invoke cFlag 
;                 --- If zreg=[0, 0], calculate Q at the photosphere
;                     invoke z0Flag 
;                 --- If csFlag is set, see below
;                     invoke cFlag
;   
;   csFlag:     to calculate Q in a cross section defined by three points; default is 0
;               point0=[xreg[0],yreg[0],zreg[0]] ; the origin of output
;               point1=[xreg[1],yreg[1],zreg[0]] ; point0 -> point1, first axis
;               point2=[xreg[0],yreg[0],zreg[1]] ; point0 -> point2, second axis
;               x{yz}reg[0] does not have to be smaller than x{yz}reg[1] in this case
;
;   factor:     to bloat up the original resolution, i.e. grid spacing of output = 1/factor; default is 4 
;
;   delta:      grid spacing of output; default is 1/factor (when stretchFlag=0B); 
;               if it is set, the keyword factor will be ignored
;
;   RK4Flag:    to trace bline by RK4; default is 0B (RKF45)
;   			   
;   step:       step size in tracing field lines for RK4; default=1.0
;
;   tol:        tolerance of a step in RKF45; default is 10.^(-4)
;
;   maxsteps:   maxium steps for stracing a field line at one direction; default is 4*(nx+ny+nz)/step; 
;               if highly twisted field lines exist and traced by RK4, this value should be larger; suggested by Jiang, Chaowei
;
;   scottFlag:  calculate Q and Q_perp by the method of Scott_2017_ApJ_848_117; default is 0B (method 3 of Pariat_2012_A&A_541_A78)
;		
;   twistFlag:  to calculate twist number Tw; see Liu_2016_ApJ_818_148; default is 0
;
;   curlB_out:  to save curlB at odir+'curlB.sav'; curlBx, curlBy, curlBz have same dimensions as Bx, By, Bz; default is 0
;
;   odir:       directory to save the results
;		
;   fstr:       filename of the results; file_sav=odir+fstr+'.sav'
;
;   nbridges:   number of processors to engage; default is !CPU.HW_NCPU - 2
; 
;   no_preview: don't produce PNG images for preview; default is 0B
;
;   RAMtmp:     use RAM (/dev/shm/tmp/) to speed up the data transmission, please only invoke it in Linux; default is 0B; 
;		if invoked, please run only one task of qfactor.pro simultaneously on one machine
;
;   tmpB:       apply temporary() to Bx, By, Bz to reduce the memory occupation of 3D magnetic field in IDL after writing b3d.bin to tmp_dir; default is 0B
;
;  memory occupation in qfactor.x: 
;       3D magnetic field + some 2D arrays + 3D curlB field (same as the occupation of the 3D magnetic field, when twistFlag=1B) +
;       grad_unit_vec_Bfield (3 times as the occupation of the 3D magnetic field)
;
;
; OUTPUTS:
;   q/qcs/q3d:  squashing factor 
;  
;   q_perp/q_perp3d:     q_perp in Titov 2007, ApJ, 660, 863
;
;   slogq/slogq_perp:    sign(Bz) x log_10(Q/Q_\perp); Calculated with all the field lines
;   slogq_orig/slogq_perp_orig:  only include field lines with both footpoints on the photoshere
;  
;   twist/twist3d:       twist number = \int \mu_0 J \cdot B /(4\pi) B^2 ds
;
;   rboundary/rsboundary/reboundary: (r:remote, s:start point, e:end point)
;                nature of the ends of field lines, see 'subroutine vp_rboundary' in trace_bline.f90
;                0 - inside
;                1 - end at zmin
;                2 - end at zmax
;                3 - end at ymin
;                4 - end at ymax
;                5 - end at xmin
;                6 - end at xmax
;                7 - others
;   rboundary3d=rsboundary3d+8*reboundary3d, for saving storage; So if rboundary3d[i, j, k] eq 9B, both two mapping surfaces of q3d[i, j, k] are the photosphere
;
;   rF/rsF/reF: coordinates of mapping points (F:foot); suggested by Jiang, Chaowei
; 
;   length: length of field lines
;
;   Bnr: abs(Bn_local/Bn_target) at the photosphere
;
; MODIFICATION HISTORY:
;   Developed by R. Liu, J. Chen and Peijing, Zhang @ USTC
;   
;   Jun 30,2014 R. Liu @ USTC, IDL edition
;   Apr 21,2015 R. Liu and J. Chen, deal with field lines pass through the boundary other than bottom
;   Apr 29,2015 R. Liu and J. Chen, further optimization on shared memory
;   Apr 27,2015 R. Liu, qcs
;   Jun 15,2015 J. Chen, Fortran Edition, correct foot points with RK4_Boundary
;   Jul  8,2015 J. Chen, q3d
;   Oct 29,2015 J. Chen, deal with field lines touching the cut plane: use the plane quasi-perp to the field line  
;   Nov 1, 2015 J. Chen,
;		(1) fuse qcs and qfactor in qfactor.f90;  
;		(2) the cross section can be paralleled to the photosphere;
;		(3) set tmp_dir;	  
;   Nov 4, 2015 J. Chen, introduce rboundary3d				
;   Jun 22,2016 J. Chen, add the map of field line length
;   Oct 30,2017 J. Chen, add the map of Bnr
;   Aug 28,2018 J. Chen, supplement Q at maginal points
;   May  1,2021 Peijing Zhang and J. Chen, trace field line with RKF45; RK4 is remainded, modify classic RK4 to the 3/8-rule RK4
;   May  1,2021 J. Chen, adapted to gfortran compiler
;   Jun  1,2021 J. Chen, forcibly convert the input Bx, By, Bz to float arrays in IDL (Real(4) in Fortran)
;   Jun 10,2021 J. Chen, provide a keyword of maxsteps for the case that has extremely long field lines
;   Jun 11,2021 J. Chen, add the coordinates of mapping points to '*.sav' data
;   Jun 13,2021 J. Chen, mark 0 for inside and 7 for others in 'subroutine vp_rboundary'
;   Jul  5,2021 J. Chen, switch the order of indexes of Bfield in trace_bline.f90 for a better efficiency
;   Jul  9,2021 J. Chen, provide the option with the method of Scott_2017_ApJ_848_117, and q_perp as an output
;   Dec  1,2021 J. Chen, optimize the correction of foot point
;   Dec 13,2021 J. Chen, adjust tol or step by incline
;   Dec 25,2021 J. Chen, remove the reliance of Solar SoftWare
;   Jan 26,2022 J. Chen, fix a bug of RKF45 in case the difference between RK4 and RK5 is 0
;   Jan 30,2022 J. Chen, remove doppler_color_mix, due to the poor recognizability of green-white-yellow
;   Feb 16,2022 J. Chen, reduce the memory occupation for 3D case in qfactor.x
;   Apr 27,2022 J. Chen, 
;		(1) fix a bug of "segmentation fault" during the output of qfactor.x;
;		(2) calculate slogq in qfactor.x, sign2d.pro is not more necessary;
;		(3) forcibly assign the minimum value of Q to 2, the theoretical minimum;
;		(4) zreg[0] can be non-zero for 3D case;
;		(5) add a keyword of tmpB;
; 		(6) cut_str can be a format of '(i0)', '(f0.1)' or '(f0.2)', according to the value of input
;   Apr 29,2022 J. Chen, extract subroutine round_weigh() for interpolation
;   May 10,2022 J. Chen, determine qx, qy, qz, z0flag, cflag, vflag in qfactor.x
;   May 19,2022 J. Chen, check the existence of nulls on grids; add a keyword of RAMtmp
;   Jun 10,2022 J. Chen, adapt to stretched grids
;   Oct 11,2022 J. Chen, adapt to Windows
;   Oct 13,2022 J. Chen, check the existence of infinite or NaN values on grids
;   Nov 25,2022 J. Chen, add a keyword of curlB_out to save curlB
;   Jan 17,2023 J. Chen, integrate doppler color in qfactor.pro, doppler_color.pro is not more necessary;
;                        to aviod an error for a remote server: % TVLCT: Unable to open X Windows display
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted
;   provided this disclaimer and information are included unchanged.
;-
;----------------------------------------------------------------------------------------------
; check input; region of output
sbx=size(Bx)
sby=size(By)
sbz=size(Bz)

if sbx[0] ne 3 or sby[0] ne 3 or sbz[0] ne 3 then message, 'Bx, By and Bz must be 3D arrays!'
if sbx[1] ne sby[1] or sbx[1] ne sbz[1] or $
   sbx[2] ne sby[2] or sbx[2] ne sbz[2] or $
   sbx[3] ne sby[3] or sbx[3] ne sbz[3] then message, 'Bx, By and Bz must have the same dimensions!' 
nx=sbz[1] & ny=sbz[2] & nz=sbz[3]

; check the existence of nulls on grids, to avoid "segmentation fault occurred"
ss=where((bx eq 0.0) and (by eq 0.0) and (bz eq 0.0))
n_ss=n_elements(ss)
if (n_ss gt 3) then message, string(n_ss)+" nulls (B=0) on grids! Too many!"

; check the existence of infinite or NaN values on grids, to avoid "segmentation fault occurred"
ss=where(~finite(bx) or ~finite(by) or ~finite(bz))
if (ss[0] ne -1) then message, "there are some infinite or NaN vaules on grids"


stretchFlag=keyword_set(xa) and keyword_set(ya) and keyword_set(za)

if (stretchFlag) then begin
	nxa=n_elements(xa) & nya=n_elements(ya) & nza=n_elements(za)
	if (nx ne nxa) or (ny ne nya) or (nz ne nza) then $
	message, 'the size of xa, ya and za must be consistant with the size of the 3D magnetic field'
	
	if ~keyword_set(xreg) then xreg=[xa(0), xa(nx-1)]
	if ~keyword_set(yreg) then yreg=[ya(0), ya(ny-1)]
	if ~keyword_set(zreg) then zreg=[za(0), za(0)]
endif else begin
	if ~keyword_set(xreg) then xreg=[0, nx-1]
	if ~keyword_set(yreg) then yreg=[0, ny-1]
	if ~keyword_set(zreg) then zreg=[0,    0]
	magnetogram=Bz[*,*,0]
endelse


if n_elements(xreg) ne 2 or n_elements(yreg) ne 2 or n_elements(zreg) ne 2 then $
message, 'xreg, yreg and zreg must be 2-element arrays!'

if total([xreg[0] eq xreg[1], yreg[0] eq yreg[1], zreg[0] eq zreg[1]]) ge 2 then $
message, 'Something is wrong with the cross section .......'

if keyword_set(csFlag) then csFlag=1B else csFlag=0B

if csFlag and (zreg[0] eq zreg[1]) then message, 'Something is wrong with the cross section .......'
;----------------------------------------------------------------------------------------------
max_threads=!CPU.HW_NCPU
if ~keyword_set(nbridges)   then nbridges=max_threads-2
nbridges=nbridges < max_threads

if ~keyword_set(delta) then begin
	if ~keyword_set(factor) then factor=4L
	if (stretchFlag) then delta=(xa[nx-1]-xa[0])/((nx-1.0)*float(factor)) else delta=1.0/factor
endif

if  keyword_set(twistFlag)  then twistFlag =1B else twistFlag =0B
if  keyword_set(RK4Flag)    then RK4Flag   =1B else RK4Flag   =0B
if  keyword_set(scottFlag)  then scottFlag =1B else scottFlag =0B
if  keyword_set(curlB_out)  then curlB_out =1B else curlB_out =0B
if  keyword_set(no_preview) then preview   =0B else preview   =1B
if  keyword_set(tmpB)       then tmpB      =1B else tmpB      =0B
if  keyword_set(RAMtmp)     then RAMtmp    =1B else RAMtmp    =0B
if ~keyword_set(tol)        then tol=10.0^(-4.)
if ~keyword_set(step)       then step=1.
if ~keyword_set(maxsteps)   then maxsteps=long(4*(nx+ny+nz)/step)
;----------------------------------------------------------------------------------------------
; the directory for output
cd, current = cdir
IF STRMID(cdir, STRLEN(cdir)-1) NE '/' THEN cdir=cdir+'/'

if  keyword_set(odir) then begin 
	preset_odir=1B
	IF STRMID(odir, STRLEN(odir)-1) NE '/' THEN odir=odir+'/'
endif else begin 
	preset_odir=0B
	odir= cdir+'qfactor/'
endelse
if ~file_test(odir) then file_mkdir, odir
;----------------------------------------------------------------------------------------------
; check the existence of curlB.sav
file_curlB=odir+'curlB.sav'
curlB_exist=file_test(file_curlB)
if curlB_out and (not curlB_exist) then curlB_out_int=1L else curlB_out_int=0L
;----------------------------------------------------------------------------------------------
; the temporary directory for the data transmission between Fortran and IDL
if RAMtmp then tmp_dir='/dev/shm/tmp/' else tmp_dir= odir+'tmp/'
if ~file_test(tmp_dir) then file_mkdir, tmp_dir
dummy=file_search(tmp_dir,'*.{txt,bin}',count=nf)
if nf ne 0 then file_delete, dummy
;----------------------------------------------------------------------------------------------
;  transmit data
get_lun,unit

openw,  unit, tmp_dir+'head.txt'
printf, unit, long(nx), long(ny), long(nz), long(nbridges), float(delta), long(maxsteps)
printf, unit, float(xreg), float(yreg), float(zreg), float(step), float(tol)
printf, unit, long(twistFlag), long(RK4flag), long(scottFlag), long(csflag), curlB_out_int
close,  unit

openw,  unit, tmp_dir+'b3d.bin'
if keyword_set(tmpB) then $
     writeu, unit, float(temporary(Bx)), float(temporary(By)), float(temporary(Bz)) $
else writeu, unit, float(Bx), float(By), float(Bz)
close,  unit

if (stretchFlag) then begin
;txt file could introduce errors to values while binary file doesn't
	openw,  unit, tmp_dir+'xa.bin'
	writeu, unit, float(xa)
	close,  unit

	openw,  unit, tmp_dir+'ya.bin'
	writeu, unit, float(ya)
	close,  unit

	openw,  unit, tmp_dir+'za.bin'
	writeu, unit, float(za)
	close,  unit
endif
;----------------------------------------------------------------------------------------------
; calculate in Fortran
cd, tmp_dir
tnow=systime(1)
spawn, 'qfactor.x' ; if not known by the system, specify the path
tend=systime(1)
cd, cdir
tcalc=tend-tnow

if (tcalc ge 3600.0) then begin
	time_elapsed=string(tcalc/3600.0,format='(f0.2)')+' hours'
endif else begin 
	if (tcalc ge 60.0) then time_elapsed=string(tcalc/60.0,format='(f0.2)')+' minutes' $
	                   else time_elapsed=string(tcalc,     format='(f0.2)')+' seconds'
endelse
print, time_elapsed+' elapsed for calculation' 

; ################################### retrieving results ###################################################### 
qx=0L & qy=0L & qz=0L & q1=0L & q2=0L
z0Flag=0 & vFlag=0 & cFlag=0
openr,  unit, tmp_dir+'tail.txt'
readf,  unit, qx, qy, qz, q1, q2
readf,  unit, z0Flag, cFlag, vFlag
close,  unit
;----------------------------------------------------------------------------------------------
; the name of .sav file
if keyword_set(fstr) then begin
	preset_fstr=1B
endif else begin 

	preset_fstr=0B
	
	if abs(delta - round(delta)) lt 0.0001 then begin 
		delta_str='delta'+string(delta,'(i0)')
	endif else begin
		if abs(delta*10 - round(delta*10)) lt 0.0001 then begin
			delta_str='delta'+string(delta,'(f0.1)') 
		endif else begin
			if abs(delta*100 - round(delta*100)) lt 0.0001 then $			
				delta_str='delta'+string(delta,'(f0.2)') $
			else	delta_str='delta'+string(delta,'(f0.3)')
		endelse
	endelse
	
	cut_str=''
	
	if vflag then begin
		head_str='q3d_' 		
	endif else  begin
		head_str='qcs_'
		 
		xFlag=qx eq 1 ; calculate Q in a cross section parallel to x=0  
		yFlag=qy eq 1 ; calculate Q in a cross section parallel to y=0
		zFlag=qz eq 1 ; calculate Q in a cross section parallel to the photosphere
		
		if xFlag then cut_str0='_x'
		if yFlag then cut_str0='_y'
		if zFlag then cut_str0='_z'
		
		if xFlag then cut_coordinate=xreg[0]
		if yFlag then cut_coordinate=yreg[0]
		if zFlag then cut_coordinate=zreg[0]

		if (xFlag or yFlag or zFlag) then begin
			if abs(cut_coordinate - round(cut_coordinate)) lt 0.0001 then begin 
				cut_str=cut_str0+string(cut_coordinate,'(i0)')
			endif else begin
				if abs(cut_coordinate*10 - round(cut_coordinate*10)) lt 0.0001 then $
					cut_str=cut_str0+string(cut_coordinate,'(f0.1)') $
				else 	cut_str=cut_str0+string(cut_coordinate,'(f0.2)')
			endelse
		endif
	endelse
	
	fstr = head_str + delta_str + cut_str
endelse

file_sav=odir+fstr+'.sav'
;----------------------------------------------------------------------------------------------
; save curlB
if curlB_out and curlB_exist then print, "'"+file_curlB+"' exist already"
if curlB_out_int then begin
	curlB=fltarr(3, nx, ny, nz)
	openr, unit, tmp_dir+'curlB.bin'
	readu, unit, curlB
	close, unit
		
	curlBx=reform(curlB[0,*,*,*])
	curlBy=reform(curlB[1,*,*,*])
	curlBz=reform(curlB[2,*,*,*])
	
	save, filename=file_curlB, curlBx, curlBy, curlBz
		
	; release the memory of curlB	
	dummy=(temporary(curlB))[0]+(temporary(curlBx))[0]+(temporary(curlBy))[0]+(temporary(curlBz))[0]
	
	print, "curlB is saved in '"+file_curlB+"'"
endif
;----------------------------------------------------------------------------------------------
; mark the area for calculation on the magnetogram
if preview then begin
	
	if (stretchFlag) then begin
		nx_mag=0L & ny_mag=0L & delta_mag=0.0
		openr, unit, tmp_dir+'mag_info.txt'
		readf, unit, nx_mag, ny_mag, delta_mag
		close, unit
		
		magnetogram=fltarr(nx_mag, ny_mag)
		openr, unit, tmp_dir+'magnetogram.bin'
		readu, unit, magnetogram
		close, unit
	endif else begin
		nx_mag=nx
		ny_mag=ny
	endelse

	cur_device=!D.name
	SET_PLOT, 'Z' 	
	DEVICE, SET_RESOLUTION=[nx_mag, ny_mag]

	scale_top= max(abs(magnetogram))/2.0 < 1000.0
	tv,bytscl(temporary(magnetogram),min=-scale_top,max=scale_top,/nan)
		
	if (stretchFlag) then begin
		if csflag then plots, (xreg-xa(0))/delta_mag, (yreg-ya(0))/delta_mag, /dev $
		          else plots, ([xreg[0],xreg[1],xreg[1],xreg[0],xreg[0]]-xa(0))/delta_mag, $
			              ([yreg[0],yreg[0],yreg[1],yreg[1],yreg[0]]-ya(0))/delta_mag, /dev
	endif else begin
		if csflag then plots, xreg, yreg, /dev $
		          else plots, [xreg[0],xreg[1],xreg[1],xreg[0],xreg[0]], [yreg[0],yreg[0],yreg[1],yreg[1],yreg[0]], /dev
	endelse
	
	im=TVRD(0,0,nx_mag,ny_mag)
	write_png, odir+fstr+'_bz.png', im
	set_plot, cur_device

	if (stretchFlag) then length_top=sqrt((xa(nx-1)-xa(0))^2.0+(ya(ny-1)-ya(0))^2.0+(za(nz-1)-za(0))^2.0) $
	                 else length_top=sqrt(nx^2.0+ny^2.0+nz^2.0)
		         
; load doppler color table
	r_doppler=[bindgen(127)*2B, REPLICATE(255B, 129)]
	b_doppler=reverse(r_doppler)
	g_doppler=[r_doppler[0:127], b_doppler[128:255]]

endif


; Q at the photosphere ----------------------------------------------------------------------------------------------
IF z0Flag THEN BEGIN
; read the output of qfactor.x
	q=fltarr(q1,q2)
	slogq=fltarr(q1,q2)
  	length=fltarr(q1,q2)
	Bnr=fltarr(q1,q2)
	rF=fltarr(3,q1,q2)
	rboundary=bytarr(q1,q2)
	
	openr, unit, tmp_dir+'qfactor0.bin'
	readu, unit, q, length, Bnr, rF, rboundary
	close, unit
	
	openr, unit, tmp_dir+'slogq.bin'
	readu, unit, slogq
	close, unit
  
  	ss_rb=where(rboundary ne 1)
	slogq_orig=slogq
	if(ss_rb[0] ne -1) then slogq_orig[ss_rb]=0.0

; preview images
  	if preview then begin  		
		slogq_tmp=slogq_orig
		ss=WHERE(FINITE(slogq_tmp, /NAN))
		if (ss[0] ne -1) then slogq_tmp[ss]=0.0  ; value of NaN should be colored with white here
		im=bytscl(slogq_tmp,min=-5,max=5,/nan)
		write_png, odir+fstr+'_slogq_orig.png', im, r_doppler, g_doppler, b_doppler
		
		slogq_tmp=slogq
		ss=WHERE(FINITE(slogq_tmp, /NAN))
		if (ss[0] ne -1) then slogq_tmp[ss]=0.0  ; value of NaN should be colored with white here
		im=bytscl(slogq_tmp, min=-5,max=5,/nan)
		write_png, odir+fstr+'_slogq.png', im, r_doppler, g_doppler, b_doppler	
			
		im=bytscl(length,min=0,max=length_top,/nan)  
		write_png, odir+fstr+'_length.png', im
		
		ss=where(Bnr eq 0.0)
		Bnr_tmp=Bnr
		if (ss[0] ne -1) then Bnr_tmp[ss]=1.   ; to avoid "Program caused arithmetic error: Floating divide by 0"
		im=bytscl(alog10(Bnr_tmp),min=-2,max=2,/nan) 
		write_png, odir+fstr+'_lg(Bnr).png', im, r_doppler, g_doppler, b_doppler
	endif
		
; q_perp map
	if scottFlag then begin
		q_perp=fltarr(q1,q2)
		slogq_perp=fltarr(q1,q2)
		openr, unit, tmp_dir+'q_perp.bin'
		readu, unit, q_perp
		close, unit
		
		openr, unit, tmp_dir+'slogq_perp.bin'
		readu, unit, slogq_perp
		close, unit
		
		slogq_perp_orig=slogq_perp
		if (ss_rb[0] ne -1) then slogq_perp_orig[ss_rb]=0.0
		
		if preview then begin
			slogq_tmp=slogq_perp_orig
			ss=WHERE(FINITE(slogq_tmp, /NAN))
			if(ss[0] ne -1) then slogq_tmp[ss]=0.0  ; value of NaN should be colored with white here
			im=bytscl(slogq_tmp,min=-5,max=5,/nan)
			write_png, odir+fstr+'_slogq_perp_orig.png', im, r_doppler, g_doppler, b_doppler	
			
			slogq_tmp=slogq_perp
			ss=WHERE(FINITE(slogq_tmp, /NAN))
			if(ss[0] ne -1) then slogq_tmp[ss]=0.0  ; value of NaN should be colored with white here
			im=bytscl(slogq_tmp,min=-5,max=5,/nan)
			write_png, odir+fstr+'_slogq_perp.png', im, r_doppler, g_doppler, b_doppler
		endif
	endif  
  
; twist map
	if twistFlag then begin
		twist=fltarr(q1,q2)
		openr, unit, tmp_dir+'twist.bin'
		readu, unit, twist
		close, unit
		if preview then begin
			im=bytscl(twist,min=-2,max=2,/nan)
			write_png, odir+fstr+'_twist.png', im, r_doppler, g_doppler, b_doppler
		endif
	endif

; save results 
	if scottFlag then begin
	 	if twistFlag then save, filename=file_sav, slogq, slogq_orig, q, length, Bnr, rboundary, xreg, yreg, zreg, delta, $
	 	                  rF, q_perp, slogq_perp, slogq_perp_orig, twist $ 
		             else save, filename=file_sav, slogq, slogq_orig, q, length, Bnr, rboundary, xreg, yreg, zreg, delta, $
		                  rF, q_perp, slogq_perp, slogq_perp_orig
	endif else begin
		if twistFlag then save, filename=file_sav, slogq, slogq_orig, q, length, Bnr, rboundary, xreg, yreg, zreg, delta, rF, twist $ 
		             else save, filename=file_sav, slogq, slogq_orig, q, length, Bnr, rboundary, xreg, yreg, zreg, delta, rF
	endelse
ENDIF 

; Q at the cross section ----------------------------------------------------------------------------------------------
IF cFlag THEN BEGIN
; read the output of qfactor.x
	qcs=fltarr(q1,q2)
	length=fltarr(q1,q2)
	rsF=fltarr(3,q1,q2)
	reF=fltarr(3,q1,q2)
	rsboundary=bytarr(q1,q2)		
	reboundary=bytarr(q1,q2)
	openr, unit, tmp_dir+'qcs.bin'
	readu, unit, qcs, length, rsF, reF, rsboundary, reboundary
	close, unit
	
	logq=alog10(qcs>1.)
	qcs_orig=qcs
	ss1=where(rsboundary ne 1)
	ss2=where(reboundary ne 1)
	if (ss1[0] ne -1) then qcs_orig[ss1]=!values.F_NAN
	if (ss2[0] ne -1) then qcs_orig[ss2]=!values.F_NAN
	
	if preview then begin
		im=bytscl(alog10(qcs_orig>1.),min=0,max=5,/nan)
		write_png, odir+fstr+'_logq_orig.png', im
		
		im=bytscl(logq,min=0,max=5,/nan)
		write_png, odir+fstr+'_logq.png', im
		im=bytscl(length,min=0,max=length_top,/nan)
		write_png, odir+fstr+'_length.png', im
	endif

	if scottFlag then begin
		q_perp=fltarr(q1,q2)
		openr, unit, tmp_dir+'q_perp.bin'
		readu, unit, q_perp
		close, unit
		
		q_perp_orig=q_perp
		if (ss1[0] ne -1) then q_perp_orig[ss1]=!values.F_NAN
		if (ss2[0] ne -1) then q_perp_orig[ss2]=!values.F_NAN
		logq_perp_orig=alog10(q_perp_orig>1)
		logq_perp=alog10(q_perp>1.)		
		
		if preview then begin
			im=bytscl(logq_perp_orig,min=0,max=5,/nan)
			write_png, odir+fstr+'_logq_perp_orig.png', im
			im=bytscl(logq_perp,min=0,max=5,/nan)
			write_png, odir+fstr+'_logq_perp.png', im
		endif
	endif

	if twistFlag then begin	
		twist=fltarr(q1,q2)
		openr, unit, tmp_dir+'twist.bin'
		readu, unit, twist
		close, unit
		if preview then begin
			im=bytscl(twist,min=-2,max=2,/nan)
			write_png, odir+fstr+'_twist.png', im, r_doppler, g_doppler, b_doppler
		endif
	endif	   

; save results 
	if scottFlag then begin
		if twistFlag then save, filename=file_sav, qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, delta, csFlag, rsF, reF, q_perp, q_perp_orig, twist $
		             else save, filename=file_sav, qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, delta, csFlag, rsF, reF, q_perp, q_perp_orig	
	endif else begin
		if twistFlag then save, filename=file_sav, qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, delta, csFlag, rsF, reF, twist $
		             else save, filename=file_sav, qcs, qcs_orig, length, rsboundary, reboundary, xreg, yreg, zreg, delta, csFlag, rsF, reF
	endelse
	
ENDIF

; Q in the 3D box volume ----------------------------------------------------------------------------------------------
IF vflag THEN BEGIN
; read the output of qfactor.x
	q3d=fltarr(qx,qy,qz)
	openr, unit, tmp_dir+'q3d.bin'
	readu, unit, q3d
	close, unit
	
	rboundary3d=bytarr(qx,qy,qz)	
	openr, unit, tmp_dir+'rboundary3d.bin'
	readu, unit, rboundary3d
	close, unit
	
	if scottFlag then begin
		q_perp3d=fltarr(qx,qy,qz)
		openr, unit, tmp_dir+'q_perp3d.bin'
		readu, unit, q_perp3d
		close, unit
	endif  

	if twistFlag then begin
		twist3d=fltarr(qx,qy,qz)
		openr, unit, tmp_dir+'twist3d.bin'
		readu, unit, twist3d
		close, unit
	endif
	
	if (zreg[0] eq 0.0) and preview then begin
		slogq=fltarr(q1,q2)	
		openr, unit, tmp_dir+'slogq.bin'
		readu, unit, slogq
		close, unit		
		
		slogq_tmp=slogq
		ss=WHERE(FINITE(slogq, /NAN))		
		if(ss[0] ne -1) then slogq_tmp[ss]=0.0		
		im=bytscl(slogq_tmp,min=-5,max=5,/nan)
		write_png, odir+fstr+'_z0_slogq.png', im, r_doppler, g_doppler, b_doppler
	endif 
		
; save results 
	if scottFlag then begin
		if twistFlag then save, filename=file_sav, q3d, rboundary3d, xreg, yreg, zreg, delta, q_perp3d, twist3d $
		             else save, filename=file_sav, q3d, rboundary3d, xreg, yreg, zreg, delta, q_perp3d
	endif else begin
		if twistFlag then save, filename=file_sav, q3d, rboundary3d, xreg, yreg, zreg, delta, twist3d $
		             else save, filename=file_sav, q3d, rboundary3d, xreg, yreg, zreg, delta	
	endelse
ENDIF 

;----------------------------------------------------------------------------------------------	
; house keeping
file_delete, tmp_dir, /recursive
free_lun, unit, /force
if ~preset_odir then dummy=temporary(odir)
if ~preset_fstr then dummy=temporary(fstr)

print, "Results are saved in '"+file_sav+"'"
END
