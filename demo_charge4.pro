pro bfield_charge4, bfile, scale, top

spawn, 'ls -d '+bfile, out, error_out
if (out ne '') then return

h_scale=scale/2.
length=4.5
h_length=length/2.
B_delta=length/scale

Bx=fltarr(scale+1,scale+1,top+1)
By=fltarr(scale+1,scale+1,top+1)
Bz=fltarr(scale+1,scale+1,top+1)

x1=[-1.5,0.0,-0.5]
x2=[-0.5,0.0,-0.5]
x3=[ 0.5,0.0,-0.5]
x4=[ 1.5,0.0,-0.5]

xi=fltarr(3,4)

xi[*,0]=x1
xi[*,1]=x4
xi[*,2]=x3
xi[*,3]=x2

epsilon=[1.,-1.,1.,-1.]


for i=-h_scale,0 do begin
for j=-h_scale,0 do begin
for k=0, top do begin
	x=[i,j,k]*B_delta
	Bp=fltarr(3)
	for s=0,3 do begin
		dx=x-xi[*,s]
		Bp=Bp+epsilon[s]*dx/(total(dx^2.))^1.5
	endfor

	Bx[i+h_scale,j+h_scale,k]=Bp[0]
	By[i+h_scale,j+h_scale,k]=Bp[1]
	Bz[i+h_scale,j+h_scale,k]=Bp[2]
endfor
endfor
endfor

Bx[h_scale+1:scale,0:h_scale,*]= reverse(Bx[0:h_scale-1,0:h_scale,*],1)
By[h_scale+1:scale,0:h_scale,*]=-reverse(By[0:h_scale-1,0:h_scale,*],1)
Bz[h_scale+1:scale,0:h_scale,*]=-reverse(Bz[0:h_scale-1,0:h_scale,*],1)

Bx[*,h_scale+1:scale,*]= reverse(Bx[*,0:h_scale-1,*],2)
By[*,h_scale+1:scale,*]=-reverse(By[*,0:h_scale-1,*],2)
Bz[*,h_scale+1:scale,*]= reverse(Bz[*,0:h_scale-1,*],2)

save, filename=bfile, Bx, By, Bz
print, bfile+' is saved'
end


;create a quadrupole field file
bfile='charge4.sav'
scale=90
top=60
bfield_charge4, bfile, scale, top
restore, bfile

;images of Figure 4 in Zhang, P., Chen, J.*, Liu, R. and Wang, C., 2022, ApJ, 937, 26
qfactor,Bx,By,Bz, fstr='method1_z0'
qfactor,Bx,By,Bz, fstr='method2_z0',/scott
qfactor,Bx,By,Bz, xreg=[0,scale],yreg=[scale/2,scale/2],zreg=[0,top/2],fstr='method1_y0'
qfactor,Bx,By,Bz, xreg=[0,scale],yreg=[scale/2,scale/2],zreg=[0,top/2],fstr='method2_y0',/scott

;An example of calculating in a cross section be tilted to x-axis and y-axis, 
;and with streched (actually uniformed) grids
length=4.5
B_delta=length/scale

xa=findgen(scale+1)*B_delta-length/2.0
ya=findgen(scale+1)*B_delta-length/2.0
za=findgen(top+1)*B_delta

qfactor, Bx, By, Bz, xa=xa, ya=ya, za=za, delta=B_delta/3, $
xreg=[-2,0], yreg=[1,0], zreg=[0,2], /csflag, $
fstr='tilted_cs', /rk4, step=2.0, odir= 'qfactor/', /twist, nbridges=4

;An example of calculating in a box volume, and exporting curlB arrays
qfactor, Bx, By, Bz, xreg=[scale/4,scale/2], yreg=[scale/9,scale/3], zreg=[top/4,top/2],$
delta=0.8, tol=1.0e-3, odir= 'qfactor',/curlB_out

end

