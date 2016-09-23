module fvlm
        implicit none
        real :: eps=1e-8
        real :: pi=3.141592653589793d0
        contains



subroutine cross(pv,v0,v1)
implicit none
real,dimension(3),intent(in) :: v0,v1
real,dimension(3),intent(out) :: pv
pv(1) = v0(2)*v1(3)-v0(3)*v1(2)
pv(2) = v0(3)*v1(1)-v0(1)*v1(3)
pv(3) = v0(1)*v1(2)-v0(2)*v1(1)
end subroutine cross

subroutine dot(ps,v0,v1)
        implicit none
        real,dimension(3),intent(in) :: v0,v1
        real,intent(out) :: ps
        ps=v0(1)*v1(1)+v0(2)*v1(2)+v0(3)*v1(3)
end subroutine  dot

subroutine vortxl(v,x,y,z,a,b,vort)
implicit none
real,dimension(3),intent(out) :: v
real,dimension(3),intent(in) :: a,b
real,intent(in) :: x,y,z,vort
real,dimension(3) :: r0,r1,r2,pv,pos
real :: ps1,ps2,normPV2,nrm1,nrm2,K
integer :: i

pos(1)=x
pos(2)=y
pos(3)=z

r0 = b-a
r1 = pos-a
r2 = pos-b
!r0(2) = yb-ya
!r0(3) = zb-za
!r1(1) = x-xa
!r1(2) = y-ya
!r1(3) = z-za
!r2(1) = x-xb
!r2(2) = y-yb
!r2(3) = z-zb



call cross(pv,r0,r1)
call dot(ps1,r0,r1)
call dot(ps2,r0,r2)

normPV2 = sum(pv(:)*pv(:))
nrm1 = sqrt(sum(r1(:)*r1(:)))
nrm2 = sqrt(sum(r2(:)*r2(:)))


if ((nrm1.lt.eps) .or. (nrm2.lt.eps) .or. (normPV2.lt.eps)) then
   K=0.
else
   K = (vort/(4.*pi*normPV2))*((ps1/nrm1)-(ps2/nrm2))
end if
v(1) = K*pv(1)
v(2) = K*pv(2)
v(3) = K*pv(3)

end subroutine vortxl


subroutine voring(v,x,y,z,a,b,c,d,gam)
implicit none
real,dimension(3),intent(out) :: v
real,dimension(3),intent(in)::a,b,c,d
real,intent(in) :: x,y,z
real,intent(in) :: gam
real,dimension(3) :: v1,v2,v3,v4

call vortxl(v1,x,y,z,a,b,gam)
call vortxl(v2,x,y,z,b,c,gam)
call vortxl(v3,x,y,z,c,d,gam)
call vortxl(v4,x,y,z,d,a,gam)

v=v1+v2+v3+v4

end subroutine voring


end module fvlm

