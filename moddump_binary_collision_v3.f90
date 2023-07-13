!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
! Editor: Sudat Khan

module moddump
    implicit none

    real :: mh ! central star mass
    real :: rs ! star accretion radius
    real :: mp1, mp2 ! mass of planet 1,2
    real :: rp1, rp2 ! radius of planet 1,2
    real :: theta
    real :: phi
    real :: a1, a2 ! semi-major axis of planet 1,2
    real :: ecc1, ecc2 ! eccentricity of planet 1,2 
    real :: x1, x2 ! x-position of planet 1,2
    real :: y1, y2 ! y-position of planet 1,2
    real :: v_x1, v_x2 ! x-component of planet 1,2 velocity
    real :: v_y1, v_y2  ! y-component of planet 1,2 velocity
    integer :: n_id !number of particles in the dump file

    contains

    subroutine get_centreofmass_part(xcom, vcom, npart_start, npart_end, xyzh, vxyzu)
        use io, only: id, master
        use dim, only: maxphase, maxp
        use part, only: massoftype, iamtype, iphase, igas, isdead_or_accreted
        use mpiutils, only: reduceall_mpi
        real, intent(out) :: xcom(3), vcom(3)
        integer, intent(in) :: npart_start, npart_end
        real, intent(in) :: xyzh(:,:), vxyzu(:,:)

        integer :: i, itype
        real :: xi, yi, zi, hi
        real(kind = 8) :: xpos, ypos, zpos, vxpos, vypos, vzpos
        real(kind = 8) :: dm, pmassi, totmass

        xpos = 0.d0
        ypos = 0.d0
        zpos = 0.d0
        vxpos = 0.d0
        vypos = 0.d0
        vzpos = 0.d0
        totmass = 0.d0
        pmassi = massoftype(igas)
        
        do i = npart_start, npart_end
            xi = xyzh(1,i)
            yi = xyzh(2,i)
            zi = xyzh(3,i)
            hi = xyzh(4,i)
            if (.not.isdead_or_accreted(hi)) then
                if (maxphase == maxp) then 
                    itype =  iamtype(iphase(i))
                    if (itype > 0) then
                        pmassi = massoftype(itype)
                    else
                        pmassi = massoftype(igas)
                    endif
                endif
                totmass = totmass + pmassi
                xpos = xpos + pmassi*xi
                ypos = ypos + pmassi*yi
                zpos = zpos + pmassi*zi
                vxpos = vxpos + pmassi*vxyzu(1,i)
                vypos = vypos + pmassi*vxyzu(2,i)
                vzpos = vzpos + pmassi*vxyzu(3,i)
            endif
        enddo

        xcom = (/xpos,ypos,zpos/)
        vcom = (/vxpos,vypos,vzpos/)
        xcom = reduceall_mpi('+',xcom)
        vcom = reduceall_mpi('+',vcom)
        totmass = reduceall_mpi('+',totmass)
        
        if (totmass > tiny(totmass)) then
            dm = 1.d0/totmass
        else
            dm = 0.d0
        enddo
        xcom = xcom*dm
        vcom = vcom*dm
        
        return
    end subroutine get_centreofmass_part

    subroutine reset_centreofmass_part(npart_start, npart_end, xyzh, vxyzu)
        use io, only:iprint
        integer, intent(in) :: npart_start, npart_end
        real, intent(inout) :: xyzh(:,:), vxyzu(:,:)
        real :: xcom(3), vcom(3), xcomold(3), vcomold(3)
        integer :: i

            call get_centreofmass_part(xcom,vcom,npart_start, npart_end, xyzh,vxyzu)