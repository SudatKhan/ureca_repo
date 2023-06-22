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
    real :: mp1 ! mass of planet 1
    real :: mp2 ! mass of planet 2
    real :: rp1 ! radius of planet 1
    real :: rp2 ! radius of planet 2
    real :: theta 
    real :: phi
    real :: x1 ! x-position of planet 1
    real :: y1 ! y-position of planet 1
    real :: v_x1 ! x-component of planet 1 velocity
    real :: v_y1  ! y-component of planet 1 velocity
    real :: x2 ! x-position of planet 2
    real :: y2 ! y-position of planet 2
    real :: v_x2 ! x-component of planet 2 velocity
    real :: v_y2  ! y-component of planet 2 velocity
    integer :: n_id !number of particles in the dump file
 
   contains
   
   subroutine get_centreofmass_part(xcom,vcom,npart_start,npart_end,xyzh,vxyzu)
    use io,       only:id,master
    use dim,      only:maxphase,maxp
    use part,     only:massoftype,iamtype,iphase,igas,isdead_or_accreted
    use mpiutils, only:reduceall_mpi
    real,         intent(out) :: xcom(3),vcom(3)
    integer,      intent(in)  :: npart_start, npart_end
    real,         intent(in)  :: xyzh(:,:),vxyzu(:,:)
 
    integer :: i,itype
    real :: xi,yi,zi,hi
    real(kind=8) :: xpos,ypos,zpos,vxpos,vypos,vzpos
    real(kind=8) :: dm,pmassi,totmass
   
    xpos  = 0.d0
    ypos  = 0.d0
    zpos  = 0.d0
    vxpos = 0.d0
    vypos = 0.d0
    vzpos = 0.d0
    totmass = 0.d0
    pmassi = massoftype(igas)
    do i=npart_start,npart_end
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       hi = xyzh(4,i)
       if (.not.isdead_or_accreted(hi)) then
          if (maxphase==maxp) then
             itype = iamtype(iphase(i))
             if (itype > 0) then ! avoid problems if called from ICs
                pmassi = massoftype(itype)
             else
                pmassi = massoftype(igas)
             endif
          endif
          totmass = totmass + pmassi
          xpos    = xpos  + pmassi*xi
          ypos    = ypos  + pmassi*yi
          zpos    = zpos  + pmassi*zi
          vxpos   = vxpos + pmassi*vxyzu(1,i)
          vypos   = vypos + pmassi*vxyzu(2,i)
          vzpos   = vzpos + pmassi*vxyzu(3,i)
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
    endif
    xcom = xcom*dm
    vcom = vcom*dm
   
    return
   end subroutine get_centreofmass_part
 
   subroutine reset_centreofmass_part(npart_start,npart_end,xyzh,vxyzu)
    use io, only:iprint
    integer, intent(in)    :: npart_start, npart_end
    real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 
    real :: xcom(3),vcom(3),xcomold(3),vcomold(3)
    integer :: i
   
       call get_centreofmass_part(xcom,vcom,npart_start, npart_end, xyzh,vxyzu)
 
    xcomold = xcom
    vcomold = vcom
    do i=npart_start,npart_end
       xyzh(1:3,i) = xyzh(1:3,i) - xcom(1:3)
       vxyzu(1:3,i) = vxyzu(1:3,i) - vcom(1:3)
    enddo
    
       call get_centreofmass_part(xcom,vcom,npart_start, npart_end, xyzh,vxyzu)
    
    write(iprint,"(' reset CofM: (',3(es9.2,1x),') -> (',3(es9.2,1x),')')") xcomold,xcom
   
    return
   end subroutine reset_centreofmass_part
 
   subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
    use centreofmass
    use part,           only:xyzmh_ptmass,vxyz_ptmass,nptmass, ihacc, ihsoft, imacc, ispinx, ispiny, ispinz
    use units,     only:umass,udist
    use prompting,      only:prompt
    use physcon,        only:pi
    integer,  intent(inout) :: npart
    integer,  intent(inout) :: npartoftype(:)
    real,     intent(inout) :: massoftype(:)
    real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
    character(len=120)      :: filename
    integer                 :: i,ierr
    logical                 :: iexist
    real                    :: Lx,Ly,Lz,L,Lp,Ltot(3)
    real                    :: rp_p1,rp_p2
    real                    :: x,y,z,vx,vy,vz
    real                    :: x01,y01,vx01,vy01,r01,U1
    real                    :: x02,y02,vx02,vy02,r02, U2
   
   ! Default parameters
    Mh    = 1.e1   ! star mass
    rs    = 1
    theta = 0.     ! stellar tilting along x
    phi   = 0.     ! stellar tilting along y
 
    Mp1    = 1.     ! 
    rp1    = 1.     ! 
    ecc1   = 0.                    
    a1    = 10 
    x1    = 0
    y1 = 0
    v_x1 = 0
    v_y1 = 0
    
    Mp2    = 1.     ! 
    rp2    = 1.     ! 
    ecc2   = 0.                    
    a2    = 10 
    x2 = 0
    y2 = 0
    v_x2 = 0
    v_y2 = 0
 
    n_id = 0
    
 
    U1 =  Mh+Mp1
    rp_p1 = a1*(1-ecc1)                       ! pericenter distance
 
    U2 =  Mh+Mp2
    rp_p2 = a2*(1-ecc2)                       ! pericenter distance
    
    !-- Read runtime parameters from binaryparams file
    filename = 'binary'//'.binaryparams'                                ! moddump should really know about the output file prefix...
    inquire(file=filename,exist=iexist)
    if (iexist) call read_setupfile(filename,ierr)
    if (.not. iexist .or. ierr /= 0) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantommoddump'
       stop
    endif
    call reset_centreofmass_part(1,n_id,xyzh,vxyzu)
 
    call reset_centreofmass_part(n_id+1,npart,xyzh,vxyzu)
   
    phi   = 0.
    theta = 0.
   
    call get_angmom(ltot,npart,xyzh,vxyzu)
    Lx = ltot(1)
    Ly = ltot(2)
    Lz = ltot(3)
    Lp = sqrt(Lx**2.0+Lz**2.0)
    if (Lx > 0.) then
       phi=acos(Lz/Lp)
    elseif (Lx < 0.) then
       phi=-acos(Lz/Lp)
    endif

   !--Rotate the star so the momentum lies in the yz plan
    print*,'tilting along y axis: ',(phi*180/pi),'degrees'
    do i=1,npart
       x=xyzh(1,i)
       z=xyzh(3,i)
       xyzh(1,i)=x*cos(-phi)+z*sin(-phi)
       xyzh(3,i)=-x*sin(-phi)+z*cos(-phi)
       vx=vxyzu(1,i)
       vz=vxyzu(3,i)
       vxyzu(1,i)=vx*cos(-phi)+vz*sin(-phi)
       vxyzu(3,i)=-vx*sin(-phi)+vz*cos(-phi)
    enddo
   
   !--Recheck the stellar angular momentum
    call get_angmom(ltot,npart,xyzh,vxyzu)
    lx = ltot(1)
    ly = ltot(2)
    lz = ltot(3)
    L  = sqrt(Lx**2.0+Ly**2.0+Lz**2.0)
    if (Ly < 0.) then
       theta=acos(Lz/L)
    elseif (Ly > 0.) then
       theta=-acos(Lz/L)
    endif

   !--Rotate the star so the momentum lies along the z axis
    print*, 'tilting along x axis: ',(theta*180/pi),'degrees'
    do i=1,npart
       y=xyzh(2,i)
       z=xyzh(3,i)
       xyzh(2,i)=y*cos(-theta)-z*sin(-theta)
       xyzh(3,i)=y*sin(-theta)+z*cos(-theta)
       vy=vxyzu(2,i)
       vz=vxyzu(3,i)
       vxyzu(2,i)=vy*cos(-theta)-vz*sin(-theta)
       vxyzu(3,i)=vy*sin(-theta)+vz*cos(-theta)
    enddo

   !--Recheck the stellar angular momentum
    call get_angmom(ltot,npart,xyzh,vxyzu)
    print*,'Stellar spin should now be along the z axis.'
   
    x01    = x1
    y01    = y1
    vx01   = v_x1
    vy01   = v_y1
   
   
    x02    = x2
    y02    = y2
    vx02   = v_x2
    vy02   = v_y2
   
    !--Tilting the star
    theta=theta*pi/180.0
    phi=phi*pi/180.0
   
    if (theta  /=  0.) then
       do i=1,npart
          y=xyzh(2,i)
          z=xyzh(3,i)
          xyzh(2,i)=y*cos(theta)-z*sin(theta)
          xyzh(3,i)=y*sin(theta)+z*cos(theta)
          vy=vxyzu(2,i)
          vz=vxyzu(3,i)
          vxyzu(2,i)=vy*cos(theta)-vz*sin(theta)
          vxyzu(3,i)=vy*sin(theta)+vz*cos(theta)
       enddo
    endif
    if (phi  /=  0.) then
       do i=1,npart
          x=xyzh(1,i)
          z=xyzh(3,i)
          xyzh(1,i)=x*cos(phi)+z*sin(phi)
          xyzh(3,i)=-x*sin(phi)+z*cos(phi)
          vx=vxyzu(1,i)
          vz=vxyzu(3,i)
          vxyzu(1,i)=vx*cos(phi)+vz*sin(phi)
          vxyzu(3,i)=-vx*sin(phi)+vz*cos(phi)
       enddo
    endif
   
    !--Putting planet into orbit
    do i = 1, n_id
       xyzh(1,i)  = xyzh(1,i)  + x01
       xyzh(2,i)  = xyzh(2,i)  + y01
       vxyzu(1,i) = vxyzu(1,i) + vx01
       vxyzu(2,i) = vxyzu(2,i) + vy01
    enddo
 
    do i = n_id+1, npart
       xyzh(1,i)  = xyzh(1,i)  + x02
       xyzh(2,i)  = xyzh(2,i)  + y02
       vxyzu(1,i) = vxyzu(1,i) + vx02
       vxyzu(2,i) = vxyzu(2,i) + vy02
    enddo
   
    theta = theta*pi/180.
    phi   = phi*pi/180.
   
   ! -- create central star as sink particle
    nptmass = 1
   
    xyzmh_ptmass(1:3,1) = 0
    vxyz_ptmass(1:3,1) = 0
    xyzmh_ptmass(4,1) = Mh
    xyzmh_ptmass(ihacc,1) = rs
    xyzmh_ptmass(ihsoft,1) = 0.0
    xyzmh_ptmass(imacc,1) = 0.0
    xyzmh_ptmass(ispinx,1) = 0.0
    xyzmh_ptmass(ispiny,1) = 0.0
    xyzmh_ptmass(ispinz,1) = 0.0
 
   
    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
   
    write(*,'(a)') "======================================================================"
    write(*,'(a,Es12.5,a)') ' Tilting along x     = ',theta,' degrees'
    write(*,'(a,Es12.5,a)') ' Tilting along y     = ',phi,' degrees'
    write(*,'(a,Es12.5,a)') ' Pericenter distance = ',rp_p1,' R_sun'
    write(*,'(a,Es12.5,a)') ' Pericenter distance = ',rp_p2,' R_sun'
 
    write(*,'(a,Es12.5,a)') ' Radius of star      = ',rp1,' R_sun'
    write(*,'(a,Es12.5,a)') ' Initial SMA         = ',a1,' R_sun'
    write(*,'(a,Es12.5,a)') ' Stellar mass        = ',Mp1,' M_sun'
 
    write(*,'(a,Es12.5,a)') ' Radius of star      = ',rp2,' R_sun'
    write(*,'(a,Es12.5,a)') ' Initial SMA         = ',a2,' R_sun'
    write(*,'(a,Es12.5,a)') ' Stellar mass        = ',Mp2,' M_sun'
   
    write(*,'(a)') "======================================================================"
   
    return
   end subroutine modify_dump
   
   !---Read/write setup file
   subroutine write_setupfile(filename)
    use infile_utils, only:write_inopt
    character(len=*), intent(in) :: filename
    integer, parameter :: iunit = 20
   
    print "(a)",' writing moddump params file '//trim(filename)
    open(unit=iunit,file=filename,status='replace',form='formatted')
    write(iunit,"(a)") '# parameters file for a TDE phantommodump'
    call write_inopt(n_id,    'n_id',    'last id of sph particle in planet 1',                     iunit)
    call write_inopt(mh,    'mh',    'mass of star (code units)',     iunit)
    call write_inopt(rs,    'rs',    'accretion radius (code units)', iunit)
 
    call write_inopt(ecc1,   'ecc1',   'eccentricity of planet 1',                                        iunit)
    call write_inopt(mp1,    'mp1',    'mass of planet 1      (code units)',                     iunit)
    call write_inopt(rp1,    'rp1',    'radius of planet 1    (code units)',                     iunit)
    call write_inopt(a1,      'a1',    'initial SMA of planet 1',                                         iunit)
    call write_inopt(x1,      'x1',    'x-position of planet 1',                                 iunit)
    call write_inopt(y1,      'y1',    'y-position of planet 1',                                 iunit)
    call write_inopt(v_x1,    'v_x1',  'x-component of planet 1 velocity', iunit)
    call write_inopt(v_y1, 'v_y1', 'y-component of planet 1 velocity', iunit)
 
    call write_inopt(ecc2,   'ecc2',   'eccentricity of planet 2',                                        iunit)
    call write_inopt(nu2,    'nu2',    'true anomaly of planet 2',                                        iunit)
    call write_inopt(mp2,    'mp2',    'mass of planet 2     (code units)',                     iunit)
    call write_inopt(rp2,    'rp2',    'radius of planet 2    (code units)',                     iunit)
    call write_inopt(a2,    'a2',    'initial SMA of planet 2',                                         iunit)
    call write_inopt(x2, 'x2', 'x-position of planet 2', iunit)
    call write_inopt(y2, 'y2', 'y-position of planet 2', iunit)
    call write_inopt(v_x2, 'v_x2', 'x-component of planet 2 velocity', iunit)
    call write_inopt(v_y2, 'v_y2', 'y-component of planet 2 velocity', iunit)
 
    close(iunit)
   
   end subroutine write_setupfile
   
   subroutine read_setupfile(filename,ierr)
    use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
    use io,           only:error
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: ierr
    integer, parameter :: iunit = 21
    integer :: nerr
    type(inopts), allocatable :: db(:)
   
    print "(a)",'reading setup options from '//trim(filename)
    nerr = 0
    ierr = 0
    call open_db_from_file(db,filename,iunit,ierr)
    call read_inopt(n_id,    'n_id',    db,min=0,errcount=nerr)
    call read_inopt(mh,    'mh',    db,min=0.,errcount=nerr)
    call read_inopt(rs,    'rs',    db,min=0.,errcount=nerr)
 
    call read_inopt(ecc1,   'ecc1',   db,min=0.,errcount=nerr)
    call read_inopt(mp1,    'mp1',    db,min=0.,errcount=nerr)
    call read_inopt(rp1,    'rp1',    db,min=0.,errcount=nerr)
    call read_inopt(a1,    'a1',    db,min=0.,errcount=nerr)
    call read_inopt(x1, 'x1', db,min=0.,errcount=nerr)
    call read_inopt(y1, 'y1', db,min=0.,errcount=nerr)
    call read_inopt(v_x1, 'v_x1', db,min=0.,errcount=nerr)
    call read_inopt(v_y1, 'v_y1', db,min=0.,errcount=nerr)

 
    call read_inopt(ecc2,   'ecc2',   db,min=0.,errcount=nerr)
    call read_inopt(nu2,    'nu2',    db,min=0.,errcount=nerr)
    call read_inopt(mp2,    'mp2',    db,min=0.,errcount=nerr)
    call read_inopt(rp2,    'rp2',    db,min=0.,errcount=nerr)
    call read_inopt(a2,    'a2',    db,min=0.,errcount=nerr)
    call read_inopt(x2, 'x2', db,min=0.,errcount=nerr)
    call read_inopt(y2, 'y2', db,min=0.,errcount=nerr)
    call read_inopt(v_x2, 'v_x2', db,min=0.,errcount=nerr)
    call read_inopt(v_y2, 'v_y2', db,min=0.,errcount=nerr)  
   
    call close_db(db)
    if (nerr > 0) then
       print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
       ierr = nerr
    endif
   
   end subroutine read_setupfile
   
   subroutine get_angmom(ltot,npart,xyzh,vxyzu)
    real, intent(out)   :: ltot(3)
    integer, intent(in) :: npart
    real, intent(in)    :: xyzh(:,:), vxyzu(:,:)
    integer :: i
    real    :: L
   
    ltot = 0.
    do i=1,npart
       ltot(1) = ltot(1)+xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i)
       ltot(2) = ltot(2)+xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i)
       ltot(3) = ltot(3)+xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i)
    enddo
   
    L = sqrt(dot_product(ltot,ltot))
   
    print*,''
    print*,'Checking angular momentum orientation and magnitude...'
    print*,'Angular momentum is L = (',ltot(1),ltot(2),ltot(3),')'
    print*,'Angular momentum modulus is |L| = ',L
    print*,''
   
   end subroutine get_angmom
   
   end module moddump