      program probe_template
c***********************************************************************
c     Hand-editable IR-structure probe (see probe-me-ir-structure skill).
c     Structure verified in-session on a 5-parton epem RR channel; all
c     UPPERCASE placeholders in the MARKED blocks must be replaced for
c     your ME (derive names/orders - me-naming-convention - or fit both
c     candidate orders: the wrong one is not rational). Solver refuses
c     degenerate systems and reports the residual.
c
c     Soft mode shown. For COLLINEAR limits: basis element must be
c     antenna*reducedME un-divided, each with its own set_map; use
c     rotp7 azimuthal averaging for gluon parents; scan x and judge by
c     coefficient stability (see SKILL.md).
c***********************************************************************
      use Process_mod
      use KinData_mod
      use Scales_mod
      use Mapping_mod
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
c --- EDIT: process commons (copy from the process's check program) ----
      common /koppel/ nf,nc,b0
      common /jetdef/ etminj,etmaxj,delrjj,rapmaxj,rapminj,jetalg
      common /inppar/ njets,nloop,norder
      common /nscale/ nscale
      common /CZFlav/ nf1, nf2
      common /BZFlav/ nfB1
      common /Zmass/ emz,ezwidth
      common/method/imethod
c ----------------------------------------------------------------------
      common/plotmode/iplot
      logical degen
      parameter (nb=6, nxs=3)
      dimension b(nb), am(nb,nb), rhs(nb), c(nb)
      dimension ipr(nb,2), xslist(nxs)

c --- EDIT: process init ----------------------------------------------
      call init_proc("MYPROC")
      call init_map()
      call setSqrts_proc(1000d0)
      iplot=2              ! nonzero or null.f ecuts stops the program
      imethod=2
      nscale=1
      emz=91.1876d0
      nf=5
      call setScales()
      call init_kin(5,10)
      njets=3
      nfB1=1
      nf1=1
      nf2=2
      etminj=0.05d0*sqrts_proc
      etmaxj=sqrts_proc/2d0
      rapmaxj=2d0
      rapminj=0d0
      delrjj=0.6d0
      jetalg=1
c --- EDIT: dipole candidates -----------------------------------------
c     one row per candidate (radiator, radiator) pair around the soft
c     parton ISOFT; list ALL candidates, the fit zeroes the spurious
      ipr(1,1)=I1
      ipr(1,2)=I2
      ipr(2,1)=I3
      ipr(2,2)=I4
      ipr(3,1)=I1
      ipr(3,2)=I3
      ipr(4,1)=I2
      ipr(4,2)=I4
      ipr(5,1)=I1
      ipr(5,2)=I4
      ipr(6,1)=I2
      ipr(6,2)=I3
      data xslist /1d-8, 1d-9, 1d-10/
      npt = 400
c ----------------------------------------------------------------------

      do ix=1,nxs
        xs = xslist(ix)
        am  = 0d0
        rhs = 0d0
        nacc = 0
        do ii=1,npt
c --- EDIT: limit generator + cuts + full ME + shared soft map --------
c     ISOFT = the soft parton; the set_map cluster surrounds it; the
c     reduced ME uses the mapped indices j1..jN from Mapping_mod
          em1=sqrts_proc*dsqrt(1d0-xs)
          call get_ss7(sqrts_proc,ISOFT,em1,I1,I2,I3,I4)
          ipass=1
          call ecuts_epem(1,7,ipass)
          if (ipass.ne.1) cycle
          ame = FULLME(3,4,5,6,7,2,1)
          call unset_map()
          ipass=1
          call set_map(7,6, (/IA,ISOFT,IB/), (/2,1,IREST/), ipass)
          if (ipass.ne.1) cycle
          red = REDME(j3,j4,j5,j6,j1,j2)
c ----------------------------------------------------------------------
          if (abs(red).lt.1d-30) cycle
          y = ame/red
          do kk=1,nb
            b(kk) = FullA30FF(ipr(kk,1),ISOFT,ipr(kk,2),7)
          end do
          scal = 0d0
          do kk=1,nb
            scal = scal + b(kk)**2
          end do
          if (scal.lt.1d-300) cycle
          scal = 1d0/scal
          do kk=1,nb
            do ll=1,nb
              am(kk,ll) = am(kk,ll) + scal*b(kk)*b(ll)
            end do
            rhs(kk) = rhs(kk) + scal*b(kk)*y
          end do
          nacc = nacc + 1
        end do

        write(*,*)
        write(*,*) 'x = ', xs, '  accepted points: ', nacc
        call gsolve(am,rhs,c,nb,degen,resid)
        if (degen) then
          write(*,*) 'DEGENERATE SYSTEM (rank < nb): refusing'
          write(*,*) 'coefficients. See SKILL.md (collinear basis).'
        else
          write(*,*) 'coefficients (residual = ', resid, '):'
          do kk=1,nb
            write(*,'(a,i2,a,i2,a,f16.8)')
     .        '   dipole (',ipr(kk,1),',',ipr(kk,2),') : ',c(kk)
          end do
        end if
      end do

      call destroy_kin()
      stop
      end program probe_template

************************************************************************
*     Gaussian elimination with degeneracy REFUSAL (no pivot clamping)
************************************************************************
      subroutine gsolve(a,rhs,x,n,degen,resid)
      implicit real*8(a-h,o-z)
      logical degen
      dimension a(n,n), rhs(n), x(n)
      dimension w(50,51), a0(50,50), r0(50)
      degen=.false.
      amax=0d0
      do i=1,n
        do j=1,n
          w(i,j)=a(i,j)
          a0(i,j)=a(i,j)
          if(abs(a(i,j)).gt.amax) amax=abs(a(i,j))
        end do
        w(i,n+1)=rhs(i)
        r0(i)=rhs(i)
      end do
      if(amax.lt.1d-300)then
        degen=.true.
        return
      endif
      tol = 1d-10*amax
      do k=1,n
        p=0d0
        ip=k
        do i=k,n
          if(abs(w(i,k)).gt.p)then
            p=abs(w(i,k))
            ip=i
          endif
        end do
        do j=k,n+1
          t=w(k,j)
          w(k,j)=w(ip,j)
          w(ip,j)=t
        end do
        if(abs(w(k,k)).lt.tol) then
          degen=.true.
          return
        endif
        do i=k+1,n
          f=w(i,k)/w(k,k)
          do j=k,n+1
            w(i,j)=w(i,j)-f*w(k,j)
          end do
        end do
      end do
      do i=n,1,-1
        s=w(i,n+1)
        do j=i+1,n
          s=s-w(i,j)*x(j)
        end do
        x(i)=s/w(i,i)
      end do
      resid=0d0
      do i=1,n
        s=-r0(i)
        do j=1,n
          s=s+a0(i,j)*x(j)
        end do
        resid=resid+s*s
      end do
      resid=dsqrt(resid)
      return
      end
