!----------------------------------------------------------------

      subroutine fnProposeParamUpdates(nparam,curval,
     $     valmin,valmax,step,logflag,
     $     parupdt,iseed,nopt,iopt)

      implicit none
      integer nparam, nopt, iopt(nopt)
      real*8 curval(nparam),valmin(nparam),valmax(nparam)
      real*8 parupdt(nparam),step(nparam)
      real*8 x, rv, rtn
      real*8 ran1
      real*8 SR_to_unit, SR_from_unit
      integer iseed
      integer i,j
      
      logical logflag(nparam)
      external SR_to_unit, SR_from_unit, ran1

       do j = 1, nopt
          i = iopt(j)
          parupdt(i)=curval(i)

          rv = ran1(iseed)

          rv = (rv - 0.50d0)*step(i)

! convert to a zero - one scale

          x = SR_to_unit(curval(i),valmin(i),valmax(i),
     $         logflag(i))

          x = x + rv

      if (x .lt. 0.0d0) x = 1.0d0 + x
      if (x .gt. 1.0d0) x = x - 1.0d0

! Do not use period boundary conditions here but rather re-smaple
c$$$      if (x .le. 0.0d0 .or. x .ge. 1.0d0) go to 101


! bring value back to original scale
      
         rtn = SR_from_unit(x,valmin(i),valmax(i),
     $        logflag(i))

         parupdt(i) = rtn

      enddo

      return
      end subroutine fnProposeParamUpdates


c----------------------------------------------------------------

      function SR_to_unit(y,ymin,ymax,logflag)

      implicit none
      real*8 y, ymin,ymax,rtn,SR_to_Unit
      logical logflag

      if (logflag) Then
         rtn = (log10(y) - log10(ymin)) /
     $        (log10(ymax)-log10(ymin))
      else
         rtn = (y - ymin)/(ymax - ymin)
      endif

      SR_to_unit = rtn
 
      return
      end function SR_to_unit

c----------------------------------------------------------------

      function SR_from_unit(x,ymin,ymax,logflag)
      
      implicit none
   
      real*8 x,ymin,ymax,rtn,SR_from_unit
      logical logflag


      if (logflag) Then
         rtn = ymin * 
     $        10.0**(x*(log10(ymax)-log10(ymin)))

      else
         rtn = ymin + (ymax-ymin)*x
      endif

      SR_from_unit = rtn

      return
      end function SR_from_unit


!--------------------------------------------------------------------------------
           
        subroutine weekly1D(ndata,ndays,nstep,imid,dsdt,pC,e_nonflu,x)
        
        implicit none
        integer ndata,nstep, ndays
        real*8 dsdt((ndays)*nstep)
        real*8 x(ndata), pC, e_nonflu
        integer i
        integer istart, iend
        integer imid((ndata+1))

        do i=1,ndata
c$$$           x(i) = dsdt(i*nstep*iday_per_week) - 
c$$$     $     dsdt(1+(i-1)*nstep*iday_per_week)

           iend   = imid(i+1)
           istart = imid(i)
           x(i) = dsdt(iend) - dsdt(istart)

        enddo

        x= x* pC +e_nonflu

        return
        end subroutine weekly1D

C--------------------------------------------------------------------------------
        function calcFit1D(y,gamay,x,wght,ndata)

        implicit none

        integer ndata, i
        real*8 y(ndata),x(ndata), gamay(ndata)
        real*8 wght(ndata)
        real*8 xi, yi,sum,val,CalcFit1D


c x is the simulated data
c y is the base profile

C calculate the P(yi,xi)


        sum = 0.0
        do i=1,ndata
           
           yi = y(i)
           xi = x(i)
           
           val = yi * log(xi) - xi - gamay(i)
           sum = sum + val  * wght(i) 
           
        enddo
        sum = -sum   
        
        CalcFit1D = sum 
           
        return
        end function CalcFit1D


!
! -------------------------------------------------------------------
!
      subroutine MH1D(fnewLLK,curLLK,curMin,nparam,ndata,
     $              iseed,iaccept,iadapt,curpars,savepar,parBest,
     $              rtn,rtnBest)


      implicit none

      integer nparam,ndata,iseed,iaccept,iadapt
      real*8 fnewLLK, curLLK, curMin
      real*8 curpars(nparam),savepar(nparam)
      real*8 parBest(nparam)
      real*8 rtn(ndata), rtnBest(ndata)
      real*8 rnd,ran1,diff_LLK
      external ran1
      logical accept

      diff_LLK = fnewLLK - curLLK
            
      accept = .false.
 
      if (diff_LLK .le. 0.0d0) Then
         accept = .true.
      else
         rnd = ran1(iseed)
         if ((exp(-diff_LLK)) .gt.rnd) Then
            accept = .true.
         else
            accept = .false.
         endif
      endif
               
! if step is accepted 

      if (accept) Then
         iaccept = iaccept + 1
         iadapt = iadapt + 1 
         savepar = curpars
         curLLK = fnewLLK
 
         if (curLLK .le. curMin) then !keep track of the best profile we have
            parBest = curpars
            curMin = curLLK
            rtnBest = rtn
         endif
!     if step is rejected - restore saved values of parameters
      else
         
         curpars = savepar
              
      endif

      return
      end subroutine MH1D

! ------------------------------------------------------------------------

C The Main Differential Equations- SIR model

      SUBROUTINE derivSIR(pars, Pop, dPop, np, nc)

      implicit none

      integer np, nc
      REAL*8 pars(np)
      REAL*8 Pop(nc), dPop(nc)

C
C     The SIR differential equations
C

C     dS/dt =
      dPop(1) = - pars(3) * Pop(1) * Pop(2)/pars(1)
C     dI/dt =
      dPop(2) =   pars(3) * Pop(1) * Pop(2)/pars(1)
     $          - Pop(2)/pars(2)
C     dR/dt =
      dPop(3) =   Pop(2)/pars(2)
c cummulative dsdt
      dPop(4) =  - dPop(1)

      RETURN
      END subroutine derivSIR

!
! -------------------------------------------------------------------
!

C The Main Differential Equations.

      SUBROUTINE derivSEIR(Param, Pop, dPop, np, nc)

      implicit real*8(a-h,o-z)

      integer np, nc
      REAL*8 Param(np), Pop(nc), dPop(nc)

C
C     The differential equations-including source/sink term
C

C     dS/dt =
      dPop(1) = - Param(3)*Pop(1)*Pop(3)/Param(1)

C     dE/dt
      dPop(2) = Param(3)*Pop(1)*Pop(3)/Param(1) - param(4) * Pop(2)
C     dI/dt =
      dPop(3) = param(4) * Pop(2) - Pop(3)/Param(2)
C     dR/dt =
      dPop(4) = Pop(3)/Param(2)
c cummulative sigma * E
      dPop(5) = Param(4)*Pop(2)
      RETURN
      END SUBROUTINE derivSEIR

C--------------------------------------------------------------------------------

      subroutine BuildIMID(ndata,nstep,tps,imid)

      implicit none

      integer ndata,nstep
      integer imid((ndata+1))
      real*8 tps(ndata)
      integer i

! Build the vector needed for calculating the incidence
! days is needed only to calculate imid
! tps is the 'cumulative' day number
! This works correctly provided tps holds cumulative days numbers
! that start with the first month or week and not some arbitrary start

! This line works for the dengue monthly/weekly data but not for the cdc data
! For now we have this somewhat dirty solution


      imid(1) = int(tps(1)) * int(nstep * 0.50)

      do i = 1, ndata
         imid((i+1)) = imid(1) + int(tps(i) * nstep)
      enddo

      return
      end subroutine BuildIMID

C--------------------------------------------------------------------------------
