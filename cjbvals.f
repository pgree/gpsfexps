(c) Vladimir Rokhlin
c
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c         This is the end of the debugging code, and the beginning of the 
c         Bessel function code proper
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file has three user-callable subroutines: cjbvals, rjbvals,
c        and rjbvals_half.  Following is a brief description of the
c        three subroutines
c
c   cjbvals - evaluates Bessel J-functions of the complex argument z. The 
c        number nvals of functions returned is determined (inter alia) by 
c        the user-specified precision eps
c
c   rjbvals - evaluates Bessel J-functions of the real argument z. The 
c        number nvals of functions returned is determined (inter alia) by 
c        the user-specified precision eps
c
c   rjbvals_half - evaluates Bessel J-functions of half-integer order 
c        of the real argument z. The number nvals of functions returned
c        is determined (inter alia) by the user-specified precision eps
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine rjbvals(z,eps,vals,nvals,nmax)
        implicit real *8 (a-h,o-z)
        save
        real *8 z,cv1,cv2,cv3,vals(1),cd,cd2,coef
c
c        This subroutine evaluates Bessel J-functions 
c
c       J_0(z), J_1(z), J_3(z), ... , J_nvals(z)
c
c        of the real argument z. The number nvals of functions returned
c        is determined (inter alia) by the user-specified precision eps.
c        More specifically, 
c
c        abs(J_{nvals} (z)) \sim eps,                               (1)
c
c        and the value J_{nvals} (z) is itself calculated to the
c        relative precision eps. Thus, this is a fairly gold-plated
c        code; no attempt to optimize it for the CPU time has been
c        made.
c
c                      Input parameters:
c
c  z - the point at which the Bessel functions are to be evaluated
c  eps - the relative precision to which the value J_{nvals} (z) is 
c        to be calculated
c  
c                      Output parameters:
c
c  vals - the array of Bessel functions (PLEASE NOTE THAT THIS ARRAY
c        MUST CONTAIN ELEMENT NUMBER 0)
c  nvals - the number of elements in array vals that are returned with 
c        relative precision eps; also, please note that vals(nvals)
c        will be roughly of (relative) size eps.
c  nmax - the total number of elements in array vals used by this 
c        subroutine; please note that it is somewhat greater than 
c        nvals
c
c
c       . . . find how high we should be going
c
        cv1=0
        cv2=1
c
        dlarge=1/eps**2/10
        do 1200 i=1,10 000 000 
c
        nmax=i
        cv3=2*i/z*cv2-cv1
c
        dd=abs(cv3)
        if(dd .gt. dlarge) goto 1400  
c
        cv1=cv2
        cv2=cv3
c
 1200 continue
 1400 continue
c
ccc        call prinf('nmax=*',nmax,1)
c
c        having found the maximum point, construct the 
c        unscaled Bessel functions
c
        vals(nmax)=0
        vals(nmax-1)=1

ccc        call prin2_long('vals(nmax)=*',vals(nmax),1)
ccc        call prin2_long('vals(nmax-1)=*',vals(nmax-1),1)

        do 1600 i=nmax-1,1,-1
c
        vals(i-1)=2*i/z*vals(i)-vals(i+1)
 1600 continue
c
ccc        call prin2('vals=*',vals(0),nmax)
ccc        call prinf('nmax=*',nmax,1)
ccc        call prin2_long('vals=*',vals(nmax-20),20+1)
ccc
ccc        call prin2_long('vals(nmax)=*',vals(nmax),1)
ccc        call prin2_long('vals(nmax-1)=*',vals(nmax-1),1)

c
c        evaluate the normalization coefficient
c
        dd1=abs(cos(z))
        dd2=abs(sin(z))
c
        rlarge=dd1+dd2
ccc        call prin2('z=*',z,1)

c
        if(dd2 .gt. dd1) goto 2000
        i0=0
        cd=vals(i0)/2
        d=-1
        do 1800 i=2,nmax,2
c
        cd=cd+vals(i)*d
c
        d=-d        
 1800 continue
c
        cd=cd*2
c
        coef=cos(z)/cd
        goto 2300
c
 2000 continue
c
        cd2=0
        d=1
        do 2200 i=1,nmax,2
c
        cd2=cd2+vals(i)*d
c
        d=-d        
 2200 continue
c
        cd2=cd2*2
c
        coef=sin(z)/cd2
c
 2300 continue

ccc        call prin2_long('coef=*',coef,1)

c
c       . . . scale them things
c
        do 2400 i=0,nmax-2
c
        vals(i)=vals(i)*coef
 2400 continue
c
c       determine the parameter nvals
c
c
        nvals=0
        eps22=eps**2*rlarge**2/10
c
        do 2600 i=nmax-2,1,-1
c
        dd=vals(i)**2
c
        if(dd .gt. eps22) goto 2800
c
        nvals=i
 2600 continue
 2800 continue
c
        return
        end
c
c 
c 
c
c
        subroutine rjbvals_half(z,eps,vals,nvals,nmax)
        implicit real *8 (a-h,o-z)
        save
        real *8 z,cv1,cv2,cv3,vals(1),cd,cd2,coef
c
c        This subroutine evaluates Bessel J-functions of half-integer
c        order
c
c       J_{1/2}(z), J_{3/2}(z), J_{5/2}(z), ... , J_{nvals+1/2}(z)
c
c        of the real argument z. The number nvals of functions returned
c        is determined (inter alia) by the user-specified precision eps.
c        More specifically, 
c
c        abs(J_{nvals+1/2} (z)) \sim eps,                               (1)
c
c        and the value J_{nvals+1/2} (z) is itself calculated to the
c        relative precision eps. Thus, this is a fairly gold-plated
c        code; no attempt to optimize it for the CPU time has been
c        made.
c
c                      Input parameters:
c
c  z - the point at which the Bessel functions are to be evaluated
c  eps - the relative precision to which the value J_{nvals+1/2} (z) is 
c        to be calculated
c  
c                      Output parameters:
c
c  vals - the array of Bessel functions (PLEASE NOTE THAT THIS ARRAY
c        MUST CONTAIN ELEMENT NUMBER 0), where vals(n) corresponds to
c        J_{n+1/2}(z).
c  nvals - the number of elements in array vals that are returned with 
c        relative precision eps; also, please note that vals(nvals)
c        will be roughly of (relative) size eps.
c  nmax - the total number of elements in array vals used by this 
c        subroutine; please note that it is somewhat greater than 
c        nvals
c
c       . . . find how high we should be going
c
        cv1=0
        cv2=1
c
        dlarge=1/eps**2/10
        do 1200 i=1,10 000 000 
c
        nmax=i
        cv3=2*(i+0.5d0)/z*cv2-cv1
c
        dd=abs(cv3)
        if(dd .gt. dlarge) goto 1400  
c
        cv1=cv2
        cv2=cv3
c
 1200 continue
 1400 continue
c
ccc        call prinf('nmax=*',nmax,1)
c
c        having found the maximum point, construct the 
c        unscaled Bessel functions
c
        vals(nmax)=0
        vals(nmax-1)=1
        do 1600 i=nmax-1,1,-1
c
        vals(i-1)=2*(i+0.5d0)/z*vals(i)-vals(i+1)
 1600 continue
c
ccc        call prin2('vals=*',vals(0),nmax)
c
c        evaluate the normalization coefficient
c
        pi = 3.1415926535897932384626433832795028841971

        v1 = sqrt(2*z/pi)*sin(z)/z
        v2 = sqrt(2*z/pi)*(sin(z)/z**2-cos(z)/z)

        dd1 = abs(v1)
        dd2 = abs(v2)

c
c       . . . scale them things
c

        if (dd2 .ge. dd1) then

            coef = v2/vals(1)
ccc            write(7,*) 'dd2 .ge. dd1'

        else

            coef = v1/vals(0)
ccc            write(7,*) 'dd1 .ge. dd2'

        end if

ccc        call prin2_long('coef=*',coef,1)


        do 2400 i=0,nmax-2
c
        vals(i)=vals(i)*coef
 2400 continue
c
c       determine the parameter nvals
c
c
        nvals=0
        eps22=eps**2/10
c
        do 2600 i=nmax-2,1,-1
c
        dd=vals(i)**2
c
        if(dd .gt. eps22) goto 2800
c
        nvals=i
 2600 continue
 2800 continue
c
        return
        end
c 
c 
        subroutine cjbvals(z,eps,vals,nvals,nmax)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,cv1,cv2,cv3,vals(1),cd,cd2,coef
c
c        This subroutine evaluates Bessel J-functions 
c
c       J_0(z), J_1(z), J_3(z), ... , J_nvals(z)
c
c        of the complex argument z. The number nvals of functions
c        returned is determined (inter alia) by the user-specified
c        precision eps. More specifically, 
c
c        abs(J_{nvals} (z)) \sim eps,                               (1)
c
c        and the value J_{nvals} (z) is itself calculated to the
c        relative precision eps. Thus, this is a fairly gold-plated
c        code; no attempt to optimize it for the CPU time has been
c        made.
c
c                      Input parameters:
c
c  z - the point at which the Bessel functions are to be evaluated
c  eps - the relative precision to which the value J_{nvals} (z) is 
c        to be calculated
c  
c                      Output parameters:
c
c  vals - the array of Bessel functions (PLEASE NOTE THAT THIS ARRAY
c        MUST CONTAIN ELEMENT NUMBER 0)
c  nvals - the number of elements in array vals that are returned with 
c        relative precision eps; also, please note that vals(nvals)
c        will be roughly of (relative) size eps.
c  nmax - the total number of elements in array vals used by this 
c        subroutine; please note that it is somewhat greater than 
c        nvals
c
c
c       . . . find how high we should be going
c
        cv1=0
        cv2=1
c
        dlarge=1/eps**2/10
        do 1200 i=1,10 000 000 
c
        nmax=i
        cv3=2*i/z*cv2-cv1
c
        dd=abs(cv3)
        if(dd .gt. dlarge) goto 1400  
c
        cv1=cv2
        cv2=cv3
c
 1200 continue
 1400 continue
c
ccc        call prinf('nmax=*',nmax,1)
c
c        having found the maximum point, construct the 
c        unscaled Bessel functions
c
        vals(nmax)=0
        vals(nmax-1)=1
        do 1600 i=nmax-1,1,-1
c
        vals(i-1)=2*i/z*vals(i)-vals(i+1)
 1600 continue
c
ccc        call prin2('vals=*',vals(0),nmax*2)
c
c        evaluate the normalization coefficient
c
        dd1=abs(cos(z))
        dd2=abs(sin(z))
c
        rlarge=dd1+dd2
ccc        call prin2('z=*',z,2)

c
        if(dd2 .gt. dd1) goto 2000
        i0=0
        cd=vals(i0)/2
        d=-1
        do 1800 i=2,nmax,2
c
        cd=cd+vals(i)*d
c
        d=-d        
 1800 continue
c
        cd=cd*2
c
        coef=cos(z)/cd
        goto 2300
c
 2000 continue
c
        cd2=0
        d=1
        do 2200 i=1,nmax,2
c
        cd2=cd2+vals(i)*d
c
        d=-d        
 2200 continue
c
        cd2=cd2*2
c
        coef=sin(z)/cd2
c
 2300 continue
c
c       . . . scale them things
c
        do 2400 i=0,nmax-2
c
        vals(i)=vals(i)*coef
 2400 continue
c
c       determine the parameter nvals
c
c
        nvals=0
        eps22=eps**2*rlarge**2/10
c
        do 2600 i=nmax-2,1,-1
c
        dd=vals(i)*conjg(vals(i))
c
        if(dd .gt. eps22) goto 2800
c
        nvals=i
 2600 continue
 2800 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine jfungen(x,rnu,f)
        implicit real *8 (a-h,o-z)
        save
        dimension gammas(100),fact(0:100)
c
c       evaluate the required gamma functions
c
        kmax=60
c
        call prin2('x=*',x,1)
        call prin2('rnu=*',rnu,1)
cccc        call gammanew_eval_extend(rnu,gam0)
        call gammanew_eval(rnu,gam0)
c
        d=rnu
        fact(0)=1
        gammas(1)=gam0
        do 1200 k=1,kmax
c
        gammas(k+1)=gammas(k)*d
        d=d+1
c
        fact(k)=fact(k-1)*k
 1200 continue
c
        call prin2('gammas as created*',gammas,60)
        call prin2('fact=*',fact,60)
c
c        evaluate the bessel function number rnu
c
        d=0
        do 1400 k=0,kmax-2
c
        d=d+(-x**2/4)**k /fact(k)/gammas(k+2)
 1400 continue
c
        f=d*(x/2)**rnu
c
        return
        end
c
c 
c 
c 
c 
        subroutine cjfungen(x,rnu,f)
        implicit real *8 (a-h,o-z)
        save
        dimension gammas(100),fact(0:100)
        complex *16 x,f,cd
c
c       evaluate the required gamma functions
c
        kmax=60
c
        call prin2('x=*',x,2)
        call prin2('rnu=*',rnu,1)
        call gammanew_eval_extend(rnu,gam0)
c
        d=rnu
        fact(0)=1
        gammas(1)=gam0
        do 1200 k=1,kmax
c
        gammas(k+1)=gammas(k)*d
        d=d+1
c
        fact(k)=fact(k-1)*k
 1200 continue
c
c        evaluate the bessel function number rnu
c
        cd=0
        do 1400 k=0,kmax-2
c
        cd=cd+(-x**2/4)**k /fact(k)/gammas(k+2)
 1400 continue
c
        f=cd*(x/2)**rnu
c
        return
        end
