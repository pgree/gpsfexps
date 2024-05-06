c
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning
c        of the xk-polynomial expansion routines
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains a set of subroutines for the handling
c        of xk polynomial expansions. The following is a brief 
c        description of these subroutines.
c 
c   xk_exps - constructs xk polynomial nodes, and corresponding Gaussian
c         weights. Also constructs the matrix v converting the
c         coefficients of an xk polynomial expansion into its values at
c         the n Gaussian nodes, and its inverse u, converting the
c         values of a function at n Gaussian nodes into the
c         coefficients of the corresponding xk polynomial series.
c 
c   xk_roots_find - finds the n roots of the n-degree xk polynomial
c         using an O(n**2) complexity algorithm.
c
c   xk_big_root_find - finds the biggest root (on the interval [0,1])
c         of the n-degree xk polynomial. 
c
c   xk_make_mat - given nodes, constructs the n*n matrix that maps 
c         a vector of coefficients of an xk polynomial expansion 
c         to the value of the expansion at the nodes, and the inverse
c         of that matrix, which takes values at nodes returns 
c         coefficients of an expansion
c
c   xk_weights - constructes the n Gaussian weights given n Gaussian
c         nodes.
c 
c   xk_normd_pol - evaluates the normalized xk polynomial
c
c   xk_normd_pols - evaluates normalized xk polynomials of degrees 0,...,n
c
c   xk_normd_pol_deriv - evaluates the normalized xk polynomial and 
c         its derivative
c
c   xk_normd_pol_derivs - evaluates the normalized xk polynomials and 
c         their derivatives of degrees 0,...,n
c
c   xk_normd_pol_deriv2 - evaluates the normalized xk polynomial in 
c         addition to its derivative and second derivative
c
c
c
c
        subroutine xk_exps(itype,n,k,work,rs,u,v,ws)
        implicit real*8 (a-h,o-z)
        real*8 rs(*),ws(*),work(*)
        real*8 v(n,n),u(n,n)
c
c         This subroutine constructs the n nodes (xs) of the n-degree
c         xk-polynomial on the interval [0,1] and the weights (ws)
c         for the gaussian quadrature. It also constructs
c         the matrix v, which converts the coefficients
c         of an xk-polynomial expansion into its values at the n
c         gaussian nodes, and its inverse u, which converts the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding xk-polynomial expansion.
c
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c  k - the parameter of xk-polynomials corresponding to the degree
c         of the weight function. Specifically, xk-polynomials are 
c         orthogonal on the interval [0,1] with respect to weight 
c         function x**k
c
c  work - work array that must be allocated at least 9*n+10 real*8 
c         locations long
c
c                 output parameters:
c 
c  rs - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the values of a function 
c         at the n gaussian nodes into the coefficients of its
c         xk-polynomial expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term xk-polynomial expansion into its values at
c         n gaussian nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  ws - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c

        icoefs1=21
        lcoefs1=4*n+10

        icoefs2=icoefs1+lcoefs1
        lcoefs2=4*n+10

        irest=icoefs2+lcoefs2

        work(1)=icoefs1+0.1
        work(2)=icoefs2+0.1

        work(20)=irest+0.1

        call xk_exps0(itype,n,k,work,rs,u,v,ws)

        return
        end
c
c
c
c
c
        subroutine xk_exps0(itype,n,k,work,rs,u,v,ws)
        implicit real*8 (a-h,o-z)
        real*8 rs(1),ws(1),work(1)
        real*8 v(n,n),u(n,n)
c
c         This subroutine constructs the n nodes (xs) of the n-degree
c         xk-polynomial on the interval [0,1] and the weights (ws)
c         for the gaussian quadrature. It also constructs
c         the matrix v, which converts the coefficients
c         of an xk-polynomial expansion into its values at the n
c         gaussian nodes, and its inverse u, which converts the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding xk-polynomial expansion.
c
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c  k - the parameter of xk-polynomials corresponding to the degree
c         of the weight function. Specifically, xk-polynomials are 
c         orthogonal on the interval [0,1] with respect to weight 
c         function x**k
c
c  work - work array that must be allocated at least 9*n real*8 
c         allocations long
c
c                 output parameters:
c 
c  rs - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the values of a function 
c         at the n gaussian nodes into the coefficients of its
c         xk-polynomial expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term xk-polynomial expansion into its values at
c         n gaussian nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  ws - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c

c
c         first check itype to see what the user wants computed
c

        if (itype .eq. 2) goto 199

        if (itype .eq. 0) then
          call xk_roots_find(n,k,work,rs)
          return
        endif

        if (itype .eq. 1) then
          call xk_roots_find(n,k,work,rs)
          call xk_make_weights(rs,n,k,work,ws)
          return
        endif

 199    continue


c
c         user wants nodes, weights, u and v computed
c
        call xk_roots_find(n,k,work,rs)
        call xk_make_mat(rs,n,k,work,v)

ccc        call prinf('k*',k,1)
ccc        call prin2('xs*',rs,n)
ccc        call prin2('a*',a,n*n)

c
c        . . . make weights
c
        
        do 100 i=1,n
ccc        call prinf('i*',i,1)
        wi=0

        do 90 j=1,n
        wi=wi+v(i,j)**2
 90     continue

        ws(i)=1/wi

 100    continue

ccc        call prin2('weights*',ws,n)

c
c        . . . and make inverse of a
c

        do 1000 i=1,n
        do 900 j=1,n
        u(j,i)=v(i,j)*ws(i)
 900    continue
 1000   continue


        return
        end
c
c
c
c
c
        subroutine xk_make_weights(rs,n,k,work,ws)
        implicit real*8 (a-h,o-z)
        real*8 rs(1),ws(1),work(1)

c
c         this subroutine computes the n weights used for gaussian 
c         quadrature of xk-polynomials. The nodes, xs, are 
c         an input
c

        ifis=work(20)
        lfis=n+10
        
        call xk_make_weights0(rs,n,k,work,work(ifis),ws)

        return
        end
c
c
c
c
c
        subroutine xk_make_weights0(rs,n,k,work,fis,ws)
        implicit real*8 (a-h,o-z)
        real*8 rs(1),fis(1),ws(1),work(1)

c
c         this subroutine computes the n weights used for gaussian 
c         quadrature of xk-polynomials. The nodes, xs, are 
c         an input
c

        do 100 i=1,n
ccc        call prinf('i*',i,1)
        call xk_normd_pols(rs(i),n,k,work,fis)
        wi=0

        do 90 j=1,n
        wi=wi+fis(j)**2
 90     continue

        ws(i)=1/wi

 100    continue

        return
        end
c
c
c
c
c
        subroutine xk_theta_ode_panel(theta0,thetat,n,k,nsteps,r0,rt)
        implicit real*8 (a-h,o-z)
        real*8 rs(10000),thetas(10000)

c
c        on the interval (theta0,thetat) use 2nd order Runge Kutta
c        to solve the ODE,
c          dx/dtheta=(-sqrt(b)-(b1/(4*b)+a/2.0d0)*sin(2*theta))**(-1)
c        which is the ODE of the Prufer transform. Here,
c        theta0>thetat and we march in the negative direction
c

c
c        initialize
c

        done=1
        pi=atan(done)*4

ccc        theta0=pi/2.0d0
ccc        thetat=-pi/2

        h=(thetat-theta0)/(nsteps-1)

        thetas(1)=theta0

        do 100 i=1,nsteps
        thetas(i+1)=thetas(i)+h
 100    continue

ccc        call prin2('thetas*',thetas,nsteps)
ccc        call prin2('r0*',r0,1)

c
c        now march
c
        rs(1)=r0

ccc        call prin2('r0*',1-r0,1)

        do 200 i=1,nsteps-1
ccc        call prinf('i*',i,1)

        done=1
        x=1-2*rs(i)
        c1=(1-x**2)/4.0d0
        c2=(k+(k+2)*x)/2.0d0
        c3=n*(n+k+done)

        a=c2/c1
        b=c3/c1
        b1=b*(-x/c1)

        theta=thetas(i)

ccc        call prin2('rs(i)*',rs(i),1)
ccc        call prin2('b*',b,1)
ccc        call prin2('term 1*',-sqrt(b),1)
ccc        call prin2('term 2*',-(b1/(4*b)+a/2.0d0)*sin(2*theta),1)

        der=-sqrt(b)-(b1/(4.0d0*b)+a/2.0d0)*sin(2*theta)
        der=der**(-1)

ccc        call prin2('der*',der,1)
ccc        call prin2('h*',h,1)
ccc        call prin2('der x h*',der*h,1)

        rs(i+1)=rs(i)+h*der
ccc        call prin2('new xs(i+1)*',1-xs(i+1),1)

c
c        tighten via rk2
c

        xp=1-2*rs(i+1)
        c1p=(1-xp**2)/4.0d0
        c2p=-(-k-(k+2)*xp)/2.0d0
        c3p=n*(n+k+done)

        ap=c2p/c1p
        bp=c3p/c1p
        b1p=bp*(-xp/c1p)

        thetap=thetas(i+1)

ccc        call prin2('bp*',bp,1)
ccc        call prin2('rs(i+1)*',rs(i+1),1)

        derp=-sqrt(bp)-(b1p/(4*bp)+ap/2.0d0)*sin(2*thetap)
        derp=derp**(-1)

        rs(i+1)=rs(i)+h/2.0d0*(der+derp)

ccc        call prin2('rs(i+1) after RK*',rs(i+1),1)
ccc        call prin2('b*',b,1)

 200    continue

c
c        set final point as rt for output
c

        rt=rs(nsteps)


        return
        end
c
c
c
c
c
        subroutine xk_big_root_find(n,k,nsteps,work,val)
        implicit real*8 (a-h,o-z)
        real*8 work(1)

c
c        This subroutine finds the biggest root of the n-degree 
c        xk-polynomial by first picking a starting point x0 
c        that depends on n and solving dtheta/dx on a panel equal in
c        to its distance to 1. After approximation is found, we 
c        then do newton to tighten the approximation.
c
        
        done=1
        pi=atan(done)*4

c
c         find initial guess, x0 which should be greater than the
c         solution, the biggest root on (0,1)

        if (n .lt. 1.0d3) then
        r0=1-1.0d-6
        goto 9000
        endif

        if (n .lt. 1.0d4) then
        r0=1-1.0d-8
        goto 9000
        endif

        if (n .lt. 1.0d5) then
        r0=1-1.0d-10
        goto 9000
        endif

        r0=1-1.0d-12

 9000   continue

        call xk_theta_eval(r0,n,k,theta0)
        thetat=pi/2

ccc        call prinf('n*',n,1)
ccc        call prin2('theta0*',theta0,1)

c
c        march to first root by solving prufer ode, that is, 
c        find xt such that thetat=pi/2
c

        nsteps=100
        call xk_theta_ode_panel(theta0,thetat,n,k,nsteps,r0,rt)

c
c        now we tighten via newton where our initial guess is xt
c

ccc        call prin2('rt*',rt,1)
ccc        call prin2('f at xt*',f0,1)
ccc        call prin2('initial guess*',rt,1)

        do 200 i=1,50

        call xk_normd_pol_deriv(rt,n,k,work,f0,f1)
        rt=rt-f0/f1

ccc        call prin2('rt*',rt,1)
ccc        call prin2('f0*',f0,1)
ccc        call prin2('f0*',f0/f1,1)

c
c        when we notice some convergence do 5 more steps
c

        if (f0/f1 .lt. 1.0d-8) then
          do 120 ijk=1,5
          call xk_normd_pol_deriv(rt,n,k,work,f0,f1)
          rt=rt-f0/f1
 120      continue
          goto 250
        endif

 200    continue
 250    continue

ccc        call prinf('i*',i,1)

        val=rt

c
c        compare to the true last root
c

ccc        call xk_roots_dumb(n,k,droots)
ccc        dd=droots(n)-val
ccc        call prin2('droots(n)*',droots(n),1)
ccc        call prin2('new root*',val,1)
ccc        call prin2('dd*',dd,1)

c
c        evaluate xk pol at val to check that it is, in fact, a root
c

ccc        call xk_normd_pol(val,n,k,f0)
ccc        call prin2('f at root*',f0,1)

        return
        end
c
c
c
c
c
        subroutine xk_roots_find(n,k,work,vals)
        implicit real*8 (a-h,o-z)
        real*8 vals(1),work(1)

c
c        This subroutine finds the n roots of the n-degree xk-polynomial 
c        by first using Prufer Transform to find initial guesses
c        then using newton to tighten
c

c
c        find the biggest root for initial condition
c
        nsteps=100
        call xk_big_root_find(n,k,nsteps,work,r0)
ccc        vals(1)=r0
        vals(n)=r0

        done=1
        pi=atan(done)*4
        theta0=pi/2

c
c        march from in the negative direction to get other roots
c
        nsteps=20

        do 100 i=1,n-1

        thetat=i*pi+pi/2
        theta0=thetat-pi

ccc        call prin2('theta0*',theta0,1)
ccc        call prin2('thetat*',thetat,1)
ccc        call prin2('r0*',r0,1)
ccc        call prin2('rt before*',rt,1)

        call xk_theta_ode_panel(theta0,thetat,n,k,nsteps,r0,rt)

ccc        call prin2('rt after*',rt,1)
ccc        vals(i+1)=rt

        vals(n-i)=rt
        r0=rt

 100    continue

c
c        tighten each root via newton
c

        do 150 i=1,n
        do 125 j=1,50

        call xk_normd_pol_deriv(vals(i),n,k,work,f0,f1)

        vals(i)=vals(i)-f0/f1

ccc        call prinf('j*',j,1)
ccc        call prin2('f at new root*',f0,1)
ccc        call prin2('f1 at new root*',f1,1)
ccc        call prin2('new root*',vals(i),1)
ccc        call prin2('f0/f1*',f0/f1,1)

c
c        when we start to see some convergence, take 5 more steps
c
        if (f0/f1 .lt. 1.0d-8) then
          do 120 ijk=1,5
          call xk_normd_pol_deriv(vals(i),n,k,work,f0,f1)
          vals(i)=vals(i)-f0/f1
 120      continue
          goto 150
        endif

 125    continue
 150    continue

        return
        end
c
c
c
c
c
        subroutine xk_jac_pols(x,n,a,a1s,b2s,b3s,b4s,vals)
        implicit real *8 (a-h,o-z)
        real*8 vals(1)
        real*8 a1s(1),b2s(1),b3s(1),b4s(1)
        save
        data ifinit /0/
        data nold /-10/
c
c        This subroutine computes the Jacobi polynomials 
c
c      P^(a,0)_0(x), P^(a,0)_1(x), ... , P^(a,0)_n(x),
c
c        and stores the result in vals.
c

ccc        if (n .eq .nold) goto 75

ccc        if (ifinit .eq. 0) then
        do 50 i=1,n-1
        a1s(i)=2d0*(i+1)*(i+a+1)*(2*i+a)
        b2s(i)=(2d0*i+a+1)*a**2/a1s(i)
        b3s(i)=(2d0*i+a)*(2d0*i+a+1)*(2d0*i+a+2)/a1s(i)
        b4s(i)=2d0*(i+a)*i*(2*i+a+2)/a1s(i)

 50     continue

ccc        call prinf('entered loop, ifinit=*',ifinit,1)
ccc        call prinf('n*',n,1)

        ifinit=1

ccc        endif

        nold=n

 75     continue


ccc        call prin2('b4s*',b4s,12)

        vals(1) = 1
        vals(2) = 0.5d0*(a+(a+2)*x)

ccc        call prinf('n*',n,1)
ccc        call prin2('b4s*',b4s,12)
        do 200 i=1,n-1
        
ccc        vals2(i+2)=((a2s(i)+a3s(i)*x)*vals(i+1)-a4s(i)*vals(i))/a1s(i)
ccc        call prin2('dd*',vals(i+2)-vals2(i+2),1)

        vals(i+2)=(b2s(i)+b3s(i)*x)*vals(i+1)-b4s(i)*vals(i)

200     continue

        return
        end
c
c
c
c
c
        subroutine xk_jac_pol_deriv(x,n,a,a1s,b2s,b3s,b4s,val0,val1)
        implicit real *8 (a-h,o-z)
        save
        real*8 a1s(1),b2s(1),b3s(1),b4s(1)
        data ifinit /0/
        data nold /-10/
c
c         This subroutine computes the derivative of the 
c         Jacobi polynomial P^(a,0)_n(x).
c

c
c         first compute coefficients for recurrence
c

ccc        if (n .eq. nold) goto 75
ccc        if (ifinit .eq. 0) then

        do 50 i=1,n-1
        a1s(i)=2d0*(i+1)*(i+a+1)*(2*i+a)
        b2s(i)=(2d0*i+a+1)*a**2/a1s(i)
        b3s(i)=(2d0*i+a)*(2d0*i+a+1)*(2d0*i+a+2)/a1s(i)
        b4s(i)=2d0*(i+a)*i*(2*i+a+2)/a1s(i)
 50     continue

ccc        call prinf('entered loop, ifinit=*',ifinit,1)
        ifinit=1

ccc        endif

        nold=n

 75     continue


        if (n .eq. 0) then
          val1=0
          val0=1
          return
        endif

        if (n .eq. 1) then
          val0=0.5d0*(a+(a+2)*x)
          val1=0.5d0*(a+2)
          return
        endif

c
c        . . . do jacobi polynomial recurrence
c
        f0 = 1
        f1 = 0.5d0*(a+(a+2)*x)

        g0 = 0
        g1 = 0.5d0*(a+2)

        do 200 i=1,n-1

ccc            f2=((a2s(i)+a3s(i)*x)*f1 - a4s(i)*f0)/a1s(i)
ccc            g2=((a2s(i)+a3s(i)*x)*g1 +a3s(i)*f1- a4s(i)*g0)/a1s(i)

            f2=(b2s(i)+b3s(i)*x)*f1-b4s(i)*f0
            g2=(b2s(i)+b3s(i)*x)*g1+b3s(i)*f1-b4s(i)*g0

            f0=f1
            f1=f2

            g0=g1
            g1=g2

200     continue

        val1=g2
        val0=f2
        
        return
        end
c
c
c
c
c
        subroutine xk_theta_eval(r,n,k,val)
        implicit real*8 (a-h,o-z)

c
c         evaluate the function theta of the Prufer Transform
c         corresponding to the n-degree xk-polynomial
c

ccc        call prin2('r in theta_eval*',r,1)

        call xk_normd_pol_deriv2(r,n,k,f0,f1,f2)

        done=1

        x=1-2*r
ccc        call prin2('r*',r,1)
ccc        call prin2('x*',x,1)
        c1=(1-x**2)/4.0d0
        c2=-(-k-(k+2)*x)/2.0d0
        c3=n*(n+k+done)

        a=c2/c1
        b=c3/c1
        b1=b*(-x/c1)

        z=f1/f0
ccc        call prin2('c1*',c1,1)
ccc        call prin2('c3*',c3,1)
ccc        call prin2('b*',b,1)
        gam=sqrt(b)
        val=atan(z/gam)
 
ccc        call prin2('theta in theta_eval*',val,1)
ccc        call prin2('x*',x,1)
ccc        call prin2('gam*',gam,1)
ccc        call prin2('z*',z,1)
       

        return
        end
c
c
c
c
c
        subroutine xk_make_mat(xs,n,k,work,a)
        implicit real*8 (a-h,o-z)
        real*8 xs(1),a(n,n),work(1)

c
c         make the matrix a in which
c         a_{i,j}=f_j(x_i)
c         and make its inverse, ai. a takes a vector of coefficients
c         of an expansion in xk polynomials and returns a vector of 
c         tabulated values at the xis for that expansion
c

        ifis=work(20)
        lfis=n+10

        call xk_make_mat0(xs,n,k,work,work(ifis),a)

        return
        end
c
c
c
c
c
        subroutine xk_make_mat0(xs,n,k,work,fis,a)
        implicit real*8 (a-h,o-z)
        real*8 xs(1),a(n,n),work(1),fis(1)

c
c         make the matrix a in which
c         a_{i,j}=f_j(x_i)
c         and make its inverse, ai. a takes a vector of coefficients
c         of an expansion in xk polynomials and returns a vector of 
c         tabulated values at the xis for that expansion
c

        
        do 100 i=1,n

        call xk_normd_pols(xs(i),n,k,work,fis)

ccc        call prin2('fis*',fis,n)

        do 90 j=1,n
        a(i,j)=fis(j)
 90     continue

 100    continue

ccc        call prin2('a*',a,n*n)
ccc        call orthom(ai,n,work,cond)

        return
        end
c
c
c
c
c
        subroutine xk_normd_pol_deriv2(x,n,k,val0,val1,val2)
        implicit real *8 (a-h,o-z)
c
c        This subroutine evaluatues the normalized 
c        xk polynomial along with its first and second
c        derivatives. They are stored in val0, val1, and val2
c        respectively
c
        
        r=1-2*x
        dk=k
        call xk_jac_pol_deriv2(r,n,dk,val0,val1,val2)

ccc        call prin2('val2*',val2,1)
ccc        call prin2('r*',r,1)

        val2=4*val2*(dk+2*n+1)**(1/2.0d0)
        val1=-2*val1*(dk+2*n+1)**(1/2.0d0)
        val0=val0*(dk+2*n+1)**(1/2.0d0)


        return
        end
c
c
c
c
c
        subroutine xk_jac_pol_deriv2(x,n,a,val0,val1,val2)
        implicit real *8 (a-h,o-z)
c
c        This subroutine computes the Jacobi polynomial 
c        P^(a,0)_n(x), along with its first and 
c        second derivatives and stores them in 
c        val0,val1, and val2 respectively
c
        if (n .eq. 0) then
          val0=1
          val1=0
          val2=0
          return
        endif

        if (n .eq. 1) then
          val0=0.5d0*(a+(a+2)*x)
          val1=0.5d0*(a+2)
          val2=0
          return
        endif

c
c        . . . do jacobi polynomial and derivatives recurrence
c
        f0 = 1
        f1 = 0.5d0*(a+(a+2)*x)

        g0 = 0
        g1 = 0.5d0*(a+2)

        h0=0
        h1=0

        do 200 i=1,n-1

            a1i = 2d0*(i+1)*(i+a+1)*(2*i+a)
            a2i = (2d0*i+a+1)*a**2
            a3i = (2d0*i+a)*(2d0*i+a+1)*(2d0*i+a+2)
            a4i = 2d0*(i+a)*i*(2*i+a+2)

ccc            call prin2('a1i=*',a1i,1)
ccc            call prin2('a2i=*',a2i,1)
ccc            call prin2('a3i=*',a3i,1)
ccc            call prin2('a4i=*',a4i,1)

            f2 = ((a2i+a3i*x)*f1 - a4i*f0)/a1i

            g2 = ((a2i+a3i*x)*g1 +a3i*f1- a4i*g0)/a1i

            h2=  ((a2i+a3i*x)*h1 +2*a3i*g1- a4i*h0)/a1i


            f0=f1
            f1=f2

            g0=g1
            g1=g2

            h0=h1
            h1=h2

ccc            call prin2('h2*',h2,1)

200     continue

        val2=h2
        val1=g2
        val0=f2

        return
        end
c
c
c
c
c
        subroutine xk_normd_pols(r,n,k,work,vals)
        implicit real *8 (a-h,o-z)
        dimension vals(1),work(1)

c
c        This subroutine evaluatues the xk-polynomials
c        of degrees 0,1,...,n

c
c        . . . but first allocate memory
c

        ia1s=work(1)
        la1s=n+10

        ib2s=ia1s+la1s
        lb2s=n+10

        ib3s=ib2s+lb2s
        lb3s=n+10

        ib4s=ib3s+lb3s
        lb4s=n+10


        call xk_normd_pols0(r,n,k,work(ia1s),work(ib2s),work(ib3s),
     1     work(ib4s),vals)

        return
        end
c
c
c
c
c
        subroutine xk_normd_pols0(r,n,k,a1s,b2s,b3s,b4s,vals)
        implicit real *8 (a-h,o-z)
        dimension vals(1),a1s(1),b2s(1),b3s(1),b4s(1)

c
c        This subroutine evaluatues the xk-polynomials
c        of degrees 0,1,...,n

        dk=k
        x=1-2*r

        call xk_jac_pols(x,n,dk,a1s,b2s,b3s,b4s,vals)


        do 1200 i=1,n+1
        vals(i)=(dk+2*(i-1)+1)**(1/2.0d0)*vals(i)
 1200 continue

ccc        call prin2('work*',work,4*n)

        return
        end
c
c
c
c
c
        subroutine xk_normd_pol_derivs(x,n,k,vals0,vals1)
        implicit real *8 (a-h,o-z)
        dimension vals1(1),vals0(1)
c
c        This subroutine evaluatues normalized xk polynomials
c        and their derivatives for n=0,1,2,..., where 0<r<1.
c

        r=1-2*x
        dk=k
        call xk_jac_pol_derivs(r,n,dk,vals0,vals1)

        do 1200 i=1,n+1

        vals1(i)=-2*vals1(i)*(dk+2*(i-1)+1)**(1/2.0d0)
        vals0(i)=vals0(i)*(dk+2*(i-1)+1)**(1/2.0d0)

 1200 continue


        return
        end
c
c
c
c
c
        subroutine xk_normd_pol_deriv(x,n,k,work,val0,val1)
        implicit real *8 (a-h,o-z)
        real*8 work(1)
c
c        This subroutine evaluatues the n-degree normalized 
c        xk polynomial and its derivative. The values are stored 
c        in val0 and val1 respectively
c
        ia1s=work(2)
        la1s=n+10

        ib2s=ia1s+la1s
        lb2s=n+10

        ib3s=ib2s+lb2s
        lb3s=n+10

        ib4s=ib3s+lb3s
        lb4s=n+10

        
        call xk_normd_pol_deriv0(x,n,k,work(ia1s),work(ib2s),
     1     work(ib3s),work(ib4s),val0,val1)


        return
        end
c
c
c
c
c
        subroutine xk_normd_pol_deriv0(x,n,k,a1s,b2s,b3s,b4s,
     1     val0,val1)
        implicit real *8 (a-h,o-z)
        real*8 a1s(1),b2s(1),b3s(1),b4s(1)
c
c        This subroutine evaluatues the n-degree normalized 
c        xk polynomial and its derivative. The values are stored 
c        in val0 and val1 respectively
c
        
        r=1-2*x
        dk=k

        call xk_jac_pol_deriv(r,n,dk,a1s,b2s,b3s,b4s,val0,val1)

ccc        call prin2('val2*',val2,1)
ccc        call prin2('r*',r,1)

        val1=-2*val1*(dk+2*n+1)**(1/2.0d0)
        val0=val0*(dk+2*n+1)**(1/2.0d0)

        return
        end
c
c
c
c
c
        subroutine xk_normd_pol(x,n,k,val)
        implicit real *8 (a-h,o-z)

c
c         compute n-degree xk polynomial,
c            f_n^k(x)=(k+2n+1)**(1/2)*P_n^(k,0)(1-2x),
c         where P_n^(k,0) is a Jacobi polynomial.
c         These polynomials have L2 norm 1 with weight function
c         x**k

        r=1-2*x
        dk=k

ccc        call prinf('k*',k,1)
ccc        call prin2('dk*',dk,1)

        call xk_jac_pol(r,n,dk,val)

        val=val*(dk+2*n+1)**(1/2.0d0)

ccc        call prin2('val*',val,1)

        return
        end
c
c
c
c
c
        subroutine xk_jac_pol_derivs(x,n,a,vals0,vals1)
        implicit real *8 (a-h,o-z)
        real*8 vals0(1),vals1(1)
c
c        This subroutine computes the derivatives of Jacobi polynomials 
c
c      P'^(a,0)_0(x), P'^(a,0)_1(x), ... , P'^(a,0)_n(x),
c
c        and stores the result in vals.
c
        vals0(1)=1
        vals0(2)=0.5d0*(a+(a+2)*x)

        vals1(1)=0
        vals1(2)=0.5d0*(a+2)

c
c        . . . do jacobi polynomial recurrence
c
        f0 = 1
        f1 = 0.5d0*(a+(a+2)*x)

        g0 = 0
        g1 = 0.5d0*(a+2)

        do 200 i=1,n-1

            a1i = 2d0*(i+1)*(i+a+1)*(2*i+a)
            a2i = (2d0*i+a+1)*a**2
            a3i = (2d0*i+a)*(2d0*i+a+1)*(2d0*i+a+2)
            a4i = 2d0*(i+a)*i*(2*i+a+2)

ccc            call prin2('a1i=*',a1i,1)
ccc            call prin2('a2i=*',a2i,1)
ccc            call prin2('a3i=*',a3i,1)
ccc            call prin2('a4i=*',a4i,1)

            f2 = ((a2i+a3i*x)*f1 - a4i*f0)/a1i
            vals0(i+2)=f2

            g2 = ((a2i+a3i*x)*g1 +a3i*f1- a4i*g0)/a1i
            vals1(i+2)=g2

            f0=f1
            f1=f2

            g0=g1
            g1=g2

200     continue

        return
        end
c
c
c
c
c
        subroutine xk_matvec(a,x,n,y)
        implicit real*8 (a-h,o-z)
        real*8 a(n,n),x(n),y(n)


        do 100 i=1,n
        y(i)=0
        do 90 j=1,n
        y(i)=y(i)+x(j)*a(i,j)
 90     continue
 100    continue


        return
        end
c
c
c
c
c
        subroutine xk_matmult(a,b,n,c)
        implicit real*8 (a-h,o-z)
        real*8 a(n,n),b(n,n),c(n,n)


        do 100 i=1,n
        do 95 j=1,n
        c(i,j)=0
        do 90 k=1,n
        c(i,j)=c(i,j)+a(i,k)*b(k,j)
 90     continue
 95     continue
 100    continue


        return
        end
c
c
c
c
c
        subroutine xk_jac_pol(x,n,a,val)
        implicit real *8 (a-h,o-z)
c
c        This subroutine computes the Jacobi polynomial P^(a,0)_n(x).
c
c

        if (n .eq. 0) then
          val=1
          return
        endif

        if (n .eq. 1) then
          val=0.5d0*(a+(a+2)*x)
          return
        endif


        f0 = 1
        f1 = 0.5d0*(a+(a+2)*x)

        do 200 i=1,n-1

            a1i = 2d0*(i+1)*(i+a+1)*(2*i+a)
            a2i = (2d0*i+a+1)*a**2
            a3i = (2d0*i+a)*(2d0*i+a+1)*(2d0*i+a+2)
            a4i = 2d0*(i+a)*i*(2*i+a+2)

ccc            call prin2('a1i=*',a1i,1)
ccc            call prin2('a2i=*',a2i,1)
ccc            call prin2('a3i=*',a3i,1)
ccc            call prin2('a4i=*',a4i,1)

            f2 = ((a2i+a3i*x)*f1 - a4i*f0)/a1i

            f0=f1
            f1=f2
200     continue

        val = f2

        return
        end
