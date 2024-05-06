        implicit real*8 (a-h,o-z) 
  
        call prini(6,13)
        call testit(n)

        stop
        end


        subroutine testit(n)
        implicit real*8 (a-h,o-z)
        real*8 w(2 000 000)

        call test_tbar_56(w)
ccc        call test_gammas_4(w)
ccc        call test_lnc_93(w)
ccc        call test_mnc_90(w)
ccc        call plot_phi_88(w)
ccc        call test_ode_268(w)

ccc        call phi_quad_test(w)
ccc        call exp_quad_test(w)

        stop
        end
c
c
c
c
c
        subroutine plot_phi_88(w)
        implicit real *8 (a-h,o-z)
        real*8 rs(100 000),phis(100 000),w(1),phis2(10 000), 
     1     vect(100 000),ws(100 000),zerns(100 000), 
     1     big_phis(10 000)
c
c        with gnuplot, plot the radial gpsfs, phis (88) and 
c        big_phis(87) of the preprint -- 
c          https://arxiv.org/pdf/1811.02733.pdf
c       
        neig=10
        nn=20
        ip=3
        c=200
        
c
c        compute eigenfunction expansion
c
        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)

c
c         lay down gauss legendre nodes on (0,1)
c        

        nnodes=400
        call legeexps(1,nnodes,rs,u,v,ws)
        a=0.0
        b=1.0
        do i=1,nnodes
        rs(i)=rs(i)/(2.0d0/(b-a))+(b+a)/2.0d0
        enddo

c
c         evaluate integrand at nodes
c        

        do i=1,nnodes
        call dprol_t_bar_coef2exp(vect,rs(i),nmat,nn,ip,w,phii)
        phis(i)=phii

        call dprol_normd_zern_exp_eval2(vect,rs(i),nn,ip,nmat,val0,
     1     val1,val2)
        big_phis(i)=val0

        call dprol_zern_normd_pol(rs(i),neig,nn,ip,val)
        zerns(i)=val
        enddo

c
c        print tikzpicture coordinates for latex
c
ccc        call print_gpsfs(nnodes, rs, phis)

c
c        . . . plot phis
c
        iw=15
        call quagraph(iw,rs,phis,nnodes,3,'title*')
        iw=16
        call quagraph(iw,rs,big_phis,nnodes,3,'title*')

        return
        end
c
c
c
c
c
        subroutine print_gpsfs(n, rs, fs)
        implicit real *8 (a-h,o-z)
        real*8 rs(n), fs(n)

        iu = 6
        open(iu)

 210    format('(',f8.5,','f8.5,')')

        do i=1,n
        write(iu, 210) rs(i), fs(i)/sqrt(rs(i))
        enddo

        return
        end
c
c
c
c
c
        subroutine phi_quad_test(w)
        implicit real*8 (a-h,o-z)
        real*8 w(1),droots(10 000),ws(10 000),rints(10 000)
        real*8 vect(10 000),fs(10 000),rs(10 000),djs(10 000)

c
c        check for ip=0, that a quadrature that integrates the gpsfs
c        
c          \Phi_{0,0}(r)*r,...,\Phi_{0,n-1}(r)*r
c
c        will also integrate J_0(c*r*rho)*r as long as we have enough
c        nodes. enough nodes means the eigenvalue of \Phi_{0,n} is small
c

        nn=2
        ip=2
        c=20
        rho=0.8

c
c        first use n nodes and then double to check accuracy
c
        n=20
        call dprol_bigphi_exps(n,nn,ip,c,w,rs,ws)

c
c        . . . check that eigenvalue is small
c
        neig=n+1
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam,vect)
        call prin2('dgam*',dgam,1)

c
c        tabulate function to be integrated
c

        eps=1.0d-75

        do i=1,n
        r=rs(i)
        call dprol_normd_zern_exp_eval(vect,r,nn,ip,nmat,val)
        fs(i)=val

        call rjbvals(c*r*rho,eps,djs,nvals,nmax)
        fs(i)=djs(1)*r
        enddo

c
c        . . . compute integral
c        
        rint=0
        do i=1,n
        rint=rint+fs(i)*ws(i)
        enddo

        call prin2('rint*',rint,1)

c
c       . . . double nodes to compare
c

        n=2*n
        call dprol_bigphi_exps(n,nn,ip,c,w,rs,ws)

c
c        tabulate function to be integrated
c

        do 600 i=1,n
        r=rs(i)
        call dprol_normd_zern_exp_eval(vect,r,nn,ip,nmat,val)
        fs(i)=val

        call rjbvals(c*r*rho,eps,djs,nvals,nmax)
        fs(i)=djs(1)*r
 600    continue

c
c        . . . compute integral
c        
        rint2=0
        do 800 i=1,n
        rint2=rint2+fs(i)*ws(i)
 800    continue

        call prin2('rint2*',rint2,1)

        dd=(rint2-rint)/rint2
        call prin2('dd*',dd,1)

        iw=19
        call quagraph(iw,rs,fs,n,3,'title*')

        return
        end
c
c
c
c
c
        subroutine test_gammas_4(w)
        implicit real *8 (a-h,o-z)
        real*8 w(1),vect(100 000),dgams(10 000)
        real*8 rlams(10 000),vects(100 000)

c
c        compute an eigenvalue of the integral operator m, 
c        in 2 ways. first, by computing the single eigenvalue
c        method of section 4.1 of the preprint -- 
c        https://arxiv.org/pdf/1811.02733.pdf -- then 
c        by computing all of the first n eigenvalues using
c        the method of section 4.2.
c

        nn=2
        ip=1
        c=10.0d0

c
c        compute the first neigs eigenvalues using ratios
c

        neigs=24
        call dprol_gammas_eval(neigs,nn,ip,c,w,nmat,rlams,dgams,vects)
        call prin2('dgams*',dgams,neigs)

        do i=1,neigs
        call dprol_gamma_n_eval(i,nn,ip,c,w,nmat,rlam,dgam1,vect)
        dgams(i) = dgam1
        enddo
        call prin2('dgams2*', dgams, neigs)

c
c        evaluate just the neig^th eigenvalue using a_0, the first
c        component of the eigenvector
c
        neig=neigs
        do i=1,neigs
        neig = i
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam1,vect)
        enddo
        call prin2('dgam1*',dgam1,1)

        call prin2('dgam2*',dgams(neig),1)
        dd=(dgams(neig)-dgam1)/dgam1
        call prin2('relative difference*',dd,1)

c
c        . . . evaluate the neig^th  eigenvalue in a primative way, by 
c        applying the integral operator using integration via legendre 
c        nodes and checking ratio at some point

ccc        call m_n_eigen_eval_rough(neig,nn,ip,c,w,dgam3,vect)
ccc        call prin2('dgam3*',dgam3,1)
ccc        call prin2('dd*',(dgams(neig)-dgam3)/dgam3,1)


        return
        end
c
c
c
c
c
        subroutine test_tbar_56(w)
        implicit real*8 (a-h,o-z)
        real*8 w(1), vect(10 000)

c
c        check equation (56) of the preprint 
c        https://arxiv.org/pdf/1811.02733.pdf. note that the 
c        equation in the original arxiv paper from 2018 is 
c        wrong. the correct equation was added to the next version.
c       
        n=4
        nn=2
        ip=1
        c=100
        r = 0.6

c
c        set variables used in formula
c
        a = nn + ip/2.0d0
        at = -2.0d0*n*(2*n + a + 2)
        bt = 2.0d0*a*(2*n + a + 1)
        ct = 2.0d0*(n + a + 1) * (2*n + a)
        a1 = n*(2.0d0*a + 4*n - 1)*(2*n + a + 2)
        b1 = a*(2.0d0*n + a + 1) - 2*(2*n+a)*(2*n+a+1)*(2*n+a+2)
        c1 = (2.0d0*a + 4*n + 5)*(n + a + 1)*(2*n+a)

c
c        evaluate tbar functions
c        
        nm = n - 1
        np = n + 1
        call dprol_tbar_deriv(r,n,nn,ip,t0,t1)
        call dprol_tbar_deriv(r,nm,nn,ip,t0m,t1m)
        call dprol_tbar_deriv(r,np,nn,ip,t0p,t1p)

c
c        scale tbar functions to be t functions
c
        const = sqrt(2.0d0*(2.0d0*n + a + 1.0d0))
        constm = sqrt(2.0d0*(2.0d0*nm + a + 1.0d0))
        constp = sqrt(2.0d0*(2.0d0*np + a + 1.0d0))
        t0 = t0 / const
        t0m = t0m / constm
        t0p = t0p / constp
        t1 = t1 / const
        t1m = t1m / constm
        t1p = t1p / constp

c
c        check formula
c
        dlhs = at*r*t1m - bt*r*t1 + ct*r*t1p 
        rhs = a1*t0m - b1*t0 + c1*t0p
        dd = dlhs - rhs

        call prin2('dlhs*', dlhs, 1)
        call prin2('rhs*', rhs, 1)
        call prin2('dd*', dd, 1)


        return
        end
c
c
c
c
c
        subroutine test_ode_268(w)
        implicit real*8 (a-h,o-z)
        real*8 w(1), vect(10 000)

c
c        check second order ode (268) of the original arxiv
c        paper https://arxiv.org/pdf/1811.02733.pdf
c       
        neig=20
        nn=1
        ip=2
        c=100
        
c
c        compute eigenfunction expansion of one eigenfunction
c        phi as well as the corresponding eigenvalue of the 
c        differential operator L_{N, c}
c
        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)

c
c         evaluate phi and its first two derivatives
c       
        r = 0.2
        nmat1 = nmat + 1
        call dprol_normd_zern_exp_eval2(vect,r,nn,ip,nmat1,val0,
     1     val1,val2)

c
c        check second order ode (268)
c
        t2 = (1 - r**2) * r**2 * val2

        t1 = (ip + 1.0d0) * r - (ip + 3.0d0)*r**3
        t1 = t1 * val1

        t0 = -rlam*r**2 - (ip+1)*(ip+3)/4.0d0*r**2
        t0 = t0 - 1.0d0*nn*(nn + ip) - c**2*r**4
        t0 = t0 * val0 

        dd = t0 + t1 + t2

        call prin2('dd*', dd, 1)

        return
        end
c
c
c
c
c
        subroutine test_mnc_90(w)
        implicit real*8 (a-h,o-z)
        real*8 rhos(10000),vals(10000), djs(10000),deigs(10000)
        real*8 ws(10000),dds(10000),vect(10 000),w(1)
        real*8 dphis(10000),dmphis(10000),dmphis_scaled(10000)

c
c        test that \phi is an eigenfunction of the integral
c        operator M_{N, c}. in particular, test formula (90)
c        of https://arxiv.org/1811.02733
c
        neig=5
        nn=1
        ip=2
        c=100

c
c        find zernike expansion of eigenfucntion phi
c
        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)

c
c        find corresponding eigenvalue of integral operator M_{N, c}
c
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam,vect)
        call prin2('dgam*',dgam,1)

c
c        check that phi is in fact an eigenfunction by computing 
c        integral operator applied to phi using legendre quadrature
c        
        nnodes=250
        call legeexps(1,nnodes,rhos,u,v,ws)
        do i=1,nnodes
        rhos(i)=rhos(i)/2.0d0+1/2.0d0
        enddo
c
c        tabulate integrand at nodes, rhos
c        
        eps=1.0d-50

        do 4500 k=1,nnodes
        do 4000 i=1,nnodes
        call dprol_t_bar_coef2exp(vect,rhos(i),nmat,nn,ip,w,val0i)
        dphis(i)=val0i
        x=c*rhos(k)*rhos(i)

        if (mod(ip,2) .eq. 0) then
        call rjbvals(x,eps,djs,nvals,nmax)
        endif

        if (mod(ip,2) .eq. 1) then
        call rjbvals_half(x,eps,djs,nvals,nmax)
        endif

        vals(i)=djs(nn+ip/2)*val0i*sqrt(x)
 4000   continue

c
c        compute integral by multiplying tabulation at gaussian 
c        nodes by weights 
c
        rint=0
        do j=1,nnodes
        rint=rint+ws(j)*vals(j)
        enddo
        rint=rint/2.0d0

c
c        compute approximate eigenvalue
c
        call dprol_t_bar_coef2exp(vect,rhos(k),nmat,nn,ip,w,val0)
        deigs(k)=rint/val0
        dmphis(k)=rint

 4500   continue

ccc        call prin2('deigs*', deigs, nnodes)

c
c        make sure that we do in fact have an eigenfunction by
c        comparing M_{N, c}[\phi] to \gamma_{N, n} \phi
c
        do i=1,nnodes
        dmphis_scaled(i)=dmphis(i)
        dds(i)=dmphis_scaled(i) - dphis(i) * dgam
        enddo
        call prin2('dds*',dds,nnodes)

        return
        end
c
c
c
c
c
        subroutine m_n_eigen_eval_rough(neig,nn,ip,c,w,dgam,vect)
        implicit real *8 (a-h,o-z)
        real*8 rhos(10000),vals(10000), djs(10000),deigs(10000)
        real*8 ws(10000),dds(10000),vect(1),w(1)
        real*8 dphis(10000),dmphis(10000),dmphis_scaled(10000)

c
c         evaluate an eigenvalue of the integral operator m in a 
c         primative way.
c         n.b. this subroutine only works for very small n,nn,c,ip
c         and is only used to verify forumlas, not compute
c

c
c         find eigenvectors of matrix corresponding to differential
c         operator l
c
        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)

c
c        find corresponding eigenvalue of integral operator M_{N, c}
c
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam,vect)
        call prin2('dgam*',dgam,1)

c
c         lay down gauss legendre nodes on (0,1)
c        
        nnodes=250
        call legeexps(1,nnodes,rhos,u,v,ws)
        do i=1,nnodes
        rhos(i)=rhos(i)/2.0d0+1/2.0d0
        enddo

c
c        tabulate integrand at nodes, rhos
c        
        eps=1.0d-50

        do 4500 k=1,nnodes
        do 4000 i=1,nnodes
        call dprol_t_bar_coef2exp(vect,rhos(i),nmat,nn,ip,w,val0i)
        dphis(i)=val0i
        x=c*rhos(k)*rhos(i)

        if (mod(ip,2) .eq. 0) then
        call rjbvals(x,eps,djs,nvals,nmax)
        endif

        if (mod(ip,2) .eq. 1) then
        call rjbvals_half(x,eps,djs,nvals,nmax)
        endif

        vals(i)=djs(nn+ip/2)*val0i*sqrt(x)
 4000   continue

c
c        compute integral by multiplying tabulation at gaussian 
c        nodes by weights 
c
        rint=0
        do 4100 j=1,nnodes
        rint=rint+ws(j)*vals(j)
 4100   continue
        rint=rint/2.0d0

c
c        compute approximate eigenvalue
c
        call dprol_t_bar_coef2exp(vect,rhos(k),nmat,nn,ip,w,val0)
        deigs(k)=rint/val0
        dmphis(k)=rint
        
ccc        call prin2('rint2*',rint2,1)
ccc        call prin2('function at nodes*',vals,nnodes)

 4500   continue

ccc        call prin2('deigs*', deigs, nnodes)

c
c        make sure that we do in fact have an eigenfunction by
c        comparing M_{N, c}[\phi] to \gamma_{N, n} \phi
c
        do i=1,nnodes
        dmphis_scaled(i)=dmphis(i)
        dds(i)=dmphis_scaled(i) - dphis(i) * dgam
        enddo
        call prin2('dds*',dds,nnodes)

c
c        and plot
c
        iw=18
        call quagraph2(iw,rhos,dphis,nnodes,3,
     1     rhos,dmphis_scaled,nnodes,3,'title*')

        return
        end
c
c
c
c
c
        subroutine test_lnc_93(w)
        implicit real*8 (a-h,o-z)
        real*8 w(1), vect(10 000)

c
c        check that \phi is an eigenfunction of the linear 
c        operator L_{N, c} of https://arxiv.org/pdf/1811.02733.pdf
c        (see (92)) of the original arxiv version. 
c       
        neig=8
        nn=2
        ip=0
        c=0
        
c
c        compute zernike expansion of the eigenfunction
c        phi as well as the corresponding eigenvalue of the 
c        differential operator L_{N, c}
c
        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)
ccc        call prin2('vect*',vect,nmat)

c
c         evaluate phi and its first two derivatives
c       
        r = 0.2
        call dprol_t_bar_coef2exp2(vect,r,nmat,nn,ip,w,val0,val1,
     1     val2)

        t1 = (1 - r**2) * val2 + val1 * (-2*r)
        t2 = (0.25 - (nn + ip/2.0d0)**2)/r**2 - c**2*r**2
        t2 = t2 * val0
        dd = t1 + t2 - rlam * val0
        call prin2('rlam*', rlam, 1)
        call prin2('dd*', dd, 1)

        return
        end
c
c
c
c
c
        subroutine exp_quad_test(w)
        implicit real*8 (a-h,o-z)
        real*8 w(*),droots(10 000),ws(10 000),rints(10 000)
        real*8 vect(10 000),fs(10 000),rs(10 000),djs(10 000)
        real*8 dthetas(10 000)
        complex*16 ima,fs_comp(100 000),rint_comp,rint_comp2,dd
        data ima/(0.0d0,1.0d0)/

c
c        integrate a complex exponential on the unit disk using
c        prolate quadrature in the radial direction and equispaced 
c        nodes in the angular direction. check accuracy by doubling
c        nodes in both angular and radial directions and integrating
c        again
c

        nn=0
        ip=1
        c=20

        nrad=20
        nang=50

        xx=0.9d0
        xy=0.2d0

c
c        generate radial nodes and weights
c
        call dprol_bigphi_exps(nrad,nn,ip,c,w,rs,ws)
ccc        call dprol_gausquad(nrad,ip,c,w,rs,ws)
        dsum=0
        do i=1,nrad
        dsum=dsum+ws(i)
        enddo

c
c        and angular ones
c

        done=1
        pi=atan(done)*4

        do i=1,nang
        dthetas(i)=2*pi*(i-1)/nang
        enddo

        eps=1.0d-20
        call rjbvals(c*xx,eps,djs,nvals,nmax)
        val=pi*2*djs(1)/(c*xx)

c
c        . . . and start quadrature
c

        rint_comp=0
        icount=1

        do j=1,nrad
        do i=1,nang
        tx=rs(j)*cos(dthetas(i))
        ty=rs(j)*sin(dthetas(i))

        scap=tx*xx+ty*xy

        fs_comp(icount)=exp(ima*c*scap)
        rint_comp=rint_comp+ws(j)*2*pi/nang*fs_comp(icount)
        icount=icount+1

        enddo
        enddo

ccc        call prin2('fs_comp*',fs_comp,nrad*nang*2)
        call prin2('rint_comp*',rint_comp,2)

c
c        now do the same thing, except with double the number of nodes 
c        in both the radial and angular directions
c

        nrad=nrad*2
        nang=nang*2

c
c        generate radial nodes and weights
c
        
        call dprol_bigphi_exps(nrad,nn,ip,c,w,rs,ws)
ccc        call dprol_gausquad(nrad,ip,c,w,rs,ws)

c
c        and angular ones
c

        done=1
        pi=atan(done)*4

        do 1200 i=1,nang
        dthetas(i)=2*pi*(i-1)/nang
 1200    continue

c
c        . . . and start quadrature
c

        rint_comp2=0
        icount=1

        do j=1,nrad
        do i=1,nang
        tx=rs(j)*cos(dthetas(i))
        ty=rs(j)*sin(dthetas(i))

        scap=tx*xx+ty*xy
        fs_comp(icount)=exp(ima*c*scap)

        rint_comp2=rint_comp2+ws(j)*2*pi/nang*fs_comp(icount)
        icount=icount+1
        enddo
        enddo

ccc        call prin2('fs_comp*',fs_comp,nrad*nang*2)
        call prin2('rint_comp2*',rint_comp2,2)

        dd=(rint_comp2-rint_comp)/rint_comp2
        call prin2('dd*',dd,2)

        return
        end
c
c 
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning
c        of the prolate routines
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains a set of subroutines for the handling
c        of prolate spheroidal wave functions. The following is a brief 
c        description of the user-callable subroutines contained in 
c        this file.
c 
c   dprol_gausquad - generates the n nodes and weights that integrate the
c       first 2n gpsfs of an inputted bandlimit
c
c   dprol_gammas_eval - evaluates the first n eigenvalues of differential 
c       operator L, the first n eigenvalues of integral operator m,
c       and the coefficients of the t-bar expansions of the first n 
c       eigenfunctions. Note that operators m and L commute and have 
c       the same eigenfunctions.
c
c   dprol_gamma_n_eval - evaluates the n^th eigenvalue of differential 
c       operator L, the n^th eigenvalue of integral operator m, and 
c       the expansion of the n^th eigenfunction in t bar functions.
c
c   dprol_t_bar_coef2exp - given a user-provided array of coefficients, 
c       this subroutine evaluates an expansion in t bar functions at a 
c       user-provided point on the interval [0,1].
c
c   dprol_log_phi_star_eval - evaluates log(phi^*), that is, the 
c       coefficient of the r^{nn+2n} term of \phi_{nn,n}.
c
c   dprol_ln_all_eigen_eval - evaluates the first n eigenvalues and 
c       the first n eigenvectors of the symmetric tridiagonal matrix 
c       corresponding to the differential operator L using sturm 
c       bisection
c
c   dprol_ln_n_eigen_eval - evaluates the n^th eigenvalue and the n^th 
c       eigenvector of the symmetric tridiagonal matrix corresponding 
c       to the differential operator L using sturm bisection
c
c   dprol_mat_size_evalq - this subroutine determines a number of terms 
c       in the t-bar expansion of the n^th eigenfunction 
c       of integral operator m (and differential operator L) sufficient
c       to capture all terms with coefficients of magnitude greater 
c       than 1.0d-35. This subroutine is used to determine the size of 
c       the (symmetric and tridiagonal) matrix corresponding to the 
c       differential operator L whose eigenvectors form the coefficients
c       of the t-bar expansions of eigenfunctions of L and m.
c
c   dprol_ratio_eval - evaluates the ratio \gamma_{nn,i}/\gamma_{nn,j}, 
c       given the user-provided coefficients of the expansions in 
c       t bar functions of the two corresponding eigenfunctions 
c       phi_{nn,i} and phi_{nn,j}.
c
c   dprol_t_bar - evaluates a t bar function. that is, a normalized 
c       zernike polynomial multiplied by r^{(p+1)/2}
c
c   dprol_zern_pol - evaluates the (non-normalized) zernike 
c       polynomial R_{nn,n}
c
c   dprol_zern_normd_pol - evaluates a normalized zernike polynomial. 
c       that is, \sqrt{2(2n+nn+p/2+1)}*R_{nn,n}, which is l2 
c       normalized with weight function r^{p+1}
c
c   dprol_jac_pol - evaluates a jacobi polynomial with parameter \Beta=0
c       and user provided \alpha
c   
c   dprol_t_bar_derivs - evaluates the first n t bar functions and their 
c       derivatives
c
c   dprol_zern_pol_derivs - evaluates the first n non-normalized zernike 
c       polynomials along with their derivatives
c
c   dprol_jac_pol_derivs - evaluates the first n jacobi polynomials (with 
c       parameter \Beta=0) along with their derivatives
c
c   dprol_zern_pol_deriv2 - evaluates a non-normalized zernike 
c       polynomial along with its first and second derivatives
c
c   dprol_jac_pol_deriv2 - evaluates a jacobi polynomial 
c       (with parameter \Beta=0) along with its first and second 
c       derivatives
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine dprol_gausquad(n,ip,c,w,rs,ws)
        implicit real*8 (a-h,o-z)
        real*8 w(1),ws(1),rs(1),vect(10 000),a2(100 000),
     1     rints(10 000),rints_true(10 000),dds(10 000),adws(10 000),
     1     adrs(100 000),rints2(10 000),rhs(10 000),dders(10 000),
     1     rws(10 000),rnorms(10 000),a(100 000),x2(10 000),u(10 000),
     1     v(10 000),s(10 000),rs2(10 000),ws2(10 000)

c
c        generate a gaussian quadrature for the gpsfs with bandlimit c
c
c          \Phi_{0,0}, \Phi_{0,1},..., \Phi_{0,2n-1}
c
c        by first generating a chebyshev quadrature for gpsfs with 
c        bandlimit c/2 and then doing newton on the first 2n gpsfs 
c        with bandlimit c
c

c
c        generate n nodes and weights corresponding to c/2
c

        nn=0
        c2=c/2.0d0
        call dprol_bigphi_exps(n,nn,ip,c2,w,rs,ws)

        call dprol_copy(rs,rs2,n)
        call dprol_copy(ws,ws2,n)

c
c        ...and start newton
c

        do 1400 ijk=1,10

        call dprol_gaus_newtstep(n,nn,ip,c2,w,rs,ws,scap_new)

c
c        if you start to see some convergence do 3 more steps
c

        if (scap_new .lt. 1.0d-16*n) goto 2200

        if (scap_new .lt. 1.0d-5*n) then
          do 2000 ijkl=1,3
          call dprol_gaus_newtstep(n,nn,ip,c2,w,rs,ws,scap_new)
 2000     continue
          goto 2200
        endif

 1400   continue
 2200   continue


        return
        end
c
c
c
c
c
        subroutine dprol_gaus_newtstep(n,nn,ip,c,w,rs,ws,scap_new)
        implicit real*8 (a-h,o-z)
        real*8 w(1),ws(1),rs(1),vect(10 000),a2(100 000),
     1     rints(10 000),rints_true(10 000),dds(10 000),adws(10 000),
     1     adrs(100 000),rints2(10 000),rhs(10 000),dders(10 000),
     1     rws(10 000),rnorms(10 000),a(100 000),x2(10 000),u(10 000),
     1     v(10 000),s(10 000),rs2(10 000),ws2(10 000)

c
c        does one newton step in generating gpsf guassian quadratures.
c        takes in n nodes (rs) and n weights (ws) and modifies them
c        by taking one newton step
c        note: nn should be 0

c
c        first evaluate the discrepencies 
c

        call dprol_discrep_eval(rs,ws,n,nn,ip,c,w,rints,
     1     rints_true,adws,adrs,dds)

        call dprol_scap(dds,dds,2*n,scap)

c
c        combine nodes and weights matrices into matrix a and copy
c

        call dprol_copy(adrs,a,2*n*n)
        call dprol_copy(adws,a(2*n*n+1),2*n*n)

        call dprol_copy(a,a2,2*n*2*n)

c
c        . . . now do a step of newton
c

        do i=1,2*n
        x2(i)=-dds(i)
        enddo
        call qrsolv(a,2*n,x2,rcond)

c
c        now tweak rs and ws accordingly
c        
        h=1.0d0

 350    continue

        do i=1,n
        rs2(i)=rs(i)+h*x2(i)
        ws2(i)=ws(i)+h*x2(n+i)
        enddo

c
c        check new discrepencies
c

        call dprol_discrep_eval(rs2,ws2,n,nn,ip,c,w,rints,
     1     rints_true,adws,adrs,dds)
        call dprol_scap(dds,dds,2*n,scap_new)

c
c         if vector of discrepencies has smaller l2 norm then continue
c

        if (scap_new .lt. scap) then
          call dprol_copy(rs2,rs,2*n)
          call dprol_copy(ws2,ws,2*n)
          goto 1300 
        endif

c
c        . . . otherwise, divide the steplength by 2 and try again
c

        h=h/2.0d0
        call dprol_copy(rs,rs2,2*n)
        call dprol_copy(ws,ws2,2*n)
        goto 350

 1300   continue

        return
        end
c
c
c
c
c
        subroutine dprol_discrep_eval(rs,ws,n,nn,ip,c,w,rints,
     1     rints_true,adws,adrs,dds)
        implicit real*8 (a-h,o-z)
        real*8 w(1),ws(1),rints(1),rints_true(1),dds(1)
        real*8 vect(10 000),fs(10 000),rs(10 000),adws(2*n,n)
        real*8 rs2(10 000),ws2(10 000),fs2(10 000),adrs(2*n,n)
        real*8 vects(100 000),rlams(10 000),dgams(10 000)
        real*8 fs1(10 000)

c
c        with the n nodes (rs) and n weights (ws), integrate the gpsfs
c
c          \Phi_{0,0}(r)*r**(ip+1),\Phi_{0,1}(r)*r**(ip+1),
c                                  ...,\Phi_{0,2*n-1}(r)*r**(ip+1)
c
c        with bandlimit 2c and compare with true values and save errors
c        in dds. also generate matrices of partial derivatives, adrs 
c        and adws. this subroutine is used in the construction of 
c        gaussian quadratures for gpsfs
c

        c2=2*c

c
c        evaluate the first 2*n prolates with quadrature and double
c        the bandlimit
c

        neigs=2*n
        call dprol_gammas_eval(neigs,nn,ip,2*c,w,nmat,rlams,
     1     dgams,vects)


        do 1000 j=1,neigs

        neig=j
        ijk=(j-1)*nmat+1

c
c        . . . tabulate bigphi at nodes
c

        do 200 i=1,n
        r=rs(i)
        call dprol_normd_zern_exp_deriv(vects(ijk),r,nn,ip,nmat,
     1     val0,val1)
        fs(i)=val0
        fs1(i)=val1
        adrs(j,i)=fs1(i)*ws(i)
        adws(j,i)=fs(i)
 200    continue

ccc        call prin2('fs*',fs,n)

c
c        . . . compute integral
c        

        rint=0
        do i=1,n
        rint=rint+fs(i)*ws(i)
        enddo

        rints(j)=rint

 1000   continue

c
c        . . . compute true integrals 
c

        n2=neigs+2
        call dprol_bigphi_exps(n2,nn,ip,c2,w,rs2,ws2)

        do 1200 j=1,neigs

        ijk=(j-1)*nmat+1
c
c        . . . tabulate bigphi at nodes
c

        do i=1,n2
        r=rs2(i)
        call dprol_normd_zern_exp_deriv(vects(ijk),r,nn,ip,nmat,
     1     val0,val1)
        fs2(i)=val0
        enddo

c
c        . . . compute true integral
c        

        rint2=0
        do i=1,n2
        rint2=rint2+fs2(i)*ws2(i)
        enddo

        rints_true(j)=rint2

        dds(j)=(rints(j)-rints_true(j))

 1200   continue
 
        return
        end
c
c
c
c
c
        subroutine dprol_bigphi_exps(n,nn,ip,c,w,droots,ws)
        implicit real*8 (a-h,o-z)
        real*8 w(1),droots(1),rhs(10 000),a(100 000),ws(1)
        real*8 rnorms(10 000)

c
c        this subroutine finds the n roots of \Phi_{0,n} (droots).
c        it also finds the n weights (ws) that integrate
c        from 0 to 1 the functions
c
c           \Phi_{0,i}(r)r**(ip+1)
c
c        for i=0,...,n-1

        neig=n+1
c
c        find roots of phi
c
        call dprol_phi_roots_eval(neig,nn,n,ip,c,w,droots)

c
c        create interpolation matrix
c

        call dprol_bigphi_tabmat(droots,n,nn,ip,c,w,a)

c
c        solve linear system where rhs contains the true integrals
c
        nphis=n
        call dprol_bigphi_quads(nphis,ip,c,w,rhs)
        call qrsolv(a,n,rhs,rcond)
        call dprol_copy(rhs,ws,n)

        return
        end
c
c
c
c
c
        subroutine dprol_bigphi_quads(nphis,ip,c,w,rints)
        implicit real*8 (a-h,o-z)
        real*8 w(1),droots(10 000),rs(10 000),ws(10 000)
        real*8 fs(10 000),rints(*),vects(100 000),rlams(10 000)
        real*8 dgams(10 000),rints2(10 000),dds(10 000)
        real*8 rints3(10 000), u(100 000), v(100 000)

c
c        use nodes and weights that integrate polynomials with respect
c        to weight function r**(ip+1) to integrate the gpsfs
c
c           \Phi_{0,0},\Phi_{0,1},...,\Phi_{0,nphis-1)
c
c        and store values in rints
c

c
c        first, find out how many nodes are needed 
c       
        neigs=nphis
        nn=0

        call dprol_gammas_eval(neigs,nn,ip,c,w,nmat,rlams,
     1     dgams,vects)

        nnodes=nmat+1

c
c        lay down xk-pol nodes
c

        k=ip+1
        itype=1
        call xk_exps(itype,nnodes,k,w,rs,u,v,ws)

c
c        . . . integrate each phi_{0,n} from n-1=1,...,nphis
c

        nn=0
        do 1000 ijk=1,nphis
        neig=ijk

        do 600 i=1,nnodes
        r=rs(i)

        call dprol_normd_zern_exp_eval(vects((ijk-1)*nmat+1),r,nn,ip,
     1     nmat,val)

        fs(i)=val
 600    continue

        rint1=0

        do 800 i=1,nnodes
        rint1=rint1+fs(i)*ws(i)
 800    continue
        rints(ijk)=rint1

 1000   continue

        return
        end
c
c
c
c
c
        subroutine dprol_bigroot_cheb_eval(neig,nn,ip,c,w,x1)
        implicit real *8 (a-h,o-z)
        real*8 w(1),vect(10 000),droots(10 000)

c
c        check for the largest root of \Phi_{0,neig-1} by placing
c        chebyshev nodes on the interval (0,x_0) where x_0 is the
c        root of the coefficient of the 0^th order derivative in 
c        the ode satisfied by \Phi_{0,neig-1}, after finding a 
c        sign change, use newton to tighten
c        note: in this subroutine, nn should equal 0
c

        done=1
        pi=atan(done)*4

c
c        use a primitive algorithm to find the roots of phi
c        

        n=neig-1
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam,vect)

c
c        find right end 
c

        call dprol_coefzero(neig,nn,ip,c,rlam,w,bb)

c
c        lay down chebyshev nodes on (0,bb) and check for a sign
c        change immediately
c

        r=bb
        call dprol_t_bar_coef2exp2(vect,r,nmat,nn,ip,w,val0,
     1     val1,val2)

        aa=0.0

        nnodes=neig*10

        do 100 i=0,nnodes-1

        dtheta=pi*(2*(i+1)-1)/(2*nnodes)
        r1=dcos(dtheta)
        r1=r1*(bb-aa)/2.0d0+(bb-aa)/2.0d0

        call dprol_t_bar_coef2exp2(vect,r1,nmat,nn,ip,w,val0p,
     1     val1,val2)

        dprod=val0*val0p

        if (dprod .lt. 0) then 
          x0=r
          f0=val0
          goto 400
        endif

        val0=val0p
        r=r1

 100    continue
 400    continue

c
c        . . . do newton to tighten
c

        do j=1,20
        call dprol_t_bar_coef2exp2(vect,x0,nmat,nn,ip,w,val0,
     1     val1,val2)

        x1=x0-val0/val1
        x0=x1
        enddo

        val=x1

 150    continue

        return
        end
c
c
c
c
c
        subroutine dprol_coefzero(neig,nn,ip,c,rlam,w,x0)
        implicit real *8 (a-h,o-z)
        real*8 rs(100 000),phis(100 000),w(1),vect(100 000)

c
c        use bisection to find root (in x) of 
c
c            (1/4-p^2/4)/x**2-c**2*x**2-rlam  \in (0,1)
c
c        the above function is the coefficient of the 0^th order term
c        in the ode satisfied by phi (rlam is the eigenvalue of the 
c        differential operator corresponding to phi)
c

c
c        start bisection 
c
        xm=0
        xp=1

        
        fm=(1/4.0d0-ip**2/4.0d0)/xm**2-c**2*xm**2-rlam
        fp=(1/4.0d0-ip**2/4.0d0)/xp**2-c**2*xp**2-rlam

        if (fp .gt. 0) then 
          x0=1
          goto 2000
        endif

        do 150 i=1,100

        x0=(xp+xm)/2.d0
        f0=(1/4.0d0-ip**2/4.0d0)/x0**2-c**2*x0**2-rlam

        if (f0 .gt. 0) then 
          xm=x0
          fm=f0
          goto 100
        endif

        xp=x0
        fp=f0

 100    continue
 150    continue

 2000   continue

        return
        end
c
c
c
c
c
        subroutine dprol_phi_big_root(neig,nn,ip,c,w,val)
        implicit real *8 (a-h,o-z)
        real*8 rs(100 000),fs(100 000),fs2(10 000),droots(10 000)
        real*8 w(1),vect(10 000)

c
c        find the largest root of phi_{nn,n} 
c        note: nn should be 0. this has not been tested for nn>0
c

        n=neig-1
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam,vect)

c
c        if dgam is big, use chebyshev node strategy 
c        

        dgam_big=1/sqrt(c)
        if (abs(dgam) .gt. dgam_big/10) then
          call dprol_bigroot_cheb_eval(neig,nn,ip,c,w,x0)
          goto 1100
        endif

c
c        otheriwse, evaluate phi and derivs at x0, the root of 
c        the 0^th order term in the ode satisfied by phi
c

        call dprol_coefzero(neig,nn,ip,c,rlam,w,x0)

c
c        do 5 steps of meuller
c

        do 600 ijk=1,5
 550    continue

        call dprol_t_bar_coef2exp2(vect,x0,nmat,nn,ip,w,val0,
     1     val1,val2)

c
c        . . . fit to quadratic
c

        a0=val2/2.0d0
        b0=val1-val2
        c0=val0-val1+val2/2.0d0

        aa=0.1
        bb=1.0
        nnodes=100
        do 200 i=1,nnodes
        rs(i)=aa+(bb-aa)*(i-1)/(nnodes-1)

        call dprol_t_bar_coef2exp2(vect,rs(i),nmat,nn,ip,w,val0,
     1     val1,val2)
        fs(i)=val0
 200    continue

        droot1=(-b0+sqrt(b0**2-4*a0*c0))/(2*a0)
        droot2=(-b0-sqrt(b0**2-4*a0*c0))/(2*a0)

c
c        . . . select the right root as the next guess
c

        if (droot1 .gt. x0) then
          x0=droot2
          goto 400
        endif

        if (droot2 .gt. x0) then 
          x0=droot1
          goto 400
        endif

        x0=droot1
        if (droot2 .gt. droot1) then 
          x0=droot2
          goto 400
        endif

 400    continue

 600    continue

c
c        . . . now some newton
c

        do 800 j=1,50

        call dprol_t_bar_coef2exp2(vect,x0,nmat,nn,ip,w,val0,
     1     val1,val2)
        x0=x0-val0/val1
c
c        when we start to see some convergence, take 5 more steps
c
        if (val0/val1 .lt. 1.0d-8) then
          do 700 ijk=1,5
          call dprol_t_bar_coef2exp2(vect,x0,nmat,nn,ip,w,val0,
     1     val1,val2)
          x0=x0-val0/val1
 700      continue
          goto 1000
        endif

 800    continue
 1000   continue


 1100   continue

        val=x0

        return
        end
c
c
c
c
c
        subroutine dprol_bigphi_tabmat(rs,n,nn,ip,c,w,a)
        implicit real*8 (a-h,o-z)
        real*8 w(1),a(1),vects(100 000),rs(1)
        real*8 dgams(10 000),rlams(10 000)

c
c        create the n x n matrix, a, of gpsfs tabulated at n points, rs.
c        that is, create a where
c
c         a(i,j)=\Phi_{nn,i-1}(rs(j))
c

c
c        . . . first evaluate the zernike expansions of phi
c

        neigs=n
        call dprol_gammas_eval(neigs,nn,ip,c,w,nmat,rlams,dgams,vects)

        do 1000 j=1,n
        r=rs(j)

        do 800 i=1,n
        ijk=(i-1)*nmat+1

        call dprol_normd_zern_exp_eval(vects(ijk),r,nn,ip,nmat,val0)

        a((j-1)*n+i)=val0
 800    continue
 1000   continue

        return
        end
c
c
c
c
c
        subroutine dprol_normd_zern_exp_eval(cs,r,nn,ip,nzerns,val)
        implicit real*8 (a-h,o-z)
        real*8 cs(1),vals0(10 000),vals1(10 000)

c
c        evaluate an nzerns-term expansion in normed zernike 
c        polynomials for fixed nn
c

        call dprol_zern_normd_derivs(r,nzerns,nn,ip,vals0,vals1)

        val=0
        do 400 i=1,nzerns
        val=val+cs(i)*vals0(i)
 400    continue

        return
        end
c
c
c
c
c
        subroutine dprol_normd_zern_exp_eval2(cs,r,nn,ip,nzerns,val0,
     1     val1,val2)
        implicit real*8 (a-h,o-z)
        real*8 cs(1),vals0(10 000),vals1(10 000),vals2(10 000)

c
c        evaluate an nzerns-term expansion in normed zernike 
c        polynomials for fixed nn
c

        call dprol_zern_normd_derivs2(r,nzerns,nn,ip,vals0,vals1,
     1     vals2)

        val0=0
        val1=0
        val2=0
        do i=1,nzerns
        val0=val0+cs(i)*vals0(i)
        val1=val1+cs(i)*vals1(i)
        val2=val2+cs(i)*vals2(i)
        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_normd_zern_exp_deriv(cs,r,nn,ip,nzerns,val0,
     1     val1)
        implicit real*8 (a-h,o-z)
        real*8 cs(1),vals0(10 000),vals1(10 000)

c
c        evaluate an expansion in normed zernike polynomials and the
c        derivatives of normed zernike polynomials
c

        call dprol_zern_normd_derivs(r,nzerns,nn,ip,vals0,vals1)

        val0=0
        val1=0
        do 400 i=1,nzerns
        val0=val0+cs(i)*vals0(i)
        val1=val1+cs(i)*vals1(i)
 400    continue

        return
        end
c
c
c
c
c
        subroutine dprol_matvec(a,x,n,y)
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
        subroutine dprol_scap(x,y,n,scap)
        implicit real*8 (a-h,o-z)
        real*8 x(1),y(1)

        scap=0
        do 100 i=1,n
        scap=scap+x(i)*y(i)
 100    continue

        return
        end
c
c
c
c
c
        subroutine dprol_phi_roots_eval(neig,nn,n,ip,c,w,vals)
        implicit real*8 (a-h,o-z)
        real*8 vals(1),w(1),vect(10 000)

c
c        find the n-1 roots of phi_{0,n-1} by first using Prufer 
c        Transform to find initial guesses then using newton to tighten
c        note: this subroutine has only been tested with nn=0
c

c
c        find zernike expansion of phi
c
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam,vect)
        n=neig-1

c
c        . . . find the biggest root
c

        call dprol_phi_big_root(neig,nn,ip,c,w,val)
        vals(1)=val

        done=1
        pi=atan(done)*4

c
c        march to get other roots
c

        nsteps=100
        x0=vals(1)
        do 200 i=1,n-1
        
        theta0=(i-1)*pi+pi/2
        thetat=theta0+pi
        call dprol_theta_ode_panel(neig,nn,ip,c,nsteps,x0,theta0,
     1     thetat,w,xt)

        vals(i+1)=xt
        x0=xt

 200    continue

c
c        tighten each root via newton
c

        do 800 i=1,n
        do 600 j=1,50

        call dprol_t_bar_coef2exp2(vect,vals(i),nmat,nn,ip,w,val0,
     1     val1,val2)
        vals(i)=vals(i)-val0/val1

c
c        when we start to see some convergence, take 5 more steps
c
        if (val0/val1 .lt. 1.0d-8) then
          do 400 ijk=1,5
          call dprol_t_bar_coef2exp2(vect,vals(i),nmat,nn,ip,w,val0,
     1     val1,val2)
          vals(i)=vals(i)-val0/val1
 400      continue
          goto 800
        endif

 600    continue
 800    continue

c
c        reorder array, putting roots in increasing order
c

        do 1000 i=1,n/2
        val=vals(i)
        vals(i)=vals(n+1-i)
        vals(n+1-i)=val
 1000   continue

        return
        end
c
c
c
c
c
        subroutine dprol_theta_ode_panel(neig,nn,ip,c,nsteps,x0,
     1     theta0,thetat,w,xt)
        implicit real*8 (a-h,o-z)
        real*8 rs(100000),dthetas(100000),w(1),droots(1000)
        real*8 fs(10 000),xs(10 000),vect(10 000),dthetas2(10 000)

c
c        on the interval (theta0,thetat) use RK2 to solve the ODE
c        dx/dtheta=(-sqrt(b)-(b1/(4*b)+a/2.0d0)*sin(2*theta))**(-1)
c        which is the ODE of the Prufer transform. 
c

        n=neig-1

c
c        . . . compute the eigenvalue, rlam
c

        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)

c
c        do one ode panel, from theta0 to thetat w/ initial condition
c        x(theta0) is the biggest root
c

        done=1
        pi=atan(done)*4

        h=(thetat-theta0)/(nsteps-1)
        dthetas(1)=theta0

        do 100 i=1,nsteps
        dthetas(i+1)=dthetas(i)+h
 100    continue

c
c        now march
c
        xs(1)=x0

        do 200 i=1,nsteps-1

        r=xs(i)

        dthet=dthetas(i)

        beta0=((1/4.0d0-(nn+ip/2.0d0)**2))/(r**2*(1-r**2))
     1     -(c**2*r**2+rlam)/(1-r**2)

        beta1=-(1/4.0d0-(nn+ip/2.0d0)**2)*2*(1-2*r**2)
     1     /((1-r**2)**2*r**3)
     1     +((1-r**2)*(-c**2*2*r)+(-c**2*r**2-rlam)*2*r)/((1-r**2)**2)

        alpha0=(-2*r)/(1-r**2)

        dtdx=-sqrt(beta0)-(beta1/(4*beta0)+alpha0/2.0d0)*
     1     sin(2*dthet)

        der=dtdx**(-1)

        xs(i+1)=xs(i)+h*der

c
c        tighten via rk2
c

        rp=xs(i+1)
        dthetp=dthetas(i+1)

        beta0=((1/4.0d0-(nn+ip/2.0d0)**2))/(rp**2*(1-rp**2))
     1     -(c**2*rp**2+rlam)/(1-rp**2)

        beta1=-(1/4.0d0-(nn+ip/2.0d0)**2)*2*(1-2*rp**2)
     1     /((1-rp**2)**2*rp**3)
     1     +((1-rp**2)*(-c**2*2*rp)+(-c**2*rp**2-rlam)*2*rp)
     1     /((1-rp**2)**2)

        alpha0=(-2*rp)/(1-rp**2)

        dtdx=-sqrt(beta0)-(beta1/(4*beta0)+alpha0/2.0d0)*
     1     sin(2*dthetp)

        derp=dtdx**(-1)

        xs(i+1)=xs(i)+h/2.0d0*(der+derp)

 200    continue

c
c        set final point as xt for output
c
        xt=xs(nsteps)

        return
        end
c
c
c
c
c
        subroutine dprol_gamma_n_eval(neig,nn,ip,c,w,nmat,rlam,dgam,
     1     vect)
        implicit real *8 (a-h,o-z)
        real*8 w(1),vect(1)

c
c         This subroutine evaluates the neig^th eigenvalue of 
c         integral operator m,
c
c           \gamma_{nn,neig-1},
c
c         the neig^th eigenvalue of differential operator l,
c
c           \lambda_{nn,neig-1},
c
c         and the corresponding eigenfunction expansion in t bar 
c         functions
c
c           \phi_{nn,neig-1}.
c
c         This subroutine evaluates \gamma_{nn,neig-1} by first 
c         evaluating log(phi*_{nn,neig-1}), which is log of the coefficient 
c         of the r^{nn+2(neig-1)} term of phi_{nn,neig-1} and then 
c         evaluating log(k*a_0) where k depends on nn,c,ip, and a_0
c         is the first coefficient of the eigenfunction expansion (vect(1)).
c
c         The integral operator m is defined by the formula:
c
c                            1
c                           |   
c         m_{N,c}[Phi](r) = |   J_{N+ip/2}(crp)*sqrt(crp)*Phi(p) dp
c                           |   
c                          0
c         
c       input parameters:
c   
c  neig - the number of eigenfunction and eigenvalues to be computed
c  nn - the azimuthal number
c  ip - we are in dimension ip+2
c  c - the bandlimit
c
c  w - work array that must be allocated at least 14*nmat+220 real*8 
c         locations long
c
c       output parameters:
c
c  nmat-number of terms in the eigenfunction expansion of vect
c  rlam-neig^th eigenvalue of the differential operator
c  dgam-neig^th eigenvalue of the symmetric integral 
c       operator, m
c  vect-coefficients of the t bar-expansion of the neig^th
c       eigenfunction
c

c
c         compute the log(phi*_{nn,neig-1})
c
        call dprol_log_phi_star_eval(neig,nn,ip,c,w,nsign,dlog)

c
c        . . . now compute log of the integral
c
        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)

        dnn1=nn+ip/2.0d0+1
        call gammanew_eval_log_extend(dnn1,dlgam)

        dlogrint=(nn+(ip+1)/2.0d0)*log(c)+log(abs(vect(1)))
     1     -(nn+ip/2.0d0)*log(2.0d0)-log(2*nn+ip+2.0d0)/2.0d0-dlgam
        nsign=sign(1.0d0,nsign*vect(1))

c
c        . . . exponentiate to get solution
c
        dgam=exp(dlogrint-dlog)*nsign

        return
        end
c
c
c
c
c
        subroutine dprol_gammas_eval(neigs,nn,ip,c,w,nmat,rlams,
     1     dgams,vects)
        implicit real *8 (a-h,o-z)
        real*8 vects(1),rlams(1),dgams(1)

c
c         This subroutine evaluates the first neigs eigenvalues of 
c         integral operator m,
c
c           \gamma_{nn,0},\gamma_{nn,1},...,\gamma_{nn,neigs-1},
c
c         the first neigs eigenvalues of differential operator L,
c
c           \lambda_{nn,0},\lambda_{nn,1},...,\lambda_{nn,neigs-1},
c
c         and the expansions in t bar functions of the corresponding 
c         eigenfunctions,
c
c           \phi_{nn,0},\phi_{nn,1},...,\phi_{nn,neigs-1}.
c
c         In order to evaluate \phi_{nn,n} at a point, use the nmat
c         coefficients vect(n*nmat+1),...,vect(n*nmat+(nmat-1) with
c         the subroutine dprol_t_bar_coef2exp.
c
c         The integral operator m is defined by the formula:
c
c                            1
c                           |   
c         m_{N,c}[Phi](r) = |   J_{N+ip/2}(crp)*sqrt(crp)*Phi(p) dp
c                           |   
c                          0
c         
c
c       input parameters:
c   
c  neigs - the number of eigenfunctions and eigenvalues to be computed
c  nn - the azimuthal number
c  ip - we are in dimension ip+2
c  c - the bandlimit
c
c  w - work array that must be allocated at least 14*nmat+220 real*8 
c         locations long
c
c
c       output parameters:
c
c  nmat - number of terms in the eigenfunction expansions stored in vects.
c       that is, vects contains neigs expansions of length nmat one after 
c       the other.
c  rlams - eigenvalues of the differential operator
c  dgams - eigenvalues of the symmetric integral operator, m
c  vects - coefficients of the t bar-expansions of the eigenfunctions.
c       that is, vects contains, one after the other, expansions
c       of length nmat of eigenfunctions of integral operator m (and 
c       (also differential operator L)
c

c
c        compute the first eigenvalue, \gamma_{nn,0}
c

        neig=1
        call dprol_gamma_n_eval(neig,nn,ip,c,w,nmat2,rlam,dgam0,vects)
        dgams(1)=dgam0

c       
c        compute the eigenvectors of the tridiagonal matrix for l
c
        call dprol_ln_all_eigen_eval(neigs,nn,ip,c,w,nmat,rlams,vects)
        indi = (neigs - 1) * nmat + 1

c
c        compute ratios and multiply
c
        do i=1,neigs
        i1=nmat*(i-1)+1
        i2=nmat*i+1
        call dprol_ratio_eval(nn,ip,nmat,vects(i1),vects(i2),r21,w)
        dgams(i+1)=dgams(i)*r21
        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar_coef2exp(cs,r,n,nn,ip,w,val0)
        implicit real*8 (a-h,o-z)
        real*8 cs(1),w(1)

c
c         given the n+1 coefficients (cs) of an n-order expansion in 
c         t bar functions, evaluate that expansion at r
c

c 
c        allocate memory
c
        it0s=1
        lt0s=n+10

        it1s=it0s+lt0s
        lt1s=n+10

        call dprol_t_bar_coef2exp0(cs,r,n,nn,ip,w(it0s),w(it1s),val0)

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar_coef2exp0(cs,r,n,nn,ip,t0s,t1s,val0)
        implicit real*8 (a-h,o-z)
        real*8 cs(1),t0s(1),t1s(1)

c
c         given the n+1 coefficients (cs) of an n-order expansion in 
c         t bar functions, evaluate that expansion at r
c
        val0=0

        call dprol_t_bar_derivs(r,n,nn,ip,t0s,t1s)
        do i=1,n+1
        val0=val0+cs(i)*t0s(i)
        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar_coef2exp2(cs,r,n,nn,ip,w,val0,val1,val2)
        implicit real*8 (a-h,o-z)
        real*8 cs(1),w(1)

c
c         given the n+1 coefficients (cs) of an n-order expansion in 
c         t bar functions, evaluate that expansion at r in addition 
c         to its first and second derivatives
c

c 
c        allocate memory
c

        it0s=1
        lt0s=n+10
c
        it1s=it0s+lt0s
        lt1s=n+10
c
        it2s=it1s+lt1s
        lt2s=n+10

        call dprol_t_bar_coef2exp20(cs,r,n,nn,ip,w(it0s),w(it1s),
     1     w(it2s),val0,val1,val2)

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar_coef2exp20(cs,r,n,nn,ip,t0s,t1s,t2s,
     1     val0,val1,val2)
        implicit real*8 (a-h,o-z)
        real*8 cs(1),t0s(1),t1s(1),t2s(1)

c
c         given the n+1 coefficients (cs) of an n-order expansion in 
c         t bar functions, evaluate that expansion at r, in addition
c         to its first and second derivatives
c
        val0=0
        val1=0
        val2=0

        call dprol_t_bar_derivs2(r,n,nn,ip,t0s,t1s,t2s)

        do i=1,n
        val0=val0+cs(i)*t0s(i)
        val1=val1+cs(i)*t1s(i)
        val2=val2+cs(i)*t2s(i)
        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_log_phi_star_eval(neig,nn,ip,c,w,nsign,dlog)
        implicit real *8 (a-h,o-z)
        real*8 w(1)

c
c         evaluate log(phi*_{nn,neig-1}) using logs of the summands 
c         so that solution is:
c            nsign*exp(dlog)
c

c
c        w should be 17*nmat+250
c 
c        allocate memory
c

        call dprol_mat_size_evalq(neig,c,nmat)
        n=nmat+1

        iw=1
        lw=12*n+200
c
        ivect=iw+lw
        lvect=n+10
c
        idlogs=ivect+lvect
        ldlogs=n+10
c
        isigns=idlogs+ldlogs
        lsigns=n+10
c
        ialphas=isigns+lsigns
        lalphas=n+10
c
        inis=ialphas+lalphas
        lnis=n+10


        call dprol_log_phi_star_eval0(neig,nn,ip,c,w(iw),
     1     w(ivect),w(idlogs),w(isigns),w(ialphas),w(inis),
     1     nsign,dlog)


        return
        end
c
c
c
c
c
        subroutine dprol_log_phi_star_eval0(neig,nn,ip,c,w,
     1     vect,dlogs,dsigns,alphas,nis,
     1     nsign,dlog)
        implicit real *8 (a-h,o-z)
        real*8 w(1),vect(1),dlogs(1),dsigns(1),alphas(1)
        integer nis(1)

c
c         evaluate log(phi*_{nn,neig-1}) using logs of the summands 
c         so that solution is:
c           nsign*exp(dlog)
c

c
c         compute the neig^th eigenfunction
c
        call dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,rlam,vect)

c
c        . . . and determine machine precision and compute pi
c
        call dprol_mach_zero(dmach0)
        done=1
        pi=atan(done)*4

c
c        . . . start loop for computing logs of summands
c
        dcoef=0
        dlogmax=-100000
        do 100 i=0,nmat

c
c        . . . compute binomial using logs of gamma function
c

        di1=i+1
        call gammanew_eval_log_extend(di1,dlifac)

        dnp1=nn+ip/2.0d0+1
        call gammanew_eval_log_extend(dnp1,dlnpfac)

        dinp1=i+nn+ip/2.0d0+1
        call gammanew_eval_log_extend(dinp1,dlinpfac)

        dexp=dlinpfac-dlifac-dlnpfac

 1234   continue


c
c         create an array, dlogs, of logs of absolute values of terms
c

        dlogs(i+1)=log((2*(2*i+nn+ip/2.0d0+1))**(1/2.0d0)
     1     *abs(vect(i+1)))+dexp

        if (dlogmax .lt. dlogs(i+1)) then
          dlogmax=dlogs(i+1)
        endif

        dsigns(i+1)=sign(1.0d0,(-1)**i*vect(i+1))

c
c        . . . check if relative precision is achieved
c
        if (dlogs(i+1) .lt. log(dmach0/1000)+dlogmax) then
          nsum=i+1
          goto 110
        endif

 100    continue
 110    continue

c
c         now compute log(a_1+...+a_i)
c         given an array of log(abs(a_i)) and signs of a_i
c         

c
c        . . .  find alphais and nis
c
        do i=1,nsum
        nis(i)=dlogs(i)
        alphas(i)=exp(dlogs(i)-nis(i))
        enddo

c
c        . . . compute nmax
c

        nmax=dlogs(1)

        do i=1,nsum
        if (dlogs(i) .gt. nmax) then
          nmax=dlogs(i)
        endif
        enddo

c
c        . . . compute sum to be logged
c
        dsum=0
        do i=1,nsum
        dexp=nis(i)-nmax
        dsum=dsum+dsigns(i)*alphas(i)*exp(dexp)
        enddo

c
c        . . . save sign of dsum and take log of abs 
c
        nsign=sign(1.0d0,dsum)
        dlog=nmax+log(abs(dsum))

        return
        end
c
c
c
c
c
        subroutine dprol_ln_all_eigen_eval(nvects,nn,ip,c,w,nmat,
     1     rlams,vects)
        implicit real *8 (a-h,o-z)
        real*8 rlams(1),vects(1),w(1)

c
c
c         find the first nvects eigenvalues and eigenvectors 
c         of the symmetric tridiagonal matrix corresponding to 
c         the differential operator l using sturm bisection
c
c        w should be 14*nmat+220
c 
c        allocate memory
c

        call dprol_mat_size_evalq(nvects,c,nmat)
        n=nmat+1

ccc        call prinf('n*',n,1)

        iw2=1
        lw2=12*n+200
c
        ibis=iw2+lw2
        lbis=n+10
c
        icis=ibis+lbis
        lcis=n+10


        call dprol_ln_all_eigen_eval0(nvects,nn,ip,c,w(iw2),w(ibis),
     1     w(icis),nmat,rlams,vects)


        return
        end
c
c
c
c
c
        subroutine dprol_ln_all_eigen_eval0(nvects,nn,ip,c,w,bis,cis,
     1     nmat,rlams,vects)
        implicit real *8 (a-h,o-z)
        real*8 w(1),bis(1),cis(1),rlams(1),vects(1)

c
c         find the first nvects eigenvalues and eigenvectors 
c         of the symmetric tridiagonal matrix corresponding to 
c         the differential operator l using sturm bisection
c

c
c        . . . determine size of matrix
c

        call dprol_mat_size_evalq(nvects,c,nmat)

c
c        create an array corresponding to the diagonal
c        and another corresponding to the superdiagonal
c

        do 1000 i=0,nmat-1

        bi=-c**2*(nn+ip/2.0d0)**2/(2*(2*i+nn+ip/2.0d0)
     1     *(2*i+nn+ip/2.0d0+2))-c**2/2.0d0-
     1      (nn+ip/2.0d0+2*i+1/2.0d0)*(nn+ip/2.0d0+2*i+3/2.0d0)

        ci=-c**2*(i+1)*(i+nn+ip/2.0d0+1)/(sqrt(2*i+nn+ip/2.0d0+1)
     1      *(2*i+nn+ip/2.0d0+2)*sqrt(2*i+nn+ip/2.0d0+3))

        bis(i+1)=bi
        cis(i+1)=ci

 1000   continue

        if (nn+ip/2.0d0 .eq. 0) then
          bis(1)=-c**2/2.0d0-(1/2.0d0)*(3/2.0d0)
        endif
        
        call qleigen_trid(ier,nmat,bis,cis,rlams,vects,nvects,w)

c
c        . . . each vect is only unique up to sign find the first component
c        with magnitude bigger than eps and make it positive
c

        do 1400 ijk=1,nvects

        do i=1,nmat
        if (abs(vects((ijk-1)*nmat+i)) .gt. 1.0d-10) goto 400
        enddo
 400    continue

        if (vects((ijk-1)*nmat+i) .gt. 0) goto 800
        do i=1,nmat
        vects((ijk-1)*nmat+i)=-1*vects((ijk-1)*nmat+i)
        enddo

 800    continue

 1400   continue

        return
        end
c
c
c
c
c
        subroutine dprol_mach_zero(dmach0)
        implicit real *8 (a-h,o-z)
c
c        check machine precision and return in dmach0
c
        dmach0=1.0d0
        do i=1,40
        if (1 .eq. 1+dmach0) goto 1100
        dmach0=dmach0/1.0d1
        enddo
 1100   continue

        return
        end
c
c
c
c
c
        subroutine dprol_ln_n_eigen_eval(neig,nn,ip,c,w,nmat,
     1     rlam,vect)
        implicit real *8 (a-h,o-z)
        real*8 vect(1),w(1)

c
c         finds the neig^th smallest eigenvalue of the tridiagonal
c         matrix of the differential operator l using sturm 
c         bisection
c
c        w should be 12*nmat+200
c 
c        allocate memory
c

        call dprol_mat_size_evalq(neig,c,nmat)
        n=nmat+1

        ibis=1
        lbis=n+2
c
        icis=ibis+lbis
        lcis=n+2
c
        iw2=icis+lcis
        liw2=10*n+100

        call dprol_ln_n_eigen_eval0(neig,n,nn,ip,c,w(ibis)
     1     ,w(icis),w(iw2),rlam,vect)

        return
        end
c
c
c
c
c
        subroutine dprol_ln_n_eigen_eval0(neig,n,nn,ip,c,bis,
     1     cis,w,val,vect)
        implicit real *8 (a-h,o-z)
        real*8 bis(1),cis(1),vect(1),w(1)

c
c         finds the neig^th largest eigenvalue of the tridiagonal
c         matrix of the differential operator l using strum 
c         bisection
c

        nmat=n

c
c        create an array corresponding to the diagonal
c        and another corresponding to the superdiagonal
c
        do 1000 i=0,n

        bi=-c**2*(nn+ip/2.0d0)**2/(2*(2*i+nn+ip/2.0d0)
     1     *(2*i+nn+ip/2.0d0+2))-c**2/2.0d0-
     1      (nn+ip/2.0d0+2*i+1/2.0d0)*(nn+ip/2.0d0+2*i+3/2.0d0)

        ci=-c**2*(i+1)*(i+nn+ip/2.0d0+1)/(sqrt(2*i+nn+ip/2.0d0+1)
     1      *(2*i+nn+ip/2.0d0+2)*sqrt(2*i+nn+ip/2.0d0+3))

          
        bis(i+1)=bi
        cis(i+1)=ci

 1000   continue

        if (nn+ip/2.0d0 .eq. 0) then
          bis(1)=0
          bis(1)=-c**2/2.0d0-(1/2.0d0)*(3/2.0d0)
        endif

c
c        . . . now find the neig^th eigenvalue of the tridiagonal
c        symmetric matrix

        chi=(nn+ip/2.0d0+2*neig+1/2.0d0)*(nn+ip/2.0d0+2*neig+3/2.0d0)
        x1=-2*chi
        x2=-chi/10

        call dprol_ln_neigen_bisect(n,x1,x2,bis,cis,neig,w,val)
        
c
c        . . . and now that we have the eigenvalue use inverse power
c        to find the eigenvector
c
        rlam=val
        call dprol_ln_neigenvector(bis,cis,rlam,nmat,w,vect)

c
c        . . . vect is only unique up to sign find the first component
c        with magnitude bigger than eps and make it positive
c

        do i=1,nmat
        if (abs(vect(i)) .gt. 1.0d-10) goto 8960
        enddo
 8960   continue

        if (vect(i) .gt. 0) goto 9000

        do i=1,nmat
        vect(i)=-1*vect(i)
        enddo

 9000   continue

        return
        end
c
c
c
c
c
        subroutine dprol_mat_size_evalq(neig,c,nmat)
        implicit real *8 (a-h,o-z)
        real*8 cc(5),cc1(5),cc2(5),cc3(5)

c
c         find the necessary matrix size using precision ~1.0d-35
c         we use least squares to approximate delta_c where
c         delta_c=max{nmat-neig}, where nmat is the necessary matrix size.
c         that is, delta_c is the maximum number of terms needed above neig
c

c       a= 0.10000000D+03     b= 0.10000000D+04
        data cc1/
     1   0.25005641076886590000000000000000D+00,
     2   0.83543480388041740000000000000000D+02,
     3  -0.95699478083523120000000000000000D+04,
     4   0.96296736600004010000000000000000D+06,
     5  -0.38146709588571820000000000000000D+08/

c       a= 0.10000000D+04     b= 0.50000000D+04
        data cc2/
     1   0.23037791447056990000000000000000D+00,
     2   0.15601195563657330000000000000000D+03,
     3  -0.13585713602066420000000000000000D+06,
     4   0.12217763525226660000000000000000D+09,
     5  -0.51723229111680140000000000000000D+11/

c       a= 0.50000000D+04     b= 0.10000000D+05
        data cc3/
     1   0.22042720148693400000000000000000D+00,
     2   0.39031768830031520000000000000000D+03,
     3  -0.23335972283350120000000000000000D+07,
     4   0.97067372700240520000000000000000D+10,
     5  -0.16133415211940700000000000000000D+14/

c
c       find the appropriate expansion
c
        x=c
        m=5

        if(x .lt. 100) then                
            delta_c=c+100
            nmat=neig+delta_c
            return
        endif
c
        if(x .lt. 1000) then                
            call dprol_copy(cc1,cc,m)
            goto 1400
        endif
c
        if(x .lt. 5000) then                
            call dprol_copy(cc2,cc,m)
            goto 1400
        endif
c
        call dprol_copy(cc3,cc,m)
c
 1400 continue


        delta_c=0

        do 1600 i=1,m
        delta_c=delta_c+cc(i)*c**(2-i)
 1600 continue

        delta_c=delta_c+100

        nmat=delta_c+neig

        return
        end
c
c
c
c
c
        subroutine dprol_copy(x,y,n)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1)
c
ccc        call prin2('x*',x,n)
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
ccc        call prin2('y*',y,n)

        return
        end
c
c
c
c
c
        subroutine dprol_ln_neigen_bisect(n,x1,x2,bis,cis,neig,w,val)
        implicit real *8 (a-h,o-z)
        real*8 bis(1),cis(1),w(1)

c
c        . . . do sturm bisection to find the nieg^th eigenvalue
c

c 
c        allocate memory, w should be at least n+10 locations
c
c
        ipols=1
        lpols=n+10

        call dprol_ln_neigen_bisect0(n,x1,x2,bis,cis,neig,w(ipols),val)

        return
        end
c
c
c
c
c
        subroutine dprol_ln_neigen_bisect0(n,x1,x2,bis,cis,neig,
     1     pols,val)
        implicit real *8 (a-h,o-z)
        real*8 bis(1),cis(1),pols(1)

c
c        do sturm bisection to find the nieg^th eigenvalue
c

        nmat=n+1

c
c        first make sure x1 is small enough
c

        do 100 i=1,100

        call dprol_charpol_sturm(bis,cis,nmat,x1,pols,nsturm1)
        f1=nsturm1+neig-1/2.0d0-nmat

        if (f1 .lt. 0) then
          goto 110
        endif

        x1=x1-abs(x1)

 100    continue
 110    continue

c
c        . . . and that x2 is large enough
c

        do 200 i=1,100

        call dprol_charpol_sturm(bis,cis,nmat,x2,pols,nsturm1)
        f1=nsturm1+neig-1/2.0d0-nmat

        if (f1 .gt. 0) then
          goto 210
        endif

        x2=x2+100

 200    continue
 210    continue

c
c        . . . and then start bisection
c
        do 2000 i=1,150
        
        call dprol_charpol_sturm(bis,cis,nmat,x1,pols,nsturm1)
        f1=nsturm1+neig-1/2.0d0-nmat

        call dprol_charpol_sturm(bis,cis,nmat,x2,pols,nsturm2)
        f2=nsturm2+neig-1/2.0d0-nmat

        x3=(x1+x2)/2

        call dprol_charpol_sturm(bis,cis,nmat,x3,pols,nsturm3)
        f3=nsturm3+neig-1/2.0d0-nmat

        if (f3 .gt. 0) then
          x2=x3
          f2=f3
        endif

        if (f3 .le. 0) then
          x1=x3
          f1=f3
        endif

 2000   continue

        val=(x1+x2)/2

        return
        end
c 
c 
c 
c 
c 
        subroutine dprol_charpol_sturm(a,b,n,x,pols,nsturm)
        implicit real *8 (a-h,o-z)
        real *8 a(*),b(*),pols(*)
c 
c        calculate the sturm sequence
c 
        call dprol_sturmpols_eval(a,b,n,x,pols)
c 
c        calculate the sturm index
c 
        nsturm=0
        do i=1,n
        dd=pols(i+1)*pols(i)
        if(dd .le. 0) nsturm=nsturm+1
        if(pols(i+1) .eq. 0) nsturm=nsturm-1
        enddo
  
        return
        end
c
c 
c 
c 
c 
        subroutine dprol_ln_neigenvector(diag,sub,rlam,n,w,vect)
        implicit real *8 (a-h,o-z)
        real *8 diag(1),sub(1),vect(1),w(1)

c 
c        allocate memory
c
        i0=0

        irhs=1
        lrhs=n+2+10

        iww=irhs+lrhs
        lww=6*n+10

        idiag1=iww+lww
        ldiag1=n+2+10

        isub1=idiag1+ldiag1
        lsub1=n+2+10

        call dprol_ln_neigenvector0(diag,sub,rlam,w(irhs),w(iww),
     1     w(idiag1),w(isub1),n,vect)

        return
        end
c
c 
c 
c 
c 
        subroutine dprol_ln_neigenvector0(diag,sub,rlam,rhs,ww,
     1     diag1,sub1,n,vect)
        implicit real *8 (a-h,o-z)
        real *8 diag(1),sub(1),vect(1),rhs(1),ww(1),
     1      diag1(1),sub1(1),a(n,n),uout(100 000)

c
c         use inverse power to find the eigenvector corresponding 
c         to eigenvalue rlam of the n x n symmetric tridiagonal 
c         matrix corresponding to differential operator l
c

        call dprol_mach_zero(dmach0)
        delta=dmach0*100

c
c       . . . shift
c

        call dprol_copy(diag,diag1,n)
        call dprol_copy(sub,sub1,n-1)

        d=rlam

        do i=1,n
        diag1(i)=diag(i)-d
        enddo

c
c       . . .  construct the shifted inverse 
c

        call qleigen_fact(diag1,sub,n,ww,delta)
        call qleigen_rand(n,rhs)

        do 3000 ij=1,20
c
c       . . . apply the inverse
c
        call qleigen_comp_solv(rhs,n,ww)

c
c       . . . normalize
c
        call qleigen_rscap(rhs,rhs,n,d)

        d=1/sqrt(d)
        do 2600 i=1,n
        rhs(i)=rhs(i)*d
 2600 continue

 3000 continue

        do i=1,n
        vect(i)=rhs(i)
        enddo

ccc        call prin2('and vects=*',vects,n*n)
ccc        call prin2('vect in subroutine *',vect(nmat-10500),1000)

c
c        . . . vect is only unique up to sign find the first component
c        with magnitude bigger than eps and make it positive
c
        nmat=n
        do i=1,nmat
        if (abs(vect(i)) .gt. 1.0d-10) goto 8960
        enddo
 8960   continue

        if (vect(i) .gt. 0) goto 9010

        do i=1,nmat
        vect(i)=-1*vect(i)
        enddo

 9010   continue

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar_derivs(r,n,nn,ip,t0s,t1s)
        implicit real *8 (a-h,o-z)
        real*8 t0s(1),t1s(1)
c
c        This subroutine evaluatues the t-bar polynomials
c       _            _                   _
c       T_{nn,0}(r), T_{nn,1}(r), .... , T_{nn,n}(r)             (1)
c
c        and their derivatives

        call dprol_zern_pol_derivs(r,n,nn,ip,t0s,t1s)

        do 100 i=1,n+1
        yi0=r**((ip+1)/2.0d0)*t0s(i)

        yi1=r**((ip+1)/2.0d0)*t1s(i)+t0s(i)*((ip+1)/2.0d0)
     1     *r**((ip-1)/2.0d0)

        t0s(i)=yi0*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        t1s(i)=yi1*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))

 100    continue

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar_derivs2(r,n,nn,ip,t0s,t1s,t2s)
        implicit real *8 (a-h,o-z)
        real*8 t0s(1),t1s(1),t2s(1)
c
c        This subroutine evaluatues the t-bar polynomials
c       _            _                   _
c       T_{nn,0}(r), T_{nn,1}(r), .... , T_{nn,n}(r)             (1)
c
c        and their first and second derivatives

        call dprol_zern_pol_derivs2(r,n,nn,ip,t0s,t1s,t2s)

        do i=1,n+1
        yi0=r**((ip+1)/2.0d0)*t0s(i)

        yi1=r**((ip+1)/2.0d0)*t1s(i)+t0s(i)*((ip+1)/2.0d0)
     1     *r**((ip-1)/2.0d0)

        yi2=r**((ip+1)/2.0d0)*t2s(i)+t1s(i)*((ip+1)/2.0d0)*
     1     r**((ip-1)/2.0d0)
     1     +t0s(i)*(ip+1)*(ip-1)/4.0d0*r**((ip-3)/2.0d0)
     1     +((ip+1)/2.0d0)*r**((ip-1)/2.0d0)*t1s(i)

        t0s(i)=yi0*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        t1s(i)=yi1*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        t2s(i)=yi2*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        enddo

        return
        end
c 
c 
c 
c 
c 
        subroutine dprol_sturmpols_eval(a,b,n,x,pols)
        implicit real *8 (a-h,o-z)
        real*8 a(*),b(*),pols(*)
c
c        evaluate the polynomials in the sturm bisection
c        algorithm (see, e.g., theorem 2.22 of PSWFs of order
c        zero by Osipov et al. 
c
c        a - the diagonal of the tridiagonal matrix
c        b - the off-diagonal  
c        n - matrix is n x n 
c        x - argument at which to evaluate the polynomials
c


        pols(1)=1.0d0
        pols(2)=a(1) - x

c 
c        start recursion
c 
        do i=2,n

        pols(i+1) = (a(i)-x)*pols(i)-b(i-1)**2*pols(i - 1)

c
c        only signs matter, so scale to avoid under/overflow
c
        d=pols(i+1)**2+pols(i)**2

        if(d .gt. 1.0d12) then
          pols(i + 1) = pols(i + 1) * 1.0d-24
          pols(i) = pols(i) * 1.0d-24
        endif

        if(d .lt. 1.0d-12) then
          pols(i + 1) = pols(i + 1) * 1.0d24
          pols(i) = pols(i) * 1.0d24
        endif

        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_zern_pol_derivs(x,n,nn,ip,vals0,vals1)
        implicit real *8 (a-h,o-z)
        dimension vals0(1),vals1(1)
c
c        This subroutine evaluatues the (radial) Zernike polynomial
c
c       R_{nn,0}(r), R_{nn,1}(r), .... , R_{nn,n}(r)             (1)
c
c        where 0<r<1, for dimension ip+2.
c
        r=1-2*x**2
        a=nn+ip/2.0d0
        call dprol_jac_pol_derivs(r,n,a,vals0,vals1)

        do i=1,n+1
        vals1(i)=(-1)**(i+1)*(x**nn*vals1(i)*(-1)*4*x
     1                +vals0(i)*nn*x**(nn-1.0d0))
        vals0(i)=(-1)**(i+1)*x**nn*vals0(i)
        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_zern_normd_derivs(x,n,nn,ip,vals0,vals1)
        implicit real *8 (a-h,o-z)
        dimension vals0(1),vals1(1)
c
c        This subroutine evaluatues the (radial) Zernike polynomials
c
c       R_{nn,0}(r), R_{nn,1}(r), .... , R_{nn,n}(r)       
c
c        where 0<r<1, for dimension ip+2.
c
        r=1-2*x**2
        a=nn+ip/2.0d0

        call dprol_jac_pol_derivs(r,n,a,vals0,vals1)

        do i=1,n+1

        vals1(i)=(-1)**(i+1)*(x**nn*vals1(i)*(-1)*4*x
     1                +vals0(i)*nn*x**(nn-1.0d0))

        vals0(i)=(-1)**(i+1)*x**nn*vals0(i)

        vals1(i)=vals1(i)*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        vals0(i)=vals0(i)*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))

        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_zern_pol_derivs2(x,n,nn,ip,vals0,vals1,vals2)
        implicit real *8 (a-h,o-z)
        dimension vals0(1),vals1(1),vals2(1)
c
c        This subroutine evaluatues the (radial) Zernike polynomial
c
c       R_{nn,0}(r), R_{nn,1}(r), .... , R_{nn,n}(r)             (1)
c
c        where 0<r<1, for dimension ip+2.
c

        r=1-2*x**2
        a=nn+ip/2.0d0
        call dprol_jac_pol_derivs2(r,n,a,vals0,vals1,vals2)

        do i=1,n+1

        vals2(i)=(-1)**(i+1)*(16*x**(nn+2)*vals2(i)+(-4*(nn+1)*x**nn)
     1     *vals1(i)+vals0(i)*nn*(nn-1)*x**(nn-2)+nn*x**(nn-1)*vals1(i)
     1     *(-4*x))

        vals1(i)=(-1)**(i+1)*(x**nn*vals1(i)*(-4)*x
     1                +vals0(i)*nn*x**(nn-1.0d0))

        vals0(i)=(-1)**(i+1)*x**nn*vals0(i)

        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_zern_normd_derivs2(x,n,nn,ip,vals0,vals1,
     1     vals2)
        implicit real *8 (a-h,o-z)
        dimension vals0(1),vals1(1),vals2(1)
c
c        This subroutine evaluatues the normed (radial) Zernike 
c        polynomial
c
c       R_{nn,0}(r), R_{nn,1}(r), .... , R_{nn,n}(r)             (1)
c
c        where 0<r<1, for dimension ip+2.
c

        r=1-2*x**2
        a=nn+ip/2.0d0

        call dprol_jac_pol_derivs2(r,n,a,vals0,vals1,vals2)

        do i=1,n+1

        vals2(i)=(-1)**(i+1)*(16*x**(nn+2)*vals2(i)+(-4*(nn+1)*x**nn)
     1     *vals1(i)+vals0(i)*nn*(nn-1)*x**(nn-2)+nn*x**(nn-1)*vals1(i)
     1     *(-4*x))
        vals1(i)=(-1)**(i+1)*(x**nn*vals1(i)*(-4)*x
     1                +vals0(i)*nn*x**(nn-1.0d0))
        vals0(i)=(-1)**(i+1)*x**nn*vals0(i)

        vals2(i)=vals2(i)*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        vals1(i)=vals1(i)*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        vals0(i)=vals0(i)*sqrt(2*(2*(i-1)+nn+ip/2.0d0+1))
        enddo


        return
        end
c
c
c
c
c
        subroutine dprol_jac_pol_derivs(x,n,a,vals0,vals1)
        implicit real *8 (a-h,o-z)
        real*8 vals0(1),vals1(1)
c
c        This subroutine computes the Jacobi polynomials 
c
c          P^(a,0)_0(x), P^(a,0)_1(x), ... , P^(a,0)_n(x),
c
c        and their derivatives and stores the result in vals0
c        and vals1.
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

        do i=1,n-1

        a1i = 2.0d0*(i+1)*(i+a+1)*(2*i+a)
        a2i = (2.0d0*i+a+1)*a**2
        a3i = (2.0d0*i+a)*(2.0d0*i+a+1)*(2.0d0*i+a+2)
        a4i = 2.0d0*(i+a)*i*(2*i+a+2)

        f2 = ((a2i+a3i*x)*f1 - a4i*f0)/a1i
        vals0(i+2)=f2

        g2 = ((a2i+a3i*x)*g1 +a3i*f1- a4i*g0)/a1i
        vals1(i+2)=g2

        f0=f1
        f1=f2

        g0=g1
        g1=g2

        enddo

        return
        end
c
c
c
c
c
        subroutine dprol_jac_pol_derivs2(x,n,a,vals0,vals1,vals2)
        implicit real *8 (a-h,o-z)
        real*8 vals0(1),vals1(1),vals2(1)
c
c        This subroutine computes Jacobi polynomials 
c
c         P^(a,0)_0(x), P^(a,0)_1(x), ... , P^(a,0)_n(x),
c
c        and their first and second derivatives 
c        and stores the results in vals0, vals1, vals2.
c
        vals0(1)=1
        vals0(2)=0.5d0*(a+(a+2)*x)

        vals1(1)=0
        vals1(2)=0.5d0*(a+2)
        
        vals2(1)=0
        vals2(2)=0
        
c
c        . . . do jacobi polynomial recurrence
c
        f0 = 1
        f1 = 0.5d0*(a+(a+2)*x)

        g0 = 0
        g1 = 0.5d0*(a+2)

        h0=0
        h1=0

        do 200 i=1,n-1

        a1i = 2.0d0*(i+1)*(i+a+1)*(2*i+a)
        a2i = (2.0d0*i+a+1)*a**2
        a3i = (2.0d0*i+a)*(2.0d0*i+a+1)*(2.0d0*i+a+2)
        a4i = 2.0d0*(i+a)*i*(2*i+a+2)

        f2 = ((a2i+a3i*x)*f1 - a4i*f0)/a1i
        vals0(i+2)=f2

        g2 = ((a2i+a3i*x)*g1 +a3i*f1- a4i*g0)/a1i
        vals1(i+2)=g2

        h2=((a2i+a3i*x)*h1 +2*a3i*g1- a4i*h0)/a1i
        vals2(i+2)=h2

        f0=f1
        f1=f2

        g0=g1
        g1=g2

        h0=h1
        h1=h2

200     continue

        return
        end
c
c
c
c
c
        subroutine dprol_loggam_tail(x,k,val)
        implicit real*8 (a-h,o-z)
        real*8 dns(9),dds(9)

c
c        evaluates the tail of stirlings formula for ln(Gamma(dm)) 
c        up to k terms

        data dns/
     1        1.0D+00,
     2        -1.0D+00,
     3        1.0D+00,
     4        -1.0D+00,
     5        1.0D+00,
     6        691.0D+00,
     7        1.0D+00,
     8        -3617.0D+00,
     9        43867.0D+00/

        data dds/
     1        12.0D+00,
     2        360.0D+00,
     3        1260.0D+00,
     4        1680.0D+00,
     5        1188.0D+00,
     6        360360.0D+00,
     7        156.0D+00,
     8        122400.0D+00,
     9        244188.0D+00/

        val=0
        do 2000 i=1,k

        dinc=dns(i)/(x**(2*i-1)*dds(i))
        val=val+dinc

 2000   continue

        return
        end
c
c
c
c
c
        subroutine dprol_zern_normd_pol(r,n,nn,ip,val)
        implicit real *8 (a-h,o-z)
c
c        This subroutine evaluatues the normalized (radial) Zernike
c        polynomial
c          _________________          
c        \/ 2*(2*n+nn+p/2+1) * R_{nn,n}(r),
c
c        where 0<r<1, where the resulting function has L2 norm 1 with
c        weight function x^{p+1}.
c   

        call dprol_zern_pol(r,n,nn,ip,val2)

        val=sqrt(2.0d0*(2*n+nn+ip/2.0d0+1))*val2

        return
        end
c
c
c
c
c
        subroutine dprol_zern_pol(r,n,nn,ip,val)
        implicit real *8 (a-h,o-z)
c
c        This subroutine evaluatues the (radial) Zernike polynomial
c
c       R_{nn,n}(r),                                             (1)
c
c        where 0<r<1, for dimension ip+2.
c
        
        x=1-2*r**2
        a=nn+ip/2.0d0

        call dprol_jac_pol(x,n,a,val2)
        val=(-1)**n*r**nn*val2

        return
        end
c
c
c
c
c
        subroutine dprol_jac_pol(x,n,a,val)
        implicit real *8 (a-h,o-z)

c
c        This subroutine computes the Jacobi polynomial P^(a,0)_n(x).
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

        do i=1,n-1

        a1i = 2.0d0*(i+1)*(i+a+1)*(2*i+a)
        a2i = (2.0d0*i+a+1)*a**2
        a3i = (2.0d0*i+a)*(2.0d0*i+a+1)*(2.0d0*i+a+2)
        a4i = 2.0d0*(i+a)*i*(2*i+a+2)

        f2 = ((a2i+a3i*x)*f1 - a4i*f0)/a1i

        f0=f1
        f1=f2
        enddo

        val = f2

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar(r,n,nn,ip,t0)
        implicit real *8 (a-h,o-z)

c
c         evaluate the normalized tn polynomial
c                                     _
c           T_{nn,n}(r)=r**((ip+1)/2)*R_{nn,n}(r)
c

        call dprol_zern_pol_deriv2(r,n,nn,ip,val0,val1,val2)
        y0=r**((ip+1)/2.0d0)*val0
        t0=y0*sqrt(2*(2*n+nn+ip/2.0d0+1))

        return
        end
c
c
c
c
c
        subroutine dprol_t_deriv(r,n,nn,ip,t0,t1)
        implicit real *8 (a-h,o-z)

c
c        tn polynomial
c                                     
c           T_{nn,n}(r)=r**((ip+1)/2)*R_{nn,n}(r)
c
        call dprol_zern_pol_deriv2(r,n,nn,ip,val0,val1,val2)
        t0=r**((ip+1)/2.0d0)*val0
        t1=r**((ip+1)/2.0d0)*val1+val0*((ip+1)/2.0d0)*r**((ip-1)/2.0d0)

        return
        end
c
c
c
c
c
        subroutine dprol_tbar_deriv(r,n,nn,ip,t0,t1)
        implicit real *8 (a-h,o-z)

c
c         evaluate the normalized tn polynomial
c                                     _
c           T_{nn,n}(r)=r**((ip+1)/2)*R_{nn,n}(r)
c

        call dprol_zern_pol_deriv2(r,n,nn,ip,val0,val1,val2)

        t0=r**((ip+1)/2.0d0)*val0
        t0=t0*sqrt(2*(2*n+nn+ip/2.0d0+1))

        t1=r**((ip+1)/2.0d0)*val1+val0*((ip+1)/2.0d0)*r**((ip-1)/2.0d0)
        t1=t1*sqrt(2*(2*n+nn+ip/2.0d0+1))

        return
        end
c
c
c
c
c
        subroutine dprol_zern_pol_deriv2(r,n,nn,ip,val0,val1,val2)
        implicit real *8 (a-h,o-z)
c
c
c        This subroutine evaluatues the Zernike polynomial
c        R_{nn,n}, along with its first and second
c        derivatives. They are stored in val0, val1, and val2
c        respectively
c
        x=1-2*r**2
        a=nn+ip/2.0d0
        call dprol_jac_pol_deriv2(x,n,a,val0j,val1j,val2j)

        val0=(-1)**n*r**nn*val0j
        val1=(-1)**n*(-4*r**(nn+1)*val1j+val0j*nn*r**(nn-1.0d0))
        val2=(-1)**n*(16*r**(nn+2)*val2j+(-4*(nn+1)*r**nn)*val1j
     1        +val0j*nn*(nn-1)*r**(nn-2)+nn*r**(nn-1)*val1j*(-4*r))

        return
        end
c
c
c
c
c
        subroutine dprol_t_bar_deriv2(r,n,nn,ip,val0,val1,val2)
        implicit real *8 (a-h,o-z)
c
c
c        This subroutine evaluatues the Zernike polynomial
c        R_{nn,n}, along with its first and second
c        derivatives. They are stored in val0, val1, and val2
c        respectively
c
        call dprol_zern_pol_deriv2(r,n,nn,ip,r0,r1,r2)

        val0=r**((ip+1)/2.0d0)*r0

        val1=r**((ip+1)/2.0d0)*r1+r0*((ip+1)/2.0d0)
     1     *r**((ip-1)/2.0d0)

        val2=r**((ip+1)/2.0d0)*r2+r1*((ip+1)/2.0d0)*
     1     r**((ip-1)/2.0d0)
     1     +r0*(ip+1)*(ip-1)/4.0d0*r**((ip-3)/2.0d0)
     1     +((ip+1)/2.0d0)*r**((ip-1)/2.0d0)*r1

        val0=val0*sqrt(2*(2*n+nn+ip/2.0d0+1))
        val1=val1*sqrt(2*(2*n+nn+ip/2.0d0+1))
        val2=val2*sqrt(2*(2*n+nn+ip/2.0d0+1))

        return
        end
c
c
c
c
c
        subroutine dprol_jac_pol_deriv2(x,n,a,val0,val1,val2)
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

        do i=1,n-1

        a1i = 2.0d0*(i+1)*(i+a+1)*(2*i+a)
        a2i = (2.0d0*i+a+1)*a**2
        a3i = (2.0d0*i+a)*(2.0d0*i+a+1)*(2.0d0*i+a+2)
        a4i = 2.0d0*(i+a)*i*(2*i+a+2)

        f2 = ((a2i+a3i*x)*f1 - a4i*f0)/a1i
        g2 = ((a2i+a3i*x)*g1 +a3i*f1- a4i*g0)/a1i
        h2=  ((a2i+a3i*x)*h1 +2*a3i*g1- a4i*h0)/a1i

        f0=f1
        f1=f2

        g0=g1
        g1=g2

        h0=h1
        h1=h2

        enddo

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
        subroutine dprol_ratio_eval(nn,ip,nmat,c1,c2,r21,w)
        implicit real *8 (a-h,o-z)
        real*8 w(*), c1(*), c2(*)
c
c        take the zernike expansion coefficients c1, c2 of 
c        two eigenfunctions and compute the ratio between the 
c        eigenvalues beta2 and beta1.
c
        itemp = 1
        idtemp = nmat+10
        call dprol_ratio_eval0(nn,ip,nmat,c1,c2,r21,w(itemp),w(idtemp))

        return
        end
c
c
c
c
c
        subroutine dprol_ratio_eval0(nn,ip,nmat,c1,c2,r21,temp,dtemp)
        implicit real *8 (a-h,o-z)
        real*8 c1(*), c2(*), temp(*), dtemp(*)
c
c        c1 - contains nmat zernike expansion coefficients of an 
c             eigenfunction phi_n for some n
c        c2 - contains nmat zernike expansion coefficients of an 
c             eigenfunction phi_{n+1} for some n
c
        call dprol_copy(c1,temp,nmat)
        call dprol_diffcoefs(nn,ip,temp,nmat,dtemp)

c
c        dtemp now contains coefficients of r*phi_n'
c
        rnum = 0.0d0
        do i=1,nmat
        rnum = rnum + dtemp(i)*c2(i)
        enddo

        call dprol_copy(c2,temp,nmat)
        call dprol_diffcoefs(nn,ip,temp,nmat,dtemp)
c
c        dtemp now contains coefficients of r*phi_{n+1}'
c
        rdenom = 0.0d0
        do i=1,nmat
        rdenom = rdenom + dtemp(i)*c1(i)
        enddo

        r21 = rnum/rdenom

        return
        end
c
c
c
c
c
        subroutine dprol_diffcoefs(nn,ip,d,nmat,dd)
        implicit real *8 (a-h,o-z)
        real*8 d(nmat), dd(nmat)
c
c        comptute the tbar expansion coeffients (dd) of the function
c        r Phi'(r) given the tbar expansion of Phi (d)
c
        a = nn+ip/2.0d0
        pp12 = (ip + 1.0)/2.0d0
c
c        convert expansion in \overline{R}_{nn,n} to expansion in 
c        R_{nn,n}
c
        do i=1,nmat
        d(i) = sqrt(2.0d0*(2*(i-1)+a+1))*d(i)
        enddo

c
c        correct for alternating signs and initialize expansion coefs
c        to zero
c
        do i=1,nmat
        d(i) = (-1)**(i+1)*d(i)
        dd(i) = 0.0d0
        enddo

c
c        go from n=nmat-2 to n=1
c
        do 100 i=1,nmat-2
        n = nmat-i-1
        at = 2.0d0*(n+a+1)*(2*n + a)
        bt = 2.0d0*a*(2*n + a + 1)
        ct = -2.0d0*n*(2*n + a + 2)
        a1 = (2.0d0*a + 4*n + 5)*(n + a + 1)*(2*n+a)
        b1 = a*(2.0d0*n + a + 1) - 2*(2*n+a)*(2*n+a+1)*(2*n+a+2)
        c1 = n*(2.0d0*a + 4*n - 1)*(2*n + a + 2)

        d(n+1) = d(n+1) - d(n+2) * bt / at
        d(n) = d(n) - d(n+2) * ct / at

        dd(n+2) = dd(n+2) + d(n+2)*(a1 - pp12*at)/at
        dd(n+1) = dd(n+1) + d(n+2)*(b1 - pp12*bt)/at
        dd(n) = dd(n) + d(n+2)*(c1-pp12*ct)/at

100     continue

c
c        now n=0
c
        n = 0
        at = 2.0d0*(n+a+1)*(2*n + a)
        bt = 2.0d0*a*(2*n + a + 1)
        ct = -2.0d0*n*(2*n + a + 2)
        a1 = (2.0d0*a + 4*n + 5)*(n + a + 1)*(2*n+a)
        b1 = a*(2.0d0*n + a + 1) - 2*(2*n+a)*(2*n+a+1)*(2*n+a+2)
        c1 = n*(2.0d0*a + 4*n - 1)*(2*n + a + 2)

        if (a .eq. 0d0) then
          d(n+1) = d(n+1) - d(n+2)
          dd(n+2) = dd(n+2) + d(n+2)*(5d0-(ip+1d0))/2d0
          dd(n+1) = dd(n+1) + d(n+2)*(-3d0-(ip+1d0))/2d0
        end if

        if (a .gt. 0d0) then
          d(n+1) = d(n+1) - d(n+2) * bt / at
          dd(n+2) = dd(n+2) + d(n+2)*(a1 - pp12*at)/at
          dd(n+1) = dd(n+1) + d(n+2)*(b1 - pp12*bt)/at
        end if

        dd(1)=dd(1)+nn*d(1)

c
c        convert expansion in R_{nn,n} back to expansion in
c        \overline{R}_{nn,n}
c
        do i=1,nmat
        dd(i) = dd(i)/sqrt(2*(2*(i-1)+a+1))
        dd(i) = (-1)**(i+1)*dd(i)
        enddo

        return
        end
