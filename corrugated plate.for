subroutine uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, props, nprops, coords, &
               mcrd, nnode, u, du, v, a, jtype, time, dtime, kstep, kinc, jelem, params, &
               ndload, jdltyp, adlmag, predef, npredf, lflags, mlvarx, ddlmag, mdload, &
               pnewdt, jprops, njprop, period)

    Include 'ABA_PARAM.INC'

    Dimension rhs(mlvarx, *), amatrx(ndofel, ndofel), props(*), svars(*), energy(8), coords(mcrd, nnode), &
        u(ndofel), du(mlvarx, *), v(ndofel), a(ndofel), time(2), params(*), jdltyp(mdload, *), adlmag(mdload, *), &
        ddlmag(mdload, *), predef(2, npredf, nnode), lflags(*), jprops(*)

    !ndofel -total dofs of element
    !nnode  -total nodes of element
    !mcrd -number of coordinates to descripe a node, at least 3
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    integer :: i,j,k,nwall
    integer,PARAMETER::nip=4    !number of integrating point
    integer,PARAMETER::ndim=2    !number of dimension
    integer,PARAMETER::nst=14    !number of stress term
    integer,PARAMETER::nsec=500    !number of sections in trapz integration
    CHARACTER(LEN=30)::element='quadrilateral'
    CHARACTER(LEN=30)::boardtype='',scheme='',honeycomb_shape=''
    REAL(iwp)::zero=0.0_iwp, det
    real(iwp):: points(nip,ndim), weights(nip), der(ndim,nnode),jac(ndim,ndim), &
        bee(nst,ndofel), dee(nst,nst), deriv(ndim,nnode), fun(nnode), coord(nnode,ndim),&
        utran(ndofel,ndofel), md_vector(3), normal_vector(3)
    real(iwp),ALLOCATABLE::faceprop(:,:),coreprop(:,:)
    integer::first_call=1
    save first_call,dee
    !---------------------------------------------------------------------------
    !debug breakpoint
    !WRITE(*,*) "Debug interruption, enter an integer to continue"
    !read(*,*) tempRead
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !read material properties form external file, call only once assuming only one kind of paperboard is used
    if (first_call==1) then
        IF (ALLOCATED(faceprop)) THEN
            DEALLOCATE(faceprop)
        END IF
        IF (ALLOCATED(coreprop)) THEN
            DEALLOCATE(coreprop)
        END IF
        CALL readfile(boardtype,nwall,faceprop,coreprop,scheme,honeycomb_shape)
        !deemat for two paperboards
        if (boardtype=='corrugated') then
            CALL deemat_corrugated(dee,faceprop,coreprop,scheme,nsec)
        else if (boardtype=='honeycomb')  then
            CALL deemat_honeycomb(dee,faceprop,coreprop,scheme)
        else if (boardtype=='isotropic') then
            CALL deemat_isotropic(dee,props(1),props(2),props(3))
        else
            WRITE(*,*) 'UEL:material property input error' 
            CALL XIT
        end if
        first_call=1
    end if
    !-----------------------------------------------------------------------------
    !check node number and dofel
    if (nnode /=8) THEN
        WRITE(*,*) 'Element node number error'
        CALL XIT
    end if
    if (ndofel /=56) THEN
        WRITE(*,*) 'Dof per node error'
        CALL XIT
    end if
    !---------------------------------------------------------------------------
    !element orientation vector, should be in the plane of the element
    if (boardtype/='isotropic') then
        md_vector=props(1:3)
        normal_vector=props(4:6)
    end if
    !-----------------------------------------------------------------------------
    !transfer 3d-coordinate to 2d plane
    CALL planetrans(coords,coord,utran,md_vector,normal_vector,boardtype)
    !-----------------------------------------------------------------------------
    !debug break
    !WRITE(*,*) 'program exit because of debugging'
    !CALL XIT
    !====================================================================
    !main code
    CALL sample(element,points,weights)
    amatrx=zero
    DO i=1,nip
        CALL shape_der_m(der,points,i) 
        jac=MATMUL(der,coord) 
        det=determinant(jac) 
        CALL invert(jac) 
        deriv=MATMUL(jac,der)
        CALL shape_fun_m(fun,points,i)
        CALL beemat_wHSDT(bee,deriv,fun)
        amatrx=amatrx+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
    END DO
    amatrx=MATMUL(MATMUL(TRANSPOSE(utran),amatrx),utran)
    !rhs=-K*U
    DO i = 1, ndofel
        rhs(i,1)=zero
        DO j = 1, ndofel
            rhs(i,1) = rhs(i,1) - amatrx(i, j) * u(j)
        END DO
    END DO
    !write(*,*) mlvarx
    !===================================================================
    !debug break
    !WRITE(*,*) 'program exit because of debugging'
    !CALL XIT
    RETURN
    !==================================================
    contains
        SUBROUTINE planetrans(coords,coord,utran,md_vector,normal_vector,boardtype)
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            CHARACTER(LEN=30)::boardtype
            real(iwp),INTENT(IN)::coords(:,:),md_vector(:),normal_vector(:)
            real(iwp),INTENT(OUT)::coord(:,:),utran(:,:)
            real(iwp)::md(3), nv(3), coortran(3,3), cd(3)
            real(iwp),PARAMETER::thres=0.01_iwp,one=1.0_iwp,zero=0.0_iwp
            real(iwp),allocatable::temp(:,:)
            integer:: nnode,i
            nnode=UBOUND(coords,2)
            allocate(temp(nnode,3))
            coortran=zero
            if (boardtype/='isotropic') then
                !transfer 3d-coordinate to 2d plane
                md=md_vector/norm(md_vector)
                nv=normal_vector/norm(normal_vector)
                cd=cross_product3(nv,md)
                coortran(:,1)=md
                coortran(:,2)=cd
                coortran(:,3)=nv
                call invert(coortran)
            else 
                coortran(1,1)=one
                coortran(2,2)=one
                coortran(3,3)=one
            end if
            temp=TRANSPOSE(MATMUL(coortran,coords))
            coord=temp(:,1:2)
            !---------------------------
            !displacement transfer matrix
            utran=zero
            do i=1,nnode !change with DOFs
                utran(7*i-6:7*i-4,7*i-6:7*i-4)=coortran
                utran(7*i-3:7*i-2,7*i-3:7*i-2)=coortran(1:2,1:2)
                utran(7*i-1,7*i-1)=one
                utran(7*i,7*i)=one
            end do
            DEALLOCATE(temp)
        END SUBROUTINE planetrans

        SUBROUTINE readfile(boardtype,nwall,faceprop,coreprop,scheme,honeycomb_shape)
            !get material property form file
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            CHARACTER(LEN=256)::filename,OUTDIR
            INTEGER::LENOUTDIR,nwall,ios
            CHARACTER(LEN=30)::boardtype,scheme,honeycomb_shape
            real(iwp),ALLOCATABLE::faceprop(:,:),coreprop(:,:)

            CALL GETOUTDIR(OUTDIR,LENOUTDIR)
            filename = trim(OUTDIR)//'\matprop.txt'
            open(unit=11, file=filename, status='OLD', action='read',FORM='formatted',iostat=ios)
            rewind(11)
            read(11,*) boardtype
                SELECT CASE(boardtype)
                    case('corrugated')
                        read(11,*) nwall
                        if (nwall>3) then
                            WRITE(*,*) 'maximum ply for corrugated board is 3'
                            CALL XIT
                        end if
                        read(11,*) scheme
                        ALLOCATE(faceprop(nwall+1,7),coreprop(nwall,9))
                        do i=1,nwall+1
                            read(11,*) faceprop(i,:)
                            if (ios /= 0) then
                                WRITE(*,*) 'Error reading faceprop data from file'
                                call xit
                            end if
                        end do
                        do i=1,nwall
                            read(11,*) coreprop(i,:)
                            if (ios /= 0) then
                                WRITE(*,*) 'Error reading faceprop data from file'
                                call xit
                            end if
                        end do      
                    case('honeycomb')
                        ALLOCATE(faceprop(2,7),coreprop(1,9))
                        read(11,*) honeycomb_shape
                        read(11,*) scheme
                        do i=1,2
                            read(11,*) faceprop(i,:)
                            if (ios /= 0) then
                                WRITE(*,*) 'Error reading faceprop data from file'
                                call xit
                            end if
                        end do               
                        if (honeycomb_shape=='regular') then
                            read(11,*) coreprop(1,:) 
                        else
                            WRITE(*,*) 'irregular honeycomb is developing'
                            call xit
                        end if   
                    case('isotropic')
                        RETURN
                    case DEFAULT
                        WRITE(*,*) 'readfile:material property input error' 
                        WRITE(*,*) boardtype
                        call xit
                end SELECT   
            close(11)
        END SUBROUTINE readfile

        SUBROUTINE sample(element,s,wt)
            !
            ! This subroutine returns the local coordinates and weighting coefficients
            ! of the integrating points.
            !
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(OUT)::s(:,:)
            REAL(iwp),INTENT(OUT),OPTIONAL::wt(:)
            CHARACTER(*),INTENT(IN)::element
            INTEGER::nip
            REAL(iwp)::root3,r15,w(3),v(9),b,c
            root3=1.0_iwp/SQRT(3.0_iwp)
            r15=0.2_iwp*SQRT(15.0_iwp)
            nip=UBOUND(s,1)
            w=(/5.0_iwp/9.0_iwp,8.0_iwp/9.0_iwp,5.0_iwp/9.0_iwp/)
            v=(/5.0_iwp/9.0_iwp*w,8.0_iwp/9.0_iwp*w,5.0_iwp/9.0_iwp*w/)
            SELECT CASE(element)
            CASE('line')
            SELECT CASE(nip)
            CASE(1)
                s(1,1)=0.0_iwp
                wt(1) =2.0_iwp
            CASE(2)
                s(1,1)=-0.577350269189626_iwp
                s(2,1)= 0.577350269189626_iwp
                wt(1) = 1.000000000000000_iwp
                wt(2) = 1.000000000000000_iwp
            CASE(3)
                s(1,1)=-0.774596669241484_iwp
                s(2,1)= 0.000000000000000_iwp
                s(3,1)= 0.774596669241484_iwp
                wt(1) = 0.555555555555556_iwp
                wt(2) = 0.888888888888889_iwp
                wt(3) = 0.555555555555556_iwp
            CASE(4)
                s(1,1)=-0.861136311594053_iwp
                s(2,1)=-0.339981043584856_iwp
                s(3,1)= 0.339981043584856_iwp
                s(4,1)= 0.861136311594053_iwp
                wt(1) = 0.347854845137454_iwp
                wt(2) = 0.652145154862546_iwp
                wt(3) = 0.652145154862546_iwp
                wt(4) = 0.347854845137454_iwp
            CASE(5)
                s(1,1)=-0.906179845938664_iwp
                s(2,1)=-0.538469310105683_iwp
                s(3,1)= 0.000000000000000_iwp
                s(4,1)= 0.538469310105683_iwp
                s(5,1)= 0.906179845938664_iwp
                wt(1) = 0.236926885056189_iwp
                wt(2) = 0.478628670499366_iwp
                wt(3) = 0.568888888888889_iwp
                wt(4) = 0.478628670499366_iwp
                wt(5) = 0.236926885056189_iwp
            CASE(6)
                s(1,1)=-0.932469514203152_iwp
                s(2,1)=-0.661209386466265_iwp
                s(3,1)=-0.238619186083197_iwp
                s(4,1)= 0.238619186083197_iwp
                s(5,1)= 0.661209386466265_iwp
                s(6,1)= 0.932469514203152_iwp
                wt(1) = 0.171324492379170_iwp
                wt(2) = 0.360761573048139_iwp
                wt(3) = 0.467913934572691_iwp
                wt(4) = 0.467913934572691_iwp
                wt(5) = 0.360761573048139_iwp
                wt(6) = 0.171324492379170_iwp
            CASE(7)
                s(1,1)=-0.9491079123427585245261897_iwp
                s(2,1)=-0.7415311855993944398638648_iwp
                s(3,1)=-0.4058451513773971669066064_iwp
                s(4,1)= 0.000000000000000_iwp
                s(5,1)= 0.4058451513773971669066064_iwp
                s(6,1)= 0.7415311855993944398638648_iwp
                s(7,1)= 0.9491079123427585245261897_iwp
                wt(1) = 0.1294849661688696932706114_iwp
                wt(2) = 0.2797053914892766679014678_iwp
                wt(3) = 0.3818300505051189449503698_iwp
                wt(4) = 0.4179591836734693877551020_iwp
                wt(5) = 0.3818300505051189449503698_iwp
                wt(6) = 0.2797053914892766679014678_iwp
                wt(7) = 0.1294849661688696932706114_iwp
            CASE(8)
                s(1,1)=-0.9602898564975362316835609_iwp
                s(2,1)=-0.7966664774136267395915539_iwp
                s(3,1)=-0.5255324099163289858177390_iwp
                s(4,1)=-0.1834346424956498049394761_iwp
                s(5,1)= 0.1834346424956498049394761_iwp
                s(6,1)= 0.5255324099163289858177390_iwp
                s(7,1)= 0.7966664774136267395915539_iwp
                s(8,1)= 0.9602898564975362316835609_iwp
                wt(1) = 0.1012285362903762591525314_iwp
                wt(2) = 0.2223810344533744705443560_iwp
                wt(3) = 0.3137066458778872873379622_iwp
                wt(4) = 0.3626837833783619829651504_iwp
                wt(5) = 0.3626837833783619829651504_iwp
                wt(6) = 0.3137066458778872873379622_iwp
                wt(7) = 0.2223810344533744705443560_iwp
                wt(8) = 0.1012285362903762591525314_iwp
            CASE(9)
                s(1,1)=-0.9681602395076260898355762_iwp
                s(2,1)=-0.8360311073266357942994298_iwp    
                s(3,1)=-0.6133714327005903973087020_iwp
                s(4,1)=-0.3242534234038089290385380_iwp    
                s(5,1)= 0.000000000000000_iwp                            
                s(6,1)= 0.3242534234038089290385380_iwp                            
                s(7,1)= 0.6133714327005903973087020_iwp                            
                s(8,1)= 0.8360311073266357942994298_iwp                            
                s(9,1)= 0.9681602395076260898355762_iwp                            
                wt(1) = 0.0812743883615744119718922_iwp                            
                wt(2) = 0.1806481606948574040584720_iwp                            
                wt(3) = 0.2606106964029354623187429_iwp                            
                wt(4) = 0.3123470770400028400686304_iwp                            
                wt(5) = 0.3302393550012597631645251_iwp                            
                wt(6) = 0.3123470770400028400686304_iwp                            
                wt(7) = 0.2606106964029354623187429_iwp                            
                wt(8) = 0.1806481606948574040584720_iwp                            
                wt(9) = 0.0812743883615744119718922_iwp                            
            CASE(10)
                s(1,1)=-0.9739065285171717200779640_iwp            
                s(2,1)=-0.8650633666889845107320967_iwp 
                s(3,1)=-0.6794095682990244062343274_iwp 
                s(4,1)=-0.4333953941292471907992659_iwp 
                s(5,1)=-0.1488743389816312108848260_iwp 
                s(6,1)= 0.1488743389816312108848260_iwp 
                s(7,1)= 0.4333953941292471907992659_iwp 
                s(8,1)= 0.6794095682990244062343274_iwp 
                s(9,1)= 0.8650633666889845107320967_iwp 
                s(10,1)= 0.9739065285171717200779640_iwp 
                wt(1) = 0.0666713443086881375935688_iwp                     
                wt(2) = 0.1494513491505805931457763_iwp                     
                wt(3) = 0.2190863625159820439955349_iwp                     
                wt(4) = 0.2692667193099963550912269_iwp                     
                wt(5) = 0.2955242247147528701738930_iwp                     
                wt(6) = 0.2955242247147528701738930_iwp                      
                wt(7) = 0.2692667193099963550912269_iwp                     
                wt(8) = 0.2190863625159820439955349_iwp                     
                wt(9) = 0.1494513491505805931457763_iwp                     
                wt(10) = 0.0666713443086881375935688_iwp                     
            CASE DEFAULT                              
                WRITE(*,*)"Wrong number of integrating points for a line"
            END SELECT
            CASE('triangle')
            SELECT CASE(nip)
            CASE(1)
                s(1,1)= 0.333333333333333_iwp
                s(1,2)= 0.333333333333333_iwp
                wt(1) = 0.500000000000000_iwp
            CASE(3)
                s(1,1)= 0.500000000000000_iwp
                s(1,2)= 0.500000000000000_iwp
                s(2,1)= 0.500000000000000_iwp
                s(2,2)= 0.000000000000000_iwp
                s(3,1)= 0.000000000000000_iwp
                s(3,2)= 0.500000000000000_iwp
                wt(1:3)=0.333333333333333_iwp
                wt=0.5_iwp*wt
            CASE(4)
                s(1,1)= 0.6_iwp
                s(1,2)= 0.2_iwp
                s(2,1)= 0.2_iwp
                s(2,2)= 0.6_iwp
                s(3,1)= 0.2_iwp
                s(3,2)= 0.2_iwp
                s(4,1)= 0.333333333333333_iwp
                s(4,2)= 0.333333333333333_iwp
                wt(1:3)= 0.520833333333333_iwp
                wt(4)=  -0.5625_iwp
                wt=0.5_iwp*wt
            CASE(6)
                s(1,1)= 0.816847572980459_iwp
                s(1,2)= 0.091576213509771_iwp
                s(2,1)= 0.091576213509771_iwp
                s(2,2)= 0.816847572980459_iwp
                s(3,1)= 0.091576213509771_iwp
                s(3,2)= 0.091576213509771_iwp
                s(4,1)= 0.108103018168070_iwp
                s(4,2)= 0.445948490915965_iwp
                s(5,1)= 0.445948490915965_iwp
                s(5,2)= 0.108103018168070_iwp
                s(6,1)= 0.445948490915965_iwp
                s(6,2)= 0.445948490915965_iwp
                wt(1:3)=0.109951743655322_iwp
                wt(4:6)=0.223381589678011_iwp
                wt=0.5_iwp*wt
            CASE(7)
                s(1,1)= 0.333333333333333_iwp
                s(1,2)= 0.333333333333333_iwp
                s(2,1)= 0.797426985353087_iwp
                s(2,2)= 0.101286507323456_iwp
                s(3,1)= 0.101286507323456_iwp
                s(3,2)= 0.797426985353087_iwp
                s(4,1)= 0.101286507323456_iwp
                s(4,2)= 0.101286507323456_iwp
                s(5,1)= 0.470142064105115_iwp
                s(5,2)= 0.059715871789770_iwp
                s(6,1)= 0.059715871789770_iwp
                s(6,2)= 0.470142064105115_iwp
                s(7,1)= 0.470142064105115_iwp
                s(7,2)= 0.470142064105115_iwp
                wt(1) = 0.225000000000000_iwp
                wt(2:4)=0.125939180544827_iwp
                wt(5:7)=0.132394152788506_iwp
                wt=0.5_iwp*wt
            CASE(12)
                s(1,1)= 0.873821971016996_iwp
                s(1,2)= 0.063089014491502_iwp
                s(2,1)= 0.063089014491502_iwp
                s(2,2)= 0.873821971016996_iwp
                s(3,1)= 0.063089014491502_iwp
                s(3,2)= 0.063089014491502_iwp
                s(4,1)= 0.501426509658179_iwp
                s(4,2)= 0.249286745170910_iwp
                s(5,1)= 0.249286745170910_iwp
                s(5,2)= 0.501426509658179_iwp
                s(6,1)= 0.249286745170910_iwp
                s(6,2)= 0.249286745170910_iwp
                s(7,1) =0.053145049844817_iwp
                s(7,2) =0.310352451033784_iwp
                s(8,1) =0.310352451033784_iwp
                s(8,2) =0.053145049844817_iwp
                s(9,1) =0.053145049844817_iwp
                s(9,2) =0.636502499121398_iwp
                s(10,1)=0.310352451033784_iwp
                s(10,2)=0.636502499121398_iwp
                s(11,1)=0.636502499121398_iwp
                s(11,2)=0.053145049844817_iwp
                s(12,1)=0.636502499121398_iwp
                s(12,2)=0.310352451033784_iwp
                wt(1:3)=0.050844906370207_iwp
                wt(4:6)=0.116786275726379_iwp
                wt(7:12)=0.082851075618374_iwp
                wt=0.5_iwp*wt
            CASE(16)
                s(1,1)=0.333333333333333_iwp
                s(1,2)=0.333333333333333_iwp
                s(2,1)=0.658861384496478_iwp
                s(2,2)=0.170569307751761_iwp
                s(3,1)=0.170569307751761_iwp
                s(3,2)=0.658861384496478_iwp
                s(4,1)=0.170569307751761_iwp
                s(4,2)=0.170569307751761_iwp
                s(5,1)=0.898905543365938_iwp
                s(5,2)=0.050547228317031_iwp
                s(6,1)=0.050547228317031_iwp
                s(6,2)=0.898905543365938_iwp
                s(7,1)=0.050547228317031_iwp
                s(7,2)=0.050547228317031_iwp
                s(8,1)=0.081414823414554_iwp
                s(8,2)=0.459292588292723_iwp
                s(9,1)=0.459292588292723_iwp
                s(9,2)=0.081414823414554_iwp
                s(10,1)=0.459292588292723_iwp
                s(10,2)=0.459292588292723_iwp
                s(11,1)=0.008394777409958_iwp
                s(11,2)=0.263112829634638_iwp
                s(12,1)=0.008394777409958_iwp
                s(12,2)=0.728492392955404_iwp
                s(13,1)=0.263112829634638_iwp
                s(13,2)=0.008394777409958_iwp
                s(14,1)=0.263112829634638_iwp
                s(14,2)=0.728492392955404_iwp
                s(15,1)=0.728492392955404_iwp
                s(15,2)=0.008394777409958_iwp
                s(16,1)=0.728492392955404_iwp
                s(16,2)=0.263112829634638_iwp
                wt(1)=0.144315607677787_iwp
                wt(2:4)=0.103217370534718_iwp
                wt(5:7)=0.032458497623198_iwp
                wt(8:10)=0.095091634267284_iwp
                wt(11:16)=0.027230314174435_iwp
                wt=0.5_iwp*wt
            CASE DEFAULT
                WRITE(*,*)"wrong number of integrating points for a triangle"
            END SELECT
            CASE('quadrilateral')
            SELECT CASE(nip)
            CASE(1)
                s(1,1)=0.0_iwp
                s(1,2)=0.0_iwp
                wt(1)=4.0_iwp
            CASE(4)
                s(1,1)=-root3
                s(1,2)= root3
                s(2,1)= root3
                s(2,2)= root3
                s(3,1)=-root3
                s(3,2)=-root3
                s(4,1)= root3
                s(4,2)=-root3
                wt=1.0_iwp
            CASE(9)
                s(1:7:3,1)=-r15
                s(2:8:3,1)=0.0_iwp
                s(3:9:3,1)=r15
                s(1:3,2)  =r15
                s(4:6,2)  =0.0_iwp
                s(7:9,2)  =-r15
                wt= v
            CASE(16)
                s(1:13:4,1)=-0.861136311594053_iwp
                s(2:14:4,1)=-0.339981043584856_iwp
                s(3:15:4,1)= 0.339981043584856_iwp
                s(4:16:4,1)= 0.861136311594053_iwp
                s(1:4,2)   = 0.861136311594053_iwp
                s(5:8,2)   = 0.339981043584856_iwp
                s(9:12,2)  =-0.339981043584856_iwp
                s(13:16,2) =-0.861136311594053_iwp
                wt(1)      = 0.121002993285602_iwp
                wt(4)      = wt(1)
                wt(13)     = wt(1)
                wt(16)     = wt(1)
                wt(2)      = 0.226851851851852_iwp
                wt(3)      = wt(2)
                wt(5)      = wt(2)
                wt(8)      = wt(2)
                wt(9)      = wt(2)
                wt(12)     = wt(2)
                wt(14)     = wt(2)
                wt(15)     = wt(2)
                wt(6)      = 0.425293303010694_iwp
                wt(7)      = wt(6)
                wt(10)     = wt(6)
                wt(11)     = wt(6)
            CASE(25)
                s(1:21:5,1)= 0.906179845938664_iwp
                s(2:22:5,1)= 0.538469310105683_iwp
                s(3:23:5,1)= 0.0_iwp
                s(4:24:5,1)=-0.538469310105683_iwp
                s(5:25:5,1)=-0.906179845938664_iwp
                s( 1: 5,2) = 0.906179845938664_iwp
                s( 6:10,2) = 0.538469310105683_iwp
                s(11:15,2) = 0.0_iwp
                s(16:20,2) =-0.538469310105683_iwp
                s(21:25,2) =-0.906179845938664_iwp
                wt(1) =0.056134348862429_iwp
                wt(2) =0.113400000000000_iwp
                wt(3) =0.134785072387521_iwp
                wt(4) =0.113400000000000_iwp
                wt(5) =0.056134348862429_iwp
                wt(6) =0.113400000000000_iwp
                wt(7) =0.229085404223991_iwp
                wt(8) =0.272286532550750_iwp
                wt(9) =0.229085404223991_iwp
                wt(10)=0.113400000000000_iwp
                wt(11)=0.134785072387521_iwp
                wt(12)=0.272286532550750_iwp
                wt(13)=0.323634567901235_iwp
                wt(14)=0.272286532550750_iwp
                wt(15)=0.134785072387521_iwp
                wt(16)=0.113400000000000_iwp
                wt(17)=0.229085404223991_iwp
                wt(18)=0.272286532550750_iwp
                wt(19)=0.229085404223991_iwp
                wt(20)=0.113400000000000_iwp
                wt(21)=0.056134348862429_iwp
                wt(22)=0.113400000000000_iwp
                wt(23)=0.134785072387521_iwp
                wt(24)=0.113400000000000_iwp
                wt(25)=0.056134348862429_iwp
            CASE DEFAULT
                WRITE(*,*)"wrong number of integrating points for a quadrilateral"
            END SELECT
            CASE('tetrahedron')
            !                       for tetrahedra weights multiplied by 1/6
            SELECT CASE(nip)
            CASE(1)
                s(1,1)=0.25_iwp
                s(1,2)=0.25_iwp
                s(1,3)=0.25_iwp
                wt(1)=1.0_iwp/6.0_iwp
            CASE(4)
                s(1,1)=0.58541020_iwp
                s(1,2)=0.13819660_iwp
                s(1,3)=s(1,2)
                s(2,2)=s(1,1)
                s(2,3)=s(1,2)
                s(2,1)=s(1,2)
                s(3,3)=s(1,1)
                s(3,1)=s(1,2)
                s(3,2)=s(1,2)
                s(4,1)=s(1,2)
                s(4,2)=s(1,2)
                s(4,3)=s(1,2)
                wt(1:4)=0.25_iwp/6.0_iwp
            CASE(5)
                s(1,1)=0.25_iwp
                s(1,2)=0.25_iwp
                s(1,3)=0.25_iwp
                s(2,1)=0.5_iwp
                s(2,2)=1.0_iwp/6.0_iwp
                s(2,3)=s(2,2)
                s(3,2)=0.5_iwp
                s(3,3)=1.0_iwp/6.0_iwp
                s(3,1)=s(3,3)
                s(4,3)=0.5_iwp
                s(4,1)=1.0_iwp/6.0_iwp
                s(4,2)=s(4,1)
                s(5,1)=1.0_iwp/6.0_iwp
                s(5,2)=s(5,1)
                s(5,3)=s(5,1)
                wt(1)=-0.8_iwp
                wt(2)=9.0_iwp/20.0_iwp
                wt(3:5)=wt(2)
                wt=wt/6.0_iwp
            CASE DEFAULT
                WRITE(*,*)"wrong number of integrating points for a tetrahedron"
            END SELECT
            CASE('hexahedron')
            SELECT CASE(nip)
            CASE(1)
                s(1,1:3)=0.0_iwp
                wt(1)=8.0_iwp
            CASE(8)
                s(1,1)= root3
                s(1,2)= root3
                s(1,3)= root3
                s(2,1)= root3
                s(2,2)= root3
                s(2,3)=-root3
                s(3,1)= root3
                s(3,2)=-root3
                s(3,3)= root3
                s(4,1)= root3
                s(4,2)=-root3
                s(4,3)=-root3
                s(5,1)=-root3
                s(5,2)= root3
                s(5,3)= root3
                s(6,1)=-root3
                s(6,2)=-root3
                s(6,3)= root3
                s(7,1)=-root3
                s(7,2)= root3
                s(7,3)=-root3
                s(8,1)=-root3
                s(8,2)=-root3
                s(8,3)=-root3
                wt=1.0_iwp
            CASE(14)
                b=0.795822426_iwp
                c=0.758786911_iwp
                wt(1:6)=0.886426593_iwp
                wt(7:14)=0.335180055_iwp
                s(1,1)=-b
                s(2,1)=b
                s(3,2)=-b
                s(4,2)=b
                s(5,3)=-b
                s(6,3)=b
                s(7:,:)=c
                s(7,1)=-c
                s(7,2)=-c
                s(7,3)=-c
                s(8,2)=-c
                s(8,3)=-c
                s(9,1)=-c
                s(9,3)=-c
                s(10,3)=-c
                s(11,1)=-c
                s(11,2)=-c
                s(12,2)=-c
                s(13,1)=-c
            CASE(15)
                b=1.0_iwp
                c       =0.674199862_iwp
                wt(1)   =1.564444444_iwp
                wt(2:7) =0.355555556_iwp
                wt(8:15)=0.537777778_iwp
                s(2,1)=-b
                s(3,1)=b
                s(4,2)=-b
                s(5,2)=b
                s(6,3)=-b
                s(7,3)=b
                s(8:,:)=c
                s(8,1)=-c
                s(8,2)=-c
                s(8,3)=-c
                s(9,2)=-c
                s(9,3)=-c
                s(10,1)=-c
                s(10,3)=-c
                s(11,3)=-c
                s(12,1)=-c
                s(12,2)=-c
                s(13,2)=-c
                s(14,1)=-c
            CASE(27)
                wt=(/5.0_iwp/9.0_iwp*v,8.0_iwp/9.0_iwp*v,5.0_iwp/9.0_iwp*v/)
                s(1:7:3,1)=-r15
                s(2:8:3,1)=0.0_iwp
                s(3:9:3,1)=r15
                s(1:3,3)=r15
                s(4:6,3)=0.0_iwp
                s(7:9,3)=-r15
                s(1:9,2)=-r15
                s(10:16:3,1)=-r15
                s(11:17:3,1)=0.0_iwp
                s(12:18:3,1)=r15
                s(10:12,3)=r15
                s(13:15,3)=0.0_iwp
                s(16:18,3)=-r15
                s(10:18,2)=0.0_iwp
                s(19:25:3,1)=-r15
                s(20:26:3,1)=0.0_iwp
                s(21:27:3,1)=r15
                s(19:21,3)= r15
                s(22:24,3)=0.0_iwp
                s(25:27,3)=-r15
                s(19:27,2)= r15
            CASE DEFAULT
                WRITE(*,*)"wrong number of integrating points for a hexahedron"
            END SELECT
            CASE DEFAULT
            WRITE(*,*)"not a valid element type"
            END SELECT
            RETURN
        END SUBROUTINE sample

        SUBROUTINE deemat_corrugated(dee,faceprop,coreprop,scheme,nsec)
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::faceprop(:,:),coreprop(:,:)
            REAL(iwp),INTENT(OUT)::dee(:,:)
            CHARACTER(LEN=30),INTENT(IN)::scheme
            INTEGER,INTENT(IN)::nsec
            INTEGER::i,j,k
            REAL,PARAMETER::zero=0.0_iwp,two_hundred=200.0_iwp,two=2.0_iwp,twelve=12.0_iwp,one=1.0_iwp,four=4.0_iwp,five=5.0_iwp,six=6.0_iwp,pi=3.14159265358979323846_iwp
            real(iwp)::matrix_q1(3,3), matrix_q2(2,2), Ez, matrix_a(3,3), matrix_b(3,3), matrix_c(3,3),matrix_d(2,2), matrix_e(2,2), matrix_f(2,2),&
             E3, v13=0.01_iwp, v23=0.01_iwp, E1, E2, v12, G12, G13, G23, Im, delta, thick, thick_c, t, period, height,K1, K2, K3, Kcorrect(14,14), tc, posai,&
             matrix_g(2,2),matrix_h(2,2)
            real(iwp),allocatable::mid_loc(:), xcor(:), c(:), s(:), Ex0(:), Ey0(:),vxy0(:), Gyz0(:), Gxz0(:), Gxy0(:),&
            vyx0(:), z(:), xhalf(:), zhalf(:), chalf(:)
            ALLOCATE(mid_loc(2*UBOUND(coreprop,1)+1), xcor(nsec), c(nsec), s(nsec), z(nsec), xhalf(nsec), zhalf(nsec), chalf(nsec), Ex0(nsec), Ey0(nsec),&
            vxy0(nsec), Gyz0(nsec), Gxz0(nsec), Gxy0(nsec), vyx0(nsec))

            SELECT case(scheme)
                case('full')
                    !get the total thickness
                    thick=zero!total thickness of a cardboard, include facesheet 
                    thick_c=zero!total thickness of corrugated core midsurface only
                    do i=1,UBOUND(faceprop,1)
                        thick=thick+faceprop(i,1)
                    end do
                    do i=1,UBOUND(coreprop,1)
                        thick=thick+coreprop(i,9)
                        thick_c=thick_c+(coreprop(i,9)-coreprop(i,1))
                    end do 
                    !get the location of midplane of each ply
                    mid_loc=zero
                    mid_loc(1)=faceprop(1,1)/two
                    do i=2,UBOUND(mid_loc,1)
                        if (MOD(i,2)/=0) then
                            mid_loc(i)=mid_loc(i-1)+faceprop((i+1)/2,1)/two+coreprop((i-1)/2,9)/two
                        else
                            mid_loc(i)=mid_loc(i-1)+coreprop(i/2,9)/two+faceprop(i/2,1)/two
                        end if
                    end do
                    mid_loc=mid_loc-thick/two

                case('simplified')
                    thick=zero
                    thick_c=zero
                    do i=1,UBOUND(coreprop,1)
                        thick=thick+coreprop(i,9)
                        thick_c=thick_c+coreprop(i,9)
                    end do 
                    mid_loc=zero
                    do i=2,UBOUND(mid_loc,1)
                        if (MOD(i,2)/=0) then
                            mid_loc(i)=coreprop((i-1)/2,9)*(i-1)/2
                        else
                            mid_loc(i)=mid_loc(i-1)+coreprop(i/2,9)/two
                        end if
                    end do
                    mid_loc=mid_loc-thick/two
                case DEFAULT
                    ERROR STOP 'calculation scheme for constitutive matrix is not correct' 
            end SELECT
            
            matrix_q1=zero
            matrix_q2=zero
            matrix_a=zero
            matrix_b=zero
            matrix_c=zero
            matrix_d=zero
            matrix_e=zero
            matrix_f=zero
            Ez=zero

            !contact author on email HenryLeung928@gmail.com for detailed code
            write(*,*) "contact author on email HenryLeung928@gmail.com for detailed code"

            !K1=sqrt(5.0_iwp/6.0_iwp)
            K1=sqrt(0.658_iwp)
            K2=K1
            K3=sqrt(1.0_iwp)
            !K3=one
            Kcorrect=zero
            do i=1,UBOUND(Kcorrect,1)
                SELECT case(i)
                case(4)
                    Kcorrect(i,i)=K3
                case(5)
                    Kcorrect(i,i)=K1
                case(6)
                    Kcorrect(i,i)=K1
                case(10)
                    Kcorrect(i,i)=K3
                case(11)
                    Kcorrect(i,i)=K2
                case(12)
                    Kcorrect(i,i)=K2
                case(13)
                    Kcorrect(i,i)=K2
                case(14)
                    Kcorrect(i,i)=K2
                case DEFAULT
                    Kcorrect(i,i)=one
                end select
            end do
            
            dee=zero
            dee(1:3,1:3)=matrix_a
            dee(1:3,7:9)=matrix_b
            dee(7:9,1:3)=matrix_b
            dee(7:9,7:9)=matrix_c

            dee(5:6,5:6)=matrix_d
            dee(5:6,11:12)=matrix_e
            dee(11:12,5:6)=matrix_e
            dee(11:12,11:12)=matrix_f

            dee(5:6,13:14)=matrix_f
            dee(13:14,5:6)=matrix_f

            dee(11:12,13:14)=matrix_g
            dee(13:14,11:12)=matrix_g

            dee(13:14,13:14)=matrix_h

            dee(4,4)=one/Ez*thick
            dee(4,1)=v23*dee(4,4)
            dee(4,2)=0.96_iwp*dee(4,4)
            dee(1,4)=dee(4,1)
            dee(2,4)=dee(4,2)

            dee(10,10)=dee(4,4)*thick_c**3/twelve
            dee(10,1)=v23*dee(10,10)
            dee(10,2)=0.96_iwp*dee(10,10)
            dee(1,10)=dee(10,1)
            dee(2,10)=dee(10,2)
            dee=MATMUL(MATMUL(TRANSPOSE(Kcorrect),dee),Kcorrect)

            RETURN
        END SUBROUTINE deemat_corrugated

        SUBROUTINE deemat_isotropic(dee,e,v,t)
            !isotropic test
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::e,v,t
            REAL(iwp),INTENT(OUT)::dee(:,:)
            REAL(iwp)::zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp,twelve=12.0_iwp
            real(iwp)::matrix_q(6,6), matrix_a(6,6), matrix_b(6,6), matrix_c(6,6), &
            matrix_d(6,6), matrix_e(6,6), matrix_f(6,6), matrix_g(6,6), matrix_h(6,6), matrix_i(6,6),K1, K2, K3, Kcorrect(20,20)
            
            matrix_q=zero
            matrix_q(1,1)=e/(one+v)/(1-two*v)*(1-v)
            matrix_q(2,2)=e/(one+v)/(1-two*v)*(1-v)
            matrix_q(1,2)=e/(one+v)/(1-two*v)*v
            matrix_q(2,1)=matrix_q(1,2)
            matrix_q(3,3)=e/(one+v)/two
            matrix_q(4,4)=e/(one+v)/(1-two*v)*(1-v)
            matrix_q(1,4)=e/(one+v)/(1-two*v)*v
            matrix_q(2,4)=e/(one+v)/(1-two*v)*v
            matrix_q(4,1)=matrix_q(1,4)
            matrix_q(4,2)=matrix_q(2,4)
            matrix_q(5,5)=e/(one+v)/two
            matrix_q(6,6)=e/(one+v)/two

            matrix_a=matrix_q*t
            matrix_b=zero 
            matrix_c=matrix_q*t**3/twelve
            matrix_d=zero
            matrix_e=matrix_q*t**5/80.0_iwp

            dee=zero
            dee(1:6,1:6)=matrix_a
            dee(1:6,7:12)=matrix_b
            dee(1:6,13:14)=matrix_c(1:6,5:6)


            dee(7:12,1:6)=matrix_b
            dee(7:12,7:12)=matrix_c
            dee(7:12,13:14)=matrix_d(1:6,5:6)


            dee(13:14,1:6)=matrix_c(5:6,1:6)
            dee(13:14,7:12)=matrix_d(5:6,1:6)
            dee(13:14,13:14)=matrix_e(5:6,5:6)


            !dee=MATMUL(MATMUL(TRANSPOSE(Kcorrect),dee),Kcorrect)
            RETURN
        END SUBROUTINE deemat_isotropic  

        SUBROUTINE shape_der_m(der,points,i)
            !
            !   This subroutine produces derivatives of shape functions withe respect
            !   to local coordinates.
            !
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            INTEGER,INTENT(IN)::i
            REAL(iwp),INTENT(IN)::points(:,:)
            REAL(iwp),INTENT(OUT)::der(:,:)
            REAL(iwp)::eta,xi,zeta,xi0,eta0,zeta0,etam,etap,xim,xip,c1,c2,c3 
            REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9,x2p1,x2m1,e2p1,e2m1,zetam,zetap
            REAL,PARAMETER::zero=0.0_iwp,pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,  &
            pt75=0.75_iwp,one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp,&
            d6=6.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d10=10.0_iwp,d11=11.0_iwp,            &
            d12=12.0_iwp,d16=16.0_iwp,d18=18.0_iwp,d27=27.0_iwp,d32=32.0_iwp,      &
            d36=36.0_iwp,d54=54.0_iwp,d64=64.0_iwp,d128=128.0_iwp
            INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
            ndim=UBOUND(der,1)
            nod= UBOUND(der,2)
            SELECT CASE(ndim)
            CASE(1)   ! one dimensional elements
            xi=points(i,1)
            SELECT CASE(nod)
            CASE(2)
                der(1,1)=-pt5 
                der(1,2)= pt5
            CASE(3)
                t1=-one-xi 
                t2=-xi  
                t3=one-xi
                der(1,1)=-(t3+t2)/two  
                der(1,2)=(t3+t1)    
                der(1,3)=-(t2+t1)/two   
            CASE(4)
                t1=-one-xi 
                t2=-one/d3-xi 
                t3=one/d3-xi 
                t4=one-xi
                der(1,1)=-(t3*t4+t2*t4+t2*t3)*d9/d16     
                der(1,2)=(t3*t4+t1*t4+t1*t3)*d27/d16 
                der(1,3)=-(t2*t4+t1*t4+t1*t2)*d27/d16 
                der(1,4)=(t2*t3+t1*t3+t1*t2)*d9/d16   
            CASE(5)
                t1=-one-xi 
                t2=-pt5-xi 
                t3=-xi 
                t4=pt5-xi 
                t5=one-xi
                der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*two/d3   
                der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*d8/d3
                der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*d4 
                der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*d8/d3
                der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*two/d3
            CASE DEFAULT
                WRITE(*,*)"wrong number of nodes in shape_der"        
            END SELECT
            CASE(2)      ! two dimensional elements
            xi=points(i,1)
            eta=points(i,2) 
            c1=xi 
            c2=eta 
            c3=one-c1-c2
            etam=pt25*(one-eta)
            etap=pt25*(one+eta)
            xim= pt25*(one-xi)
            xip= pt25*(one+xi)
            x2p1=two*xi+one 
            x2m1=two*xi-one 
            e2p1=two*eta+one 
            e2m1=two*eta-one
            SELECT CASE(nod)
            CASE(3)
                der(1,1)=one
                der(1,3)=zero
                der(1,2)=-one
                der(2,1)=zero
                der(2,3)=one
                der(2,2)=-one
            CASE(6) 
                der(1,1)=d4*c1-one 
                der(1,6)=d4*c2
                der(1,5)=zero  
                der(1,4)=-d4*c2
                der(1,3)=-(d4*c3-one)
                der(1,2)=d4*(c3-c1)
                der(2,1)=zero
                der(2,6)=d4*c1 
                der(2,5)=d4*c2-one
                der(2,4)=d4*(c3-c2)
                der(2,3)=-(d4*c3-one)  
                der(2,2)=-d4*c1
            CASE(10)                          
                der(1,1)=(d27*c1**2-d18*c1+two)/two
                der(1,9)=(d9*(d6*c1-one)*c2)/two
                der(1,8)=(d9*(d3*c2-one)*c2)/two
                der(1,7)=zero
                der(1,6)=-(d9*(d3*c2-one)*c2)/two
                der(1,5)= (d9*(d6*c1+d6*c2-d5)*c2)/two
                der(1,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
                der(1,3)= (d9*(d9*c1**2+d12*c1*c2-d10*c1+d3*c2**2-d5*c2+two))/two
                der(1,2)=-(d9*(d9*c1**2+d6*c1*c2-d8*c1-c2+one))/two
                der(1,10)=-d27*(((c2-one)+c1)+c1)*c2
                der(2,1)=zero
                der(2,9)= (d9*(d3*c1-one)*c1)/two
                der(2,8)= (d9*(d6*c2-one)*c1)/two
                der(2,7)=(d27*c2**2-d18*c2+two)/two
                der(2,6)=-(d9*((c1+c2-one)*(d6*c2-one)+(d3*c2-one)*c2))/two
                der(2,5)= (d9*(d3*c1**2+d12*c1*c2-d5*c1+d9*c2**2-d10*c2+two))/two
                der(2,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
                der(2,3)= (d9*(d6*c1+d6*c2-d5)*c1)/two
                der(2,2)=-(d9*(d3*c1-one)*c1)/two
                der(2,10)=-d27*(((c2-one)+c1)+c2)*c1
            CASE(15)                          
                t1=c1-pt25  
                t2=c1-pt5 
                t3=c1-pt75   
                t4=c2-pt25
                t5=c2-pt5   
                t6=c2-pt75 
                t7=c3-pt25  
                t8=c3-pt5 
                t9=c3-pt75
                der(1,1)=d32/d3*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
                der(1,12)=d128/d3*c2*(t2*(t1+c1)+c1*t1) 
                der(1,11)=d64*c2*t4*(t1+c1)
                der(1,10)=d128/d3*c2*t4*t5  
                der(1,9)=zero 
                der(1,8)=-d128/d3*c2*t4*t5
                der(1,7)=-d64*c2*t4*(t7+c3) 
                der(1,6)=-d128/d3*c2*(t8*(t7+c3)+c3*t7)
                der(1,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
                der(1,4)=d128/d3*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
                der(1,3)=d64*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
                der(1,2)=d128/d3*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
                der(1,13)=d128*c2*(c3*(t1+c1)-c1*t1) 
                der(1,15)=d128*c2*t4*(c3-c1)
                der(1,14)=d128*c2*(c3*t7-c1*(t7+c3))
                der(2,1)=zero 
                der(2,12)=d128/d3*c1*t1*t2
                der(2,11)=d64*c1*t1*(t4+c2)
                der(2,10)=d128/d3*c1*(t5*(t4+c2)+c2*t4)
                der(2,9)=d32/d3*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
                der(2,8)=d128/d3*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
                der(2,7)=d64*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
                der(2,6)=d128/d3*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
                der(2,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
                der(2,4)=-d128/d3*c1*(t8*(t7+c3)+c3*t7)
                der(2,3)=-d64*c1*t1*(t7+c3)  
                der(2,2)=-d128/d3*c1*t1*t2
                der(2,13)=d128*c1*t1*(c3-c2)
                der(2,15)=d128*c1*(c3*(t4+c2)-c2*t4)
                der(2,14)=d128*c1*(c3*t7-c2*(c3+t7))        
            CASE (4) 
                !the order is changed                                                             
                der(1,1)=-etam
                der(1,4)=-etap
                der(1,3)=etap
                der(1,2)=etam
                der(2,1)=-xim
                der(2,4)=xim
                der(2,3)=xip
                der(2,2)=-xip
            CASE(5)
                der(1,1)=-etam+pt5*xi*(one-eta**2)
                der(1,2)=-etap+pt5*xi*(one-eta**2)
                der(1,3)=etap+pt5*xi*(one-eta**2)
                der(1,4)=etam+pt5*xi*(one-eta**2)
                der(1,5)=-two*xi*(one-eta**2)
                der(2,1)=-xim+pt5*eta*(one-xi**2)
                der(2,2)=xim+pt5*eta*(one-xi**2)
                der(2,3)=xip+pt5*eta*(one-xi**2)
                der(2,4)=-xip+pt5*eta*(one-xi**2)
                der(2,5)=-two*eta*(one-xi**2)
            CASE(8) !order changed
                der(1,1)=etam*(two*xi+eta)
                der(1,8)=-d8*etam*etap
                der(1,4)=etap*(two*xi-eta)
                der(1,7)=-d4*etap*xi
                der(1,3)=etap*(two*xi+eta)
                der(1,6)=d8*etap*etam
                der(1,2)=etam*(two*xi-eta)
                der(1,5)=-d4*etam*xi

                der(2,1)=xim*(xi+two*eta)
                der(2,8)=-d4*xim*eta
                der(2,4)=xim*(two*eta-xi)
                der(2,7)=d8*xim*xip
                der(2,3)=xip*(xi+two*eta)
                der(2,6)=-d4*xip*eta
                der(2,2)=xip*(two*eta-xi)
                der(2,5)=-d8*xim*xip   
            CASE(9)
                etam=eta-one
                etap=eta+one
                xim=xi-one
                xip=xi+one
                der(1,1)=pt25*x2m1*eta*etam  
                der(1,2)=-pt5*x2m1*etap*etam
                der(1,3)=pt25*x2m1*eta*etap  
                der(1,4)=-xi*eta*etap
                der(1,5)=pt25*x2p1*eta*etap  
                der(1,6)=-pt5*x2p1*etap*etam
                der(1,7)=pt25*x2p1*eta*etam  
                der(1,8)=-xi*eta*etam
                der(1,9)=two*xi*etap*etam    
                der(2,1)=pt25*xi*xim*e2m1
                der(2,2)=-xi*xim*eta        
                der(2,3)=pt25*xi*xim*e2p1
                der(2,4)=-pt5*xip*xim*e2p1   
                der(2,5)=pt25*xi*xip*e2p1
                der(2,6)=-xi*xip*eta        
                der(2,7)=pt25*xi*xip*e2m1
                der(2,8)=-pt5*xip*xim*e2m1   
                der(2,9)=two*xip*xim*eta
            CASE DEFAULT
                WRITE(*,*)"wrong number of nodes in shape_der"        
            END SELECT
            CASE(3)  ! d3 dimensional elements
            xi=points(i,1)
            eta=points(i,2)
            zeta=points(i,3)
            etam=one-eta 
            xim=one-xi
            zetam=one-zeta
            etap=eta+one 
            xip=xi+one 
            zetap=zeta+one
            SELECT CASE(nod)
            CASE(4)
                der(1:3,1:4)=zero
                der(1,1)=one
                der(2,2)=one  
                der(3,3)=one
                der(1,4)=-one 
                der(2,4)=-one 
                der(3,4)=-one  
            CASE(8)
                der(1,1)=-pt125*etam*zetam    
                der(1,2)=-pt125*etam*zetap
                der(1,3)= pt125*etam*zetap     
                der(1,4)= pt125*etam*zetam
                der(1,5)=-pt125*etap*zetam    
                der(1,6)=-pt125*etap*zetap
                der(1,7)= pt125*etap*zetap     
                der(1,8)= pt125*etap*zetam
                der(2,1)=-pt125*xim*zetam     
                der(2,2)=-pt125*xim*zetap
                der(2,3)=-pt125*xip*zetap     
                der(2,4)=-pt125*xip*zetam
                der(2,5)= pt125*xim*zetam      
                der(2,6)= pt125*xim*zetap
                der(2,7)= pt125*xip*zetap      
                der(2,8)= pt125*xip*zetam
                der(3,1)=-pt125*xim*etam      
                der(3,2)= pt125*xim*etam
                der(3,3)= pt125*xip*etam       
                der(3,4)=-pt125*xip*etam
                der(3,5)=-pt125*xim*etap      
                der(3,6)= pt125*xim*etap
                der(3,7)= pt125*xip*etap       
                der(3,8)=-pt125*xip*etap  
            CASE(14) ! type 6 element
                der(1,1)= (two*xi*eta+two*xi*zeta+d4*xi+eta*zeta+eta+zeta)*          &
                (eta-one)*(zeta-one)/d8
                der(1,2)=-(two*xi*eta-two*xi*zeta+d4*xi-eta*zeta+eta-zeta)*          &
                (eta-one)*(zeta+one)/d8
                der(1,3)=-(two*xi*eta-two*xi*zeta+d4*xi+eta*zeta-eta+zeta)*          &
                (eta-one)*(zeta+one)/d8
                der(1,4)= (two*xi*eta+two*xi*zeta+d4*xi-eta*zeta-eta-zeta)*          &
                (eta-one)*(zeta-one)/d8
                der(1,5)= -(eta-one)*(zeta+one)*(zeta-one)*xi 
                der(1,6)=-(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
                der(1,7)=  (eta+one)*(eta-one)*(zeta+one)*xi
                der(1,8)= (eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
                der(1,9)= -(eta+one)*(eta-one)*(zeta-one)*xi  
                der(1,10)= (two*xi*eta-two*xi*zeta-d4*xi+eta*zeta+eta-zeta)*         &
                (eta+one)*(zeta-one)/d8
                der(1,11)=-(two*xi*eta+two*xi*zeta-d4*xi-eta*zeta+eta+zeta)*         &
                (eta+one)*(zeta+one)/d8
                der(1,12)=-(two*xi*eta+two*xi*zeta-d4*xi+eta*zeta-eta-zeta)*         &
                (eta+one)*(zeta+one)/d8
                der(1,13)= (two*xi*eta-two*xi*zeta-d4*xi-eta*zeta-eta+zeta)*         &
                (eta+one)*(zeta-one)/d8
                der(1,14)=  (eta+one)*(zeta+one)*(zeta-one)*xi
                der(2,1)= (two*xi*eta+xi*zeta+xi+two*eta*zeta+d4*eta+zeta)*          &
                (xi-one)*(zeta-one)/d8                                  
                der(2,2)=-(two*xi*eta-xi*zeta+xi-two*eta*zeta+d4*eta-zeta)*          &
                (xi-one)*(zeta+one)/d8
                der(2,3)=-(two*xi*eta-xi*zeta+xi+two*eta*zeta-d4*eta+zeta)*          &
                (xi+one)*(zeta+one)/d8
                der(2,4)= (two*xi*eta+xi*zeta+xi-two*eta*zeta-d4*eta-zeta)*          &
                (xi+one)*(zeta-one)/d8
                der(2,5)=-(xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
                der(2,6)= -(xi-one)*(zeta+one)*(zeta-one)*eta
                der(2,7)=  (xi+one)*(xi-one)*(zeta+one)*eta
                der(2,8)=  (xi+one)*(zeta+one)*(zeta-one)*eta
                der(2,9)= -(xi+one)*(xi-one)*(zeta-one)*eta
                der(2,10)= (two*xi*eta-xi*zeta-xi+two*eta*zeta+d4*eta-zeta)*         &
                (xi-one)*(zeta-one)/d8
                der(2,11)=-(two*xi*eta+xi*zeta-xi-two*eta*zeta+d4*eta+zeta)*         &
                (xi-one)*(zeta+one)/d8
                der(2,12)=-(two*xi*eta+xi*zeta-xi+two*eta*zeta-d4*eta-zeta)*         &
                (xi+one)*(zeta+one)/d8
                der(2,13)= (two*xi*eta-xi*zeta-xi-two*eta*zeta-d4*eta+zeta)*         &
                (xi+one)*(zeta-one)/d8
                der(2,14)= (xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
                der(3,1)= (xi*eta+two*xi*zeta+xi+two*eta*zeta+eta+d4*zeta)*          &
                (xi-one)*(eta-one)/d8
                der(3,2)=-(xi*eta-two*xi*zeta+xi-two*eta*zeta+eta-d4*zeta)*          &
                (xi-one)*(eta-one)/d8
                der(3,3)=-(xi*eta-two*xi*zeta+xi+two*eta*zeta-eta+d4*zeta)*          &
                (xi+one)*(eta-one)/d8
                der(3,4)= (xi*eta+two*xi*zeta+xi-two*eta*zeta-eta-d4*zeta)*          &
                (xi+one)*(eta-one)/d8
                der(3,5)= -(xi+one)*(xi-one)*(eta-one)*zeta  
                der(3,6)= -(xi-one)*(eta+one)*(eta-one)*zeta  
                der(3,7)= (xi+one)*(xi-one)*(eta+one)*(eta-one)/two
                der(3,8)=  (xi+one)*(eta+one)*(eta-one)*zeta
                der(3,9)=-(xi+one)*(xi-one)*(eta+one)*(eta-one)/two
                der(3,10)= (xi*eta-two*xi*zeta-xi+two*eta*zeta+eta-d4*zeta)*         &
                (xi-one)*(eta+one)/d8
                der(3,11)=-(xi*eta+two*xi*zeta-xi-two*eta*zeta+eta+d4*zeta)*         &
                (xi-one)*(eta+one)/d8
                der(3,12)=-(xi*eta+two*xi*zeta-xi+two*eta*zeta-eta-d4*zeta)*         &
                (xi+one)*(eta+one)/d8
                der(3,13)= (xi*eta-two*xi*zeta-xi-two*eta*zeta-eta+d4*zeta)*         &
                (xi+one)*(eta+one)/d8
                der(3,14)=  (xi+one)*(xi-one)*(eta+one)*zeta
            CASE(20)
                xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
                etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
                zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
                DO l=1,20
                xi0=xi*xii(l)
                eta0=eta*etai(l)
                zeta0=zeta*zetai(l)
                IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
                    der(1,l)=-pt5*xi*(one+eta0)*(one+zeta0)
                    der(2,l)=pt25*etai(l)*(one-xi*xi)*(one+zeta0)
                    der(3,l)=pt25*zetai(l)*(one-xi*xi)*(one+eta0)
                ELSE IF(l>=9.AND.l<=12)THEN
                    der(1,l)=pt25*xii(l)*(one-eta*eta)*(one+zeta0)
                    der(2,l)=-pt5*eta*(one+xi0)*(one+zeta0)
                    der(3,l)=pt25*zetai(l)*(one+xi0)*(one-eta*eta)
                ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
                    der(1,l)=pt25*xii(l)*(one+eta0)*(one-zeta*zeta)
                    der(2,l)=pt25*etai(l)*(one+xi0)*(one-zeta*zeta)
                    der(3,l)=-pt5*zeta*(one+xi0)*(one+eta0)
                ELSE
                    der(1,l)=pt125*xii(l)*(one+eta0)*(one+zeta0)*                    &
                    (two*xi0+eta0+zeta0-one)
                    der(2,l)=pt125*etai(l)*(one+xi0)*(one+zeta0)*                    &
                    (xi0+two*eta0+zeta0-one)
                    der(3,l)=pt125*zetai(l)*(one+xi0)*(one+eta0)*                    &
                    (xi0+eta0+two*zeta0-one)
                END IF
                END DO 
            CASE DEFAULT
                WRITE(*,*)"wrong number of nodes in shape_der"        
            END SELECT
            CASE DEFAULT
            WRITE(*,*)"wrong number of dimensions in shape_der"
            END SELECT
            RETURN
        END SUBROUTINE shape_der_m                

        SUBROUTINE invert(matrix)
            !
            ! This subroutine inverts a small square matrix onto itself.
            !
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN OUT)::matrix(:,:)
            REAL(iwp)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
            INTEGER::ndim,i,k
            ndim=UBOUND(matrix,1)
            IF(ndim==2)THEN
            det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
            j11=matrix(1,1)
            matrix(1,1)=matrix(2,2)
            matrix(2,2)=j11
            matrix(1,2)=-matrix(1,2)
            matrix(2,1)=-matrix(2,1)
            matrix=matrix/det
            ELSE IF(ndim==3)THEN
            det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
            det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
            det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
            j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
            j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
            j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
            j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
            j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
            j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
            j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
            j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
            j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
            matrix(1,1)=j11
            matrix(1,2)=j12
            matrix(1,3)=j13
            matrix(2,1)=j21
            matrix(2,2)=j22
            matrix(2,3)=j23
            matrix(3,1)=j31
            matrix(3,2)=j32
            matrix(3,3)=j33
            matrix=matrix/det
            ELSE
            DO k=1,ndim
                con=matrix(k,k)
                matrix(k,k)=1.0_iwp
                matrix(k,:)=matrix(k,:)/con
                DO i=1,ndim
                IF(i/=k)THEN
                    con=matrix(i,k)
                    matrix(i,k)=0.0_iwp
                    matrix(i,:)=matrix(i,:)-matrix(k,:)*con
                END IF
                END DO
            END DO
            END IF
            RETURN
        END SUBROUTINE invert

        SUBROUTINE beemat_wHSDT(bee,deriv,fun)
            !
            ! This subroutine forms the bee matrix for wHSDT
            !
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::deriv(:,:)!deriv(ndim,nnode)
            REAL(iwp),INTENT(IN)::fun(:)
            REAL(iwp),INTENT(OUT)::bee(:,:)!bee(nst,ndofel)
            INTEGER::j,nod
            REAL,PARAMETER::zero=0.0_iwp,one=1.0_iwp,two=2.0_iwp,three=3.0_iwp,four=4.0_iwp
            bee=zero
            nod=UBOUND(deriv,2)
            
            do j=1,nod
                bee(1,7*j-6)=deriv(1,j)
                bee(2,7*j-5)=deriv(2,j)
                bee(3,7*j-6)=deriv(2,j)
                bee(3,7*j-5)=deriv(1,j)
                bee(4,7*j-1)=fun(j)
                bee(5,7*j-4)=deriv(1,j)
                bee(5,7*j-2)=fun(j)
                bee(6,7*j-4)=deriv(2,j)
                bee(6,7*j-3)=-fun(j)
                bee(7,7*j-2)=deriv(1,j)
                bee(8,7*j-3)=-deriv(2,j)
                bee(9,7*j-2)=deriv(2,j)
                bee(9,7*j-3)=-deriv(1,j)
                bee(10,7*j)=two*fun(j)
                bee(11,7*j-1)=deriv(1,j)
                bee(12,7*j-1)=deriv(2,j)
                bee(13,7*j)=deriv(1,j)
                bee(14,7*j)=deriv(2,j)
            end do
            RETURN
        END SUBROUTINE beemat_wHSDT

        FUNCTION determinant(jac)RESULT(det)
            !
            ! This function returns the determinant of a 1x1, 2x2 or 3x3
            ! Jacobian matrix.
            !
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::jac(:,:)
            REAL(iwp)::det
            INTEGER::it 
            it=UBOUND(jac,1)  
            SELECT CASE(it)
            CASE(1)
            det=1.0_iwp
            CASE(2)
            det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
            CASE(3)
            det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
            det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
            det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
            CASE DEFAULT
            WRITE(*,*)' wrong dimension for Jacobian matrix'
            END SELECT
            RETURN
        END FUNCTION determinant

        SUBROUTINE shape_fun_m(fun,points,i)
            !
            !   This subroutine computes the values of the shape functions.
            !   to local coordinates
            !
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            INTEGER,INTENT(in)::i
            REAL(iwp),INTENT(IN)::points(:,:)
            REAL(iwp),INTENT(OUT)::fun(:)
            REAL(iwp)::eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3     
            REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9
            REAL(iwp)::zeta,xi0,eta0,zeta0
            INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
            REAL,PARAMETER::pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,pt75=0.75_iwp, &
            one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d8=8.0_iwp,d9=9.0_iwp,   &
            d16=16.0_iwp,d27=27.0_iwp,d32=32.0_iwp,d64=64.0_iwp,d128=128.0_iwp
            ndim=UBOUND(points,2)
            nod=UBOUND(fun,1)  
            SELECT CASE(ndim)
            CASE(1) ! one dimensional case
            xi=points(i,1)
            SELECT CASE(nod)
            CASE(2)
                t1=-one-xi 
                t2= one-xi
                fun(1)=t2/two 
                fun(2)=-t1/two
            CASE(3)
                t1=-one-xi 
                t2=-xi 
                t3=one-xi
                fun(1)=t2*t3/two 
                fun(2)=-t1*t3 
                fun(3)=t1*t2/two
            CASE(4)
                t1=-one-xi 
                t2=-one/d3-xi 
                t3=one/d3-xi 
                t4=one-xi
                fun(1)=t2*t3*t4*d9/d16  
                fun(2)=-t1*t3*t4*d27/d16
                fun(3)=t1*t2*t4*d27/d16 
                fun(4)=-t1*t2*t3*d9/d16
            CASE(5)
                t1=-one -xi 
                t2=-pt5-xi 
                t3=-xi 
                t4=pt5-xi 
                t5=one-xi
                fun(1)=t2*t3*t4*t5*two/d3 
                fun(2)=-t1*t3*t4*t5*d8/d3
                fun(3)=t1*t2*t4*t5*d4 
                fun(4)=-t1*t2*t3*t5*d8/d3
                fun(5)=t1*t2*t3*t4*two/d3
            CASE DEFAULT
                WRITE(*,*)"wrong number of nodes in shape_fun"
            END SELECT
            CASE(2) ! two dimensional case
            c1=points(i,1)
            c2=points(i,2)
            c3=one-c1-c2 
            xi=points(i,1)
            eta=points(i,2)
            etam=pt25*(one-eta)
            etap=pt25*(one+eta)
            xim=pt25*(one-xi)
            xip=pt25*(one+xi)
            SELECT CASE(nod)
            CASE(3)
                fun = (/c1,c3,c2/)  
            CASE(6)
                fun(1)=(two*c1-one)*c1 
                fun(2)=d4*c3*c1
                fun(3)=(two*c3-one)*c3 
                fun(4)=d4*c2*c3      
                fun(5)=(two*c2-one)*c2
                fun(6)=d4*c1*c2 
            CASE(10)
                fun(1)= ((d3*c1-one)*(d3*c1-two)*c1)/two
                fun(2)= -(d9*(d3*c1-one)*(c1+c2-one)*c1)/two
                fun(3)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c1)/two
                fun(4)=-((d3*c1+d3*c2-one)*(d3*c1+d3*c2-two)*(c1+c2-one))/two    
                fun(5)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c2)/two
                fun(6)= -(d9*(c1+c2-one)*(d3*c2-one)*c2)/two
                fun(7)= ((d3*c2-one)*(d3*c2-two)*c2)/two
                fun(8)=  (d9*(d3*c2-one)*c1*c2)/two
                fun(9)=  (d9*(d3*c1-one)*c1*c2)/two
                fun(10)=-d27*((c2-one)+c1)*c1*c2
            CASE(15)
                t1=c1-pt25  
                t2=c1-pt5 
                t3=c1-pt75   
                t4=c2-pt25
                t5=c2-pt5   
                t6=c2-pt75 
                t7=c3-pt25  
                t8=c3-pt5 
                t9=c3-pt75
                fun(1)=d32/d3*c1*t1*t2*t3   
                fun(2)=d128/d3*c3*c1*t1*t2
                fun(3)=d64*c3*c1*t1*t7      
                fun(4)=d128/d3*c3*c1*t7*t8
                fun(5)=d32/d3*c3*t7*t8*t9   
                fun(6)=d128/d3*c2*c3*t7*t8
                fun(7)=d64*c2*c3*t4*t7      
                fun(8)=d128/d3*c2*c3*t4*t5
                fun(9)=d32/d3*c2*t4*t5*t6   
                fun(10)=d128/d3*c1*c2*t4*t5
                fun(11)=d64*c1*c2*t1*t4     
                fun(12)=d128/d3*c1*c2*t1*t2
                fun(13)=d128*c1*c2*t1*c3    
                fun(15)=d128*c1*c2*c3*t4
                fun(14)=d128*c1*c2*c3*t7      
            CASE(4)
                !order changed
                fun(1)=d4*xim*etam
                fun(2)=d4*xip*etam
                fun(3)=d4*xip*etap
                fun(4)=d4*xim*etap
            CASE(5)
                fun=(/d4*xim*etam-pt25*(one-xi**2)*(one-eta**2),			  &
                    d4*xim*etap-pt25*(one-xi**2)*(one-eta**2),			  &
                    d4*xip*etap-pt25*(one-xi**2)*(one-eta**2),			  &
                    d4*xip*etam-pt25*(one-xi**2)*(one-eta**2),			  &
                    (one-xi**2)*(one-eta**2)/)
            CASE(8)!order is changed
                fun=(/d4*etam*xim*(-xi-eta-one),&
                    d4*xip*etam*(xi-eta-one),&
                    d4*etap*xip*(xi+eta-one),&
                    d4*etap*xim*(-xi+eta-one),&
                    d32*xim*xip*etam,&
                    d32*etap*xip*etam,&
                    d32*xim*xip*etap,&
                    d32*etam*xim*etap /)
            CASE(9)
                etam=eta-one
                etap=eta+one
                xim=xi-one
                xip=xi+one
                fun=(/pt25*xi*xim*eta*etam,-pt5*xi*xim*etap*etam,                    &
                    pt25*xi*xim*eta*etap,-pt5*xip*xim*eta*etap,                    &
                    pt25*xi*xip*eta*etap,-pt5*xi*xip*etap*etam,                    &
                    pt25*xi*xip*eta*etam,-pt5*xip*xim*eta*etam,                    &
                    xip*xim*etap*etam/)
            CASE DEFAULT
                WRITE(*,*)"wrong number of nodes in shape_fun"
            END SELECT
            CASE(3) ! d3 dimensional case
            xi=points(i,1)
            eta=points(i,2)
            zeta=points(i,3)
            etam=one-eta 
            xim=one-xi  
            zetam=one-zeta
            etap=eta+one 
            xip=xi+one   
            zetap=zeta+one
            SELECT CASE(nod)
            CASE(4)
                fun(1)=xi   
                fun(2)=eta 
                fun(3)=zeta 
                fun(4)=one-fun(1)-fun(2)-fun(3)
            CASE(8)
                fun=(/pt125*xim*etam*zetam,pt125*xim*etam*zetap,                     &
                    pt125*xip*etam*zetap,pt125*xip*etam*zetam,                     &
                    pt125*xim*etap*zetam,pt125*xim*etap*zetap,                     &
                    pt125*xip*etap*zetap,pt125*xip*etap*zetam/)
            CASE(14) ! type 6 element
                fun(1) = (xi*eta+xi*zeta+two*xi+eta*zeta+two*eta+two*zeta+two)*      &
                (xi-one)*(eta-one)*(zeta-one)/d8
                fun(2) =-(xi*eta-xi*zeta+two*xi-eta*zeta+two*eta-two*zeta+two)*      &
                (xi-one)*(eta-one)*(zeta+one)/d8
                fun(3) =-(xi*eta-xi*zeta+two*xi+eta*zeta-two*eta+two*zeta-two)*      &
                (xi+one)*(eta-one)*(zeta+one)/d8
                fun(4) = (xi*eta+xi*zeta+two*xi-eta*zeta-two*eta-two*zeta-two)*      &
                (xi+one)*(eta-one)*(zeta-one)/d8
                fun(5) =-(xi+one)*(xi-one)*(eta-one)*(zeta+one)*(zeta-one)/two
                fun(6) =-(xi-one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
                fun(7) = (xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta+one)/two
                fun(8) = (xi+one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
                fun(9) =-(xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta-one)/two
                fun(10)= (xi*eta-xi*zeta-two*xi+eta*zeta+two*eta-two*zeta-two)*      &
                (xi-one)*(eta+one)*(zeta-one)/d8
                fun(11)=-(xi*eta+xi*zeta-two*xi-eta*zeta+two*eta+two*zeta-two)*      &
                (xi-one)*(eta+one)*(zeta+one)/d8
                fun(12)=-(xi*eta+xi*zeta-two*xi+eta*zeta-two*eta-two*zeta+two)*      &
                (xi+one)*(eta+one)*(zeta+one)/d8
                fun(13)= (xi*eta-xi*zeta-two*xi-eta*zeta-two*eta+two*zeta+two)*      &
                (xi+one)*(eta+one)*(zeta-one)/d8
                fun(14)= (xi+one)*(xi-one)*(eta+one)*(zeta+one)*(zeta-one)/two
            CASE(20)
                xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
                etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
                zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
                DO l=1,20
                xi0=xi*xii(l)
                eta0=eta*etai(l)
                zeta0=zeta*zetai(l)
                IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
                    fun(l)=pt25*(one-xi*xi)*(one+eta0)*(one+zeta0)
                ELSE IF(l>=9.AND.l<=12)THEN
                    fun(l)=pt25*(one+xi0)*(one-eta*eta)*(one+zeta0)
                ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18)THEN
                    fun(l)=pt25*(one+xi0)*(one+eta0)*(one-zeta*zeta)
                ELSE
                    fun(l)=pt125*(one+xi0)*(one+eta0)*(one+zeta0)*(xi0+eta0+zeta0-2)
                END IF
                END DO
            CASE DEFAULT
                WRITE(*,*)"wrong number of nodes in shape_fun"
            END SELECT
            CASE DEFAULT
            WRITE(*,*)"wrong number of dimensions in shape_fun"
            END SELECT
            RETURN
        END SUBROUTINE shape_fun_m

        FUNCTION norm(x)RESULT(l2n)
            !
            ! THis function calculates the l2 norm of vector x
            !
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::x(:)
            REAL(iwp)::l2n
            l2n=SQRT(SUM(x**2))
            RETURN
        END FUNCTION norm

        FUNCTION cross_product3(a, b)RESULT(c)
            !aXb=c
            IMPLICIT NONE
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp):: a(3), b(3)
            REAL(iwp):: c(3)
            
            c(1) = a(2) * b(3) - a(3) * b(2)
            c(2) = a(3) * b(1) - a(1) * b(3)
            c(3) = a(1) * b(2) - a(2) * b(1)
            
            RETURN
        END FUNCTION cross_product3 

        FUNCTION dot_product(vec1, vec2) result(dp)
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            real(iwp):: vec1(:), vec2(:)
            real(iwp) :: dp
            integer :: i
            if (size(vec1)/=size(vec2)) THEN
                WRITE(*,*) 'Error in UEL function:dot_product'
                CALL XIT
            end if
            dp = 0.0_iwp
            do i = 1,size(vec1)
                dp = dp + vec1(i) * vec2(i)
            end do
        END FUNCTION dot_product

        SUBROUTINE fun_c(x,h,p,tc,scheme,y)
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::x(:),h,p,tc
            CHARACTER(LEN=30),INTENT(IN)::scheme 
            REAL(iwp),INTENT(OUT)::y(:)
            REAL(iwp)::hm
            REAL(iwp),PARAMETER::pi=acos(-1.0),two=2.0_iwp,one=1.0_iwp

            if (scheme=='full')THEN
                hm=h-tc
                y=one/sqrt((hm**2*pi**2/p**2*cos(two*pi*x/p)**2) + one)
            else if (scheme=='simplified')THEN
                y=one/sqrt((h**2*pi**2/p**2*cos(two*pi*x/p)**2) + one)
            end if

        END SUBROUTINE fun_c

        SUBROUTINE fun_c_ez(x,h,p,tc,scheme,y)
            !the fx for ez is cos, not sin
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::x(:),h,p,tc
            CHARACTER(LEN=30),INTENT(IN)::scheme 
            REAL(iwp),INTENT(OUT)::y(:)
            REAL(iwp)::hm
            REAL(iwp),PARAMETER::pi=acos(-1.0),two=2.0_iwp,one=1.0_iwp

            if (scheme=='full')THEN
                hm=h-tc
                y=one/sqrt((hm**2*pi**2/p**2*sin(two*pi*x/p)**2) + one)
            else if (scheme=='simplified')THEN
                y=one/sqrt((h**2*pi**2/p**2*sin(two*pi*x/p)**2) + one)
            end if

        END SUBROUTINE fun_c_ez

        SUBROUTINE fun_s(x,h,p,tc,scheme,y)
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::x(:),h,p,tc
            CHARACTER(LEN=30),INTENT(IN)::scheme 
            REAL(iwp),INTENT(OUT)::y(:)
            REAL(iwp)::hm
            REAL(iwp),PARAMETER::pi=acos(-1.0),two=2.0_iwp,one=1.0_iwp

            if (scheme=='full')THEN
                hm=h-tc
                y=(hm*pi/p*cos(two*pi*x/p))/sqrt(((hm**2*pi**2/p**2 *cos(two*x*pi/p)**2)+ one))
            else if (scheme=='simplified')THEN
                y=(h*pi/p*cos(two*pi*x/p))/sqrt(((h**2*pi**2/p**2 *cos(two*x*pi/p)**2)+ one))
            end if            
        END SUBROUTINE fun_s       

        SUBROUTINE fun_z(x,h,p,tc,offset,scheme,y)
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::x(:),h,p,tc,offset
            CHARACTER(LEN=30),INTENT(IN)::scheme 
            REAL(iwp),INTENT(OUT)::y(:)
            REAL(iwp)::hm
            REAL(iwp),PARAMETER::pi=acos(-1.0),two=2.0_iwp

            if (scheme=='full')THEN
                hm=h-tc
                y=offset+hm/two*sin(two*pi*x/p)
            else if (scheme=='simplified')THEN
                y=offset+h/two*sin(two*pi*x/p)
            end if            
        END SUBROUTINE fun_z  

        SUBROUTINE fun_arch(x,h,p,tc,scheme,y)
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::x(:),h,p,tc
            CHARACTER(LEN=30),INTENT(IN)::scheme 
            REAL(iwp),INTENT(OUT)::y(:)
            REAL(iwp)::hm
            REAL(iwp),PARAMETER::pi=acos(-1.0),two=2.0_iwp

            if (scheme=='full')THEN
                hm=h-tc
                y=hm/two+hm/two*-cos(two*pi*x/p)
            else if (scheme=='simplified')THEN
                y=h/two+h/two*-cos(two*pi*x/p)
            end if            
        END SUBROUTINE fun_arch  

        SUBROUTINE thick_modu(xhalf,h,p,tc,scheme,zhalf,chalf,E1,Im,thick_c,Ez)
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::xhalf(:),h,p,tc,zhalf(:),chalf(:),E1,Im,thick_c
            CHARACTER(LEN=30),INTENT(IN)::scheme 
            REAL(iwp),INTENT(OUT)::Ez
            REAL(iwp)::hm,delta
            REAL(iwp),PARAMETER::two=2.0_iwp,four=4.0_iwp
            if (scheme=='full')THEN
                hm=h-tc
                delta=two*trapz(xhalf,(xhalf/two-p/four/hm*zhalf)**2/chalf/E1/Im)
                Ez=Ez+delta*p/hm*hm/thick_c!reciprocal
            else if (scheme=='simplified')THEN
                delta=two*trapz(xhalf,(xhalf/two-p/four/h*zhalf)**2/chalf/E1/Im)
                Ez=Ez+delta*p/h*h/thick_c!reciprocal
            end if            
        END SUBROUTINE thick_modu

        FUNCTION linspace(start, end, num_points) result(sequence)
            implicit none
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            real(iwp), intent(in) :: start, end
            integer, intent(in) :: num_points
            real(iwp), dimension(:), allocatable :: sequence
            integer :: i
            real(iwp) :: step

            if (num_points <= 0) then
                ERROR STOP 'error in function linspace' 
            end if
            allocate(sequence(num_points))
            step = (end - start) / (num_points - 1)
            do i = 1, num_points
                sequence(i) = start + (i - 1) * step
            end do
        END FUNCTION linspace

        FUNCTION trapz(x, y) result(integral)
            implicit none
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            real(iwp), dimension(:), intent(in) :: x, y
            real(iwp) :: integral
            integer :: n, i
    
            if (size(x) /= size(y)) then
                ERROR STOP "Error in function trapz: x and y must have the same length."
            end if
            n = size(x)
    
            if (n < 2) then
                ERROR STOP "Error in function trapz: At least two points are required."
                return
            end if    
            integral = 0.0
    
            do i = 1, n - 1
                integral = integral + 0.5 * (x(i+1) - x(i)) * (y(i) + y(i+1))
            end do
        END FUNCTION trapz

        FUNCTION expandvector21(vector) result(matrix)
            implicit none
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            real(iwp):: vector(21),matrix(6,6)
            INTEGER i,j
            matrix=0.0_iwp
            do i=1,6
                do j=1,6
                    if (i<=j) then
                        matrix(i,j)=vector((j-1)*j/2+i)
                        matrix(j,i)=matrix(i,j)
                    end if
                end do
            end do
            return
        END FUNCTION

        SUBROUTINE check_symmetric(mat)
            implicit none
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            real(iwp):: mat(:,:)
            integer :: size1, size2,i, j
            size1=UBOUND(mat,1)
            size2=UBOUND(mat,2)
            if (size1/=size2) then 
                WRITE(*,*) 'matrix dimension error'
                call xit
            end if
            do i = 1, size1
               do j = 1, size2
                  if (mat(i, j) /= mat(j, i)) then
                    WRITE(*,*) 'unsymmetric matrix error'
                    WRITE(*,*) i,j
                    call XIT
                     return
                  end if
               end do
            end do
        END SUBROUTINE check_symmetric

        SUBROUTINE matprint(mat,fname)
            implicit none
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            real(iwp):: mat(:,:)
            integer :: size1,size2, i, j,LENOUTDIR
            CHARACTER(LEN=256)::filename,OUTDIR
            CHARACTER(LEN=*)::fname
            CALL GETOUTDIR(OUTDIR,LENOUTDIR)
            size1=UBOUND(mat,1)
            size2=UBOUND(mat,2)
            filename = trim(OUTDIR)//'\'//trim(fname)
            open(unit=6, file=filename, status='REPLACE', action='write')
            do i = 1, size1
               do j = 1, size2
                    write(6,*) i,' ',j,' ',mat(i, j)
               end do
            end do
            close(6)

        END SUBROUTINE matprint
    

End Subroutine uel

SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
    INCLUDE 'ABA_PARAM.INC'
    DIMENSION U(3),TIME(3),COORDS(3)
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    !user coding
    RETURN
END SUBROUTINE

SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,COORDS,JLTYP,SNAME)
     INCLUDE 'ABA_PARAM.INC'
     DIMENSION TIME(2), COORDS (3)
     CHARACTER*80 SNAME
     !user coding to define F
     RETURN
END SUBROUTINE

