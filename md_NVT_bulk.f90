
PROGRAM MD_NVT
!********************************************************************!
!                  PROGRMA DINAMICA MOLECULAR                        !
!            ENSEMBLE CANONICO - TERMOSTATO DE NOSE-HOOVER           !
!                                                                    !
! 1) leio do input 'run_poro' o valor de todos os parâmetros para    !
! a simulação                                                        !
!                                                                    !
! 2) Vetor posição das partículas de fluido                          !
! 3) Vetor posição das partículas de poro                            !
! 4) Vetor velocidade das partículas de fluido                       !
! 5) Vetor força sentido pode cada partícula                         !
! 6) Posição usada para iniciar as velocidades com o vínculo de que a!
! velocidade do centro de massa tem que ser nula                     !
! 7) Usado no termostato NH                                          !
! 8) Vetores da funcao de correlacao de pares                        !
!                                                                    !
!               g(r) = correlacao fluido-fluido                      !
!               gfp(r) = correlaca fluido-poro                       !
!                                                                    !
! 9)  Vetores usados no cálculo da difusão                           !
! 10) Vetor usado no parametro de ordem orientacional Q6             !
!                                                                    !
! 11) O correção para o tamanho da caixa quando coloco poros, note   !
! que a densidade tem que ser mantida cte, tanto sem poros como com  !
! poros, para que isso aconteca aumento o volume da caixa por um     !
! fator que e igual a soma do volume de todos os poros. Lembrando que!
! o lado da caixa de simualação na condição bulk é somente           !
!                                                                    !
!                   L=(\frac{N}{\rho})^{1/3}                         !
!                                                                    !
! 12) N^o de graus de liberdade                                      !
!                                                                    !
! 13) Metado do lado da caixa                                        !
!                                                                    !
! 14) Tempo total para equilibração                                  !
!                                                                    !
! 15) Tempo necessario para fazer as medias                          !
!                                                                    !
! 16) Posição de NOSE-HOOVER - variavel virtual usada na lagrangeana !
!                                                                    !
! 17) Velocidade de NOSE-HOOVER - ps - point-s - ds/dt               !
!                                                                    !
! * NOTA: Os INCLUDE's que coloco antes do inicio do programa serve  !
! para a definicao dos MODULE's usados para deixar os vetores,       !
! parametros globais de simulcao e parametros do potencial como      !
! variaveis do tipo comom. Note que nao preciso ficar colocando essas!
! variaveis na declaracao das subraotinas, ela sao usadas e          !
! atualizadas automaticamente                                        !
!********************************************************************!
  IMPLICIT NONE

  INTEGER::ngr,ntel,ntel2,iblm,ipr,inn,seed,switch
  DOUBLE PRECISION::s,ps,delt
  DOUBLE PRECISION::lbox,gg,hbox
  DOUBLE PRECISION::tequil,tmax,t,temp1
  DOUBLE PRECISION::delg,en,vir,etot,vezes
  DOUBLE PRECISION::press,kin,dtime,diff,pi
  DOUBLE PRECISION::sumtemp,sumen,sumetot,rclor

  INTEGER::nsamp,nequil,nmax,npart,nhis,maxtD
  INTEGER::maxt0,it0,ibmax,nblock
  INTEGER,DIMENSION(:),ALLOCATABLE::ibl
  INTEGER,DIMENSION(:,:),ALLOCATABLE::tel
  INTEGER,DIMENSION(:,:),ALLOCATABLE::list


  DOUBLE PRECISION::rho,temp,dt,sigmalj,q
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::x,y,z
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vx,vy,vz
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::fx,fy,fz
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::r2t
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::g
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::bx,by,bz
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::delr2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::vxb,vyb,vzb

  CHARACTER(len=50)::fmt,filename,x1,x2,x3

!  READ(*,*)rho       !densidade                                      !1)
!  READ(*,*)temp      !temperatura
!  READ(*,*)dt        !intevalo de tempo para integração
!  READ(*,*)sigmalj   !diametro do fluido
!  READ(*,*)q         !massa de nose-hoover
!  READ(*,*)nsamp     !intervalo entre impressão de resultados
!  READ(*,*)nequil    !numero de passos de simulação para equilibrar
!  READ(*,*)nmax      !numero total de passos da simulação
!  READ(*,*)npart     !numero de particulas
!  READ(*,*)nhis      !histrograma g(r)
!  READ(*,*)maxtD     !tamanho do vetor
!  READ(*,*)maxt0
!  READ(*,*)it0

  rho=0.25
!  temp=0.42d0
  READ(*,*)temp
  dt=0.001d0
  sigmalj=1.00d0
  q=2.d0
  nsamp=50
  nequil=1000000
  nmax=  5000000
  npart=800
  nhis=250
  maxtD=1000
  maxt0=1000
  it0=10

  fmt='(f5.3)'
  WRITE(x1,fmt)temp
  x1=ADJUSTL(x1)

  fmt='(f5.3)'
  WRITE(x2,fmt)rho
  x2=ADJUSTL(x2)

  filename='gr_rho_'//TRIM(x2)//'_T_'//TRIM(x1)//'.dat'
  OPEN(unit=1,  file=filename,       action='write')

  filename='data_rho_'//TRIM(x2)//'_T_'//TRIM(x1)//'.dat'
  OPEN(unit=7,  file=filename,        action="write")
  WRITE(7,'(6(6a,2x))')'#','temp','press','rho','E_tot','E_kin'
  WRITE(7,'(a)')REPEAT('-', 6*10)


!  OPEN(unit=2,  file="inirandom.xyz",      action="write")
!  OPEN(unit=6,  file="init_conf_read.xyz", action="write")
  OPEN(unit=5,  file="initsolid.xyz",      action="write")
!  OPEN(unit=7,  file="data.dat",           action="write")
  OPEN(unit=12, file="Q6_psi.dat",         action='write')
  OPEN(unit=8,  file="snapshot.xyz",       action="write")
  OPEN(unit=3,  file="temp.dat",           action="write")
  OPEN(unit=9,  file="r2t.dat",            action="write")
  OPEN(unit=4,  file="av_ener.dat",        action="write")

  seed=34234232
  ibmax=20
  nblock=10
  ipr=8

  lbox=(npart/rho)**(1.d0/3.d0)

  ALLOCATE(x(npart),y(npart),z(npart))                                !2)
  ALLOCATE(vx(npart),vy(npart),vz(npart))                             !4)
  ALLOCATE(fx(npart),fy(npart),fz(npart))                             !5)
  ALLOCATE(bx(npart),by(npart),bz(npart))                             !7)
  ALLOCATE(g(nhis))                                                   !8)
  ALLOCATE(vxb(ibmax,nblock,npart),vyb(ibmax,nblock,npart),&
       vzb(ibmax,nblock,npart))                                       !9)
  ALLOCATE(delr2(ibmax,nblock))                                       !9)
  ALLOCATE(ibl(ibmax),tel(ibmax,nblock))                              !9)

  pi=3.141592653589793d0
  gg=3.d0*npart                                                       !12)
  hbox=0.5d0*lbox                                                     !13)
  tequil=nequil*dt                                                    !14)
  tmax=nmax*dt                                                        !15)
  s=0.d0                                                              !16)
  ps=0.d0                                                             !17)
  rclor=3.69d0

  WRITE(*,*)
  WRITE(*,'(''               CORE-SOFTENED FLUID BULK               '')')
  WRITE(*,'(''------------------------------------------------------'')')
  WRITE(*,'(''     Number of particle              ='',i8)')npart
  WRITE(*,'(''     Fluid density                ='',f10.5)')rho
  WRITE(*,'(''     Box length                   ='',f10.5)')lbox
  WRITE(*,'(''     Time step                    ='',f10.5)')dt
  WRITE(*,'(''     Times of sample              ='',f10.5)')tmax
  WRITE(*,'(''     Times of equilibration       ='',f10.5)')tequil
  WRITE(*,'(''     Total time of simuation      ='',f10.5)')tequil+tmax
  WRITE(*,'(''     Temperature                  ='',f10.5)')temp
  WRITE(*,*)
  WRITE(*,'(''                NOSE-HOOVER THERMOSTAT                '')')
  WRITE(*,'(''------------------------------------------------------'')')
  WRITE(*,'(''     Nose Hoover mass             ='',f10.5)')q
  WRITE(*,'(''     Number of deegres of freedom ='',f10.5)')gg
  WRITE(*,*)

  CALL TIMESTAMP

  IF(rclor.GT.hbox)THEN
     PRINT*,'Your cut-off is larger than half of simulation box'
     PRINT*,''
     CALL ALERT(rclor,hbox)
  END IF

  CALL INILATPOROUS(lbox,npart,x,y,z)

  CALL SETVEL(seed,npart,x,y,z,vx,vy,vz,temp)

  CALL DIFFUSION(0,ntel2,dtime,tel,ibl,iblm,ibmax,nblock,dt,  &
     vx,vy,vz,npart,nsamp,delr2,vxb,vyb,vzb)

  CALL RADIAL(0,ngr,delg,g,lbox,npart,nhis,rho)

  CALL FORCE(g,delg,switch,ngr,lbox,en,vir,ipr,list,inn,npart,&
       x,y,z,fx,fy,fz,rclor,nhis,sigmalj)

  CALL THERMAL(g,delg,ngr,lbox,en,vir,etot,tequil,ps,gg,s,    &
       temp1,ipr,list,npart,nhis,dt,fx,fy,fz,rclor,sigmalj,   &
       temp,vx,vy,vz,x,y,z,q)

  PRINT*,'THERMALIZED SYSTEM'
  PRINT*,'CALCULATING AVERAGES'

  CALL SAMPLE(g,delg,ngr,en,etot,vir,lbox,tmax,press,kin,ntel,&
       dtime,delt,ps,gg,s,temp1,ntel2,tel,ibl,iblm,ipr,list,x,&
       y,z,vx,vy,vz,fx,fy,fz,rho,nhis,npart,dt,ibmax,nblock,  &
       nsamp,rclor,temp,sigmalj,q,delr2,vxb,vyb,vzb)

  PRINT*,'AVERAGED'

  CALL RADIAL(2,ngr,delg,g,lbox,npart,nhis,rho)

  PRINT*,'RADIAL DISTRIBUTION FUNCTION FINISHED'

  CALL DIFFUSION(2,ntel2,dtime,tel,ibl,iblm,ibmax,nblock,dt,&
     vx,vy,vz,npart,nsamp,delr2,vxb,vyb,vzb)

  PRINT*,'DIFFUSION FINISHED'

  OPEN(unit=9,file='fim',action='write')

  CLOSE(1)
  CLOSE(2)
  CLOSE(3)
  CLOSE(4)
  CLOSE(5)
  CLOSE(6)
  CLOSE(8)
  CLOSE(9)
  CLOSE(12)

  CALL TIMESTAMP

END PROGRAM MD_NVT

SUBROUTINE ALERT(rclor,hbox)

  IMPLICIT NONE

  DOUBLE PRECISION,INTENT(INOUT)::rclor,hbox

  rclor=hbox
  print*, 'The cut-off distance was changed to:',rclor
  print*, ''

END SUBROUTINE ALERT

SUBROUTINE INIRANDOM(lbox,x,y,z,npart,sigmalj)

  IMPLICIT NONE

  INTEGER::j,i,seed
  INTEGER,INTENT(IN)::npart

  DOUBLE PRECISION:: xn,yn,zn
  DOUBLE PRECISION:: r1,r2,r3
  DOUBLE PRECISION:: normxff,normyff,normzff, normff2

  DOUBLE PRECISION,INTENT(IN)::lbox
  DOUBLE PRECISION,INTENT(IN)::sigmalj
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::x,y,z

  LOGICAL::overlap

  REAL,EXTERNAL::ran2

  j=0
  seed=654

  DO WHILE(j.LT.npart)

     r1=RAN2(seed)
     r2=RAN2(seed)
     r3=RAN2(seed)

     xn=r1*lbox
     yn=r2*lbox
     zn=r3*lbox

     overlap=.FALSE.

     DO i=1,j
        normxff=ABS(x(i)-xn)
        normyff=ABS(y(i)-yn)
        normzff=ABS(z(i)-zn)
        normff2=normxff**2+normyff**2+normzff**2
        IF(normff2.LT.sigmalj*sigmalj)THEN
           overlap=.TRUE.
        END IF
     END DO

     IF(.NOT.overlap)THEN
        x(j+1)=xn
        y(j+1)=yn
        z(j+1)=zn
        j=j+1
     END IF

  END DO

  PRINT*,'SISTEMA INICIADO - RANDOM BULK'

  WRITE(2,'(i4)') npart
  WRITE(2,'(A,3f6.2)')'Atom'
  DO i=1,npart
     WRITE(2,'("N",f12.6,1x,f12.6,1x,f12.6)')x(i),y(i),z(i)
  END DO

END SUBROUTINE INIRANDOM

SUBROUTINE INILATPOROUS(lbox,npart,x,y,z)
!*********************************************************************!
! Inicio o sistema com os poros situados em uma rede cúbica, para     !
! esta configuração o numero de poros tem que ser uma raiz cubica     !
! exata                                                               !
!                                                                     !
! 1) Se o numero de porous não for uma raiz cúbica exata então pego o !
! valor do proximo inteiro. n é o numero de partículas por eixo da    !
! rede cúbia                                                          !
!                                                                     !
! 2) Espaçamento entre as partículas na rede                          !
!                                                                     !
! 3) Uma vez que o meio poros foi distribuido em uma rede cúbica,     !
! tenho que colocar as particulas do fluidos sem que haja superposição!
! logo mesmo procedimento realizado na subrotina anterior é realizado !
!                                                                     !
! 4) A logica por tras do passo abaixo, sao as mesmas da subrotina    !
! anterior                                                            !
!                                                                     !
!*********************************************************************!

  INTEGER::i,j,k,itel,n
  INTEGER,INTENT(IN)::npart

  DOUBLE PRECISION::del,dx,dy,dz
  DOUBLE PRECISION,INTENT(IN)::lbox
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::x,y,z

  n=NINT((npart)**(1.d0/3.d0))                                      !1)
  del=lbox/n
  itel=0
  dx=-del


  DO i=1,n
     dx=dx+del
     dy=-del
     DO j=1,n
        dy=dy+del
        dz=-del
        DO k=1,n
           dz=dz+del
           IF(itel.LT.npart)THEN
              itel=itel+1
              x(itel)=dx
              y(itel)=dy
              z(itel)=dz
           END IF
        END DO
     END DO
  END DO

  PRINT*,'INITIATED SYSTEM - SOLID CONFIGURATION'

  WRITE(5,'(i4)') npart
  WRITE(5,'(A,3f6.2)')'Atom'
  DO i=1,npart
     WRITE(5,'("N",f12.6,1x,f12.6,1x,f12.6)')x(i),y(i),z(i)
  END DO

END SUBROUTINE INILATPOROUS

SUBROUTINE INIT_READ(npart,x,y,z)

  IMPLICIT NONE

  INTEGER::i,npt2
  INTEGER,INTENT(IN)::npart

  DOUBLE PRECISION::r1,r2,r3
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::x,y,z
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::temp2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::val

  ALLOCATE(temp2(3,1000))

  npt2=0

  OPEN(10,file='init_Hden.dat',status='old')

  DO
     READ(10,*,end=300)temp2(1,npt2+1),temp2(2,npt2+1),temp2(3,npt2+1)
     npt2=npt2+1
  END DO
300 CONTINUE
  CLOSE(4)

  ALLOCATE(val(3,npt2))
  val=temp2(:,:npt2)

  DO i=1,npart
     r1=val(1,i)
     r2=val(2,i)
     r3=val(3,i)
     x(i)=r1
     y(i)=r2
     z(i)=r3
  END DO

  WRITE(6,*)npart
  WRITE(6,*)'Atom'
  DO i=1,npart
     WRITE(6,'("N"f12.6,1x,f12.6,1x,f12.6)')x(i),y(i),z(i)
  END DO
  CLOSE(2)

  PRINT*,'INITIATED SYSTEM - READ FILE'

END SUBROUTINE INIT_READ

SUBROUTINE SETVEL(seed,npart,x,y,z,vx,vy,vz,temp)
!*********************************************************************!
! Inicio a velocidade das partículas, respeitando o vínculo de que    !
! a velocidade do centro de massa é nula.                             !
!                                                                     !
! 1) Inicio a velocidade das partículas com uma distribuição          !
! gaussiana ,isso fará com que menos passos de integração sejam       !
! necessários para a termalização do sistema. No equilíbrio a         !
! distribuição de velocidades do fluido é igual a distribuição de     !
! Maxwell-Boltzmann.                                                  !
!                                                                     !
! 2) Componentes x,y,z da velocidade do centro de massa do sistema    !
!                                                                     !
! 3) Velocidade média de cada partícula                               !
!                                                                     !
! 4) Este é o fator de escala das velocidades (T/T(t))^{1/2}, este    !
! fator faz com que as velocidades iniciais sejam aquelas cuja        !
! a soma da a temperatura desejada                                    !
!                                                                     !
! 5) Neste passo, diminuimos o valor da velocidade previamente        !
! sorteado em rangauss pela velocidade média do centro de massa,      !
! então 'shiftamos' este valor pelo fator de escala fs                !
!                                                                     !
! 6) E a posição média das partículas pode ser obtido por meio da     !
! velocidade encontrada no passo anterior                             !
!                                                                     !
!*********************************************************************!
  IMPLICIT NONE

  INTEGER::i
  INTEGER,INTENT(IN)::seed,npart

  DOUBLE PRECISION::sumvx,sumvy,sumvz,sumv2
  DOUBLE PRECISION::fs,temp,dt,r4,r5,r6
  DOUBLE PRECISION,DIMENSION(npart),INTENT(IN)::x,y,z
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::vx,vy,vz
  DOUBLE PRECISION,DIMENSION(npart)::xm,ym,zm

  REAL,EXTERNAL::ran_gauss,ran2

  sumvx=0.d0
  sumvy=0.d0
  sumvz=0.d0

  DO i=1,npart
     vx(i)=ran_gauss(seed)-0.5d0                                      !1)
     vy(i)=ran_gauss(seed)-0.5d0
     vz(i)=ran_gauss(seed)-0.5d0
     sumvx=sumvx+vx(i)                                                !2)
     sumvy=sumvy+vy(i)
     sumvz=sumvz+vz(i)
  END DO

  sumvx=sumvx/npart                                                   !3)
  sumvy=sumvy/npart
  sumvz=sumvz/npart
  sumv2=sumvx*sumvx+sumvy*sumvy+sumvz*sumvz

  fs=SQRT(3.d0*temp/sumv2)                                            !4)

  DO i=1,npart
     vx(i)=(vx(i)-sumvx)*fs                                           !5)
     vy(i)=(vy(i)-sumvy)*fs
     vz(i)=(vz(i)-sumvz)*fs

     xm(i)=x(i)-vx(i)*dt                                              !6)
     ym(i)=y(i)-vy(i)*dt
     zm(i)=z(i)-vz(i)*dt
  END DO

END SUBROUTINE SETVEL

SUBROUTINE RADIAL(switch,ngr,delg,g,lbox,npart,nhis,rho)
!*********************************************************************!
! Aqui calculo a função de distribuição radial, ou também conhecida   !
! de função de correlação de pares, ela me da a probabilidade de      !
! encontrar uma partícula a uma distância 'r' de outra que foi        !
! colocada na origem                                                  !
!                                                                     !
! 1) Quando switch=0 vou zerar as variáveis que irei acumuluar na     !
! rotina de força e calcular o delta da g(r)                          !
!                                                                     !
! 2) Quando switch=2 começo a calcular a g(r), para isso preciso      !
! calcular a expressão \int_0^{\infty}\pho g(r)4\pir^2dr              !
!                                                                     !
! 3) Foi acumulado em g(i) na rotina de força o número de partículas  !
! por fatia do histograma, aqui simplismente normalizo todas as       !
! componentes da g(r)                                                 !
!                                                                     !
! 4) g(r) fluido-fluido                                               !
!                                                                     !
! 5) gr(r) fluido-poro : fator de normalização, lembrando que para    !
! cada partícula tenho N poros, logo npart*nporous termos             !
!                                                                     !
! 6) Aqui mantenho o termo rho em nidp já que sabemos que a integral  !
! de \rho tem que ser N, e os poros não são partículas, logo não      !
! contribuem para a integral                                          !
!                                                                     !
!*********************************************************************!
  IMPLICIT NONE

  INTEGER,INTENT(IN)::nhis,npart
  DOUBLE PRECISION,INTENT(IN)::lbox,rho
  DOUBLE PRECISION,INTENT(INOUT)::delg
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nhis)::g
!  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nhis)::gfp
  DOUBLE PRECISION::xr,yr,zr,vb,nid,pi,r,r2,nidp,lboxp
  DOUBLE PRECISION::rhop,nequiv,vaux
  INTEGER,INTENT(INOUT)::ngr
  INTEGER,INTENT(IN)::switch
  INTEGER::ig,i

  pi=3.141592653589793d0

  IF(switch.EQ.0)THEN                                                 !1)
     ngr=0
     delg=DBLE(lbox/(2.d0*nhis))
     DO i=1,nhis
        g(i)=0.d0
     END DO
  ENDIF

  IF(switch.EQ.2)THEN                                                 !2)

     DO i=1,nhis                                                      !3)
        r=delg*(i+0.5)
        vb=((i+1)**3-i**3)*delg**3
        nid=(4.d0/3.d0)*pi*vb*rho
        g(i)=g(i)/(ngr*npart*nid)                                     !4)
        WRITE(1,*)r,g(i)
     END DO
     CLOSE(1)

  END IF

END SUBROUTINE RADIAL

SUBROUTINE FORCE(g,delg,switch,ngr,lbox,en,vir,ipr,list,inn,npart,&
     x,y,z,fx,fy,fz,rclor,nhis,sigmalj)
!*********************************************************************!
! 1) Calculando a distância entre duas partículas de fluido           !
!                                                                     !
! 2) Convenção de mínima imagem, se uma dada partícua interage com    !
! outra de modo que seu xr é maior que lbox, então copiamos uma       !
! caixa se simulação ao lado e fazer a interação com uma imagem       !
! da partícula na extremidade oposta da caixa. Obser que se xr<0      !
! então o vetor entre as partículas muda de direção.                  !
!                                                                     !
! 3) No cálculo da g(r), nesta etapa após obtido o valor de r entre   !
! duas partículas, calculamos quantos histogramas da g(r) tem lá      !
! dentro, com isso tomamos o seu inteiro a usamos ele para acessar a  !
! componente do vetor g(r). Ao fazer isso várias vezes acessaremos    !
! alguma componente do vetor mais de uma vez, logo ela terá uma valor !
! maior, ou seja mais partículas foram encontradas naquele valor de   !
! distância.                                                          !
!                                                                     !
! 4) Com o valor de distância entre duas partículas estamos aptos a   !
! calcular a forca entre elas. A força é obtida com a derivada        !
! do potencial. Neste caso tenho uma função que dado r é capaz        !
! de calcular a força entre elas e então seu virial                   !
!                                                                     !
! 5) Aqui igualmente ao procedimento anterior calculamos a força em   !
! uma dada particula, porém desta vez a força é gerada pelos poros e  !
! não por outra partícula                                             !
!                                                                     !
! 6) Aqui vou realizar o mesmo procedimento da g(r) mas agora para a  !
! correlação fluido-poro.                                             !
!                                                                     !
! 7) Como as particulas do fluido e dos poros podem ter tamanhos      !
! diferentes temos que usar a regra de mistura                        !
!                                                                     !
! 8) Lembrando que o raio de corte tem que ser multiplicado por       !
! \sigma, nos caso anterior não foi necessário já que o \sigma dos    !
! fluidos é 1                                                         !
!                                                                     !
! 9) As mesmas observações feitas para as interações fluido - fluido  !
! podem ser feitas aqui                                               !
!                                                                     !
!*********************************************************************!
  IMPLICIT NONE

  INTEGER::i,j,ig
  INTEGER,INTENT(INOUT)::ngr
  INTEGER,INTENT(IN)::switch,nhis
  INTEGER,INTENT(IN)::ipr,inn,npart
  INTEGER,DIMENSION(npart)::icont
  INTEGER,DIMENSION(npart,ipr),INTENT(INOUT)::list

  DOUBLE PRECISION::r2,rc2,ff,virf
  DOUBLE PRECISION::xr,yr,zr,en,vir,r
  DOUBLE PRECISION,INTENT(IN)::lbox,delg,rclor,sigmalj
  DOUBLE PRECISION,INTENT(IN),DIMENSION(npart)::x,y,z
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(npart)::fx,fy,fz
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nhis)::g
  DOUBLE PRECISION,DIMENSION(npart,ipr)::RR

  DOUBLE PRECISION,EXTERNAL::lorentz,virlorentz

  en=0.d0
  vir=0.d0
  virf=0.d0
  fx=0.d0
  fy=0.d0
  fz=0.d0
  icont=0

  IF(switch.EQ.1)ngr=ngr+1
  DO i=1,npart-1
     DO j=i+1,npart                                                   !1)
        xr=x(i)-x(j)
        yr=y(i)-y(j)
        zr=z(i)-z(j)
        xr=xr-lbox*NINT(xr/lbox)                                      !2)
        yr=yr-lbox*NINT(yr/lbox)
        zr=zr-lbox*NINT(zr/lbox)
        r2=xr*xr+yr*yr+zr*zr
        r=SQRT(r2)                                                    !3)
        rc2=rclor*rclor

        IF(switch.EQ.1)THEN
!           IF(MOD(inn,4).EQ.0) THEN
!              CALL NEIGHBORS(icont,RR,list,ipr,r,i,j,npart)
!           END IF
           IF(r.LT.lbox/2.d0)THEN
              ig=INT(r/delg)
              g(ig)=g(ig)+2
           END IF
        END IF

        IF(r2.LT.rc2)THEN                                             !4)
           ff=VIRLORENTZ(r,sigmalj)
           vir=vir+ff*r2
           fx(i)=fx(i)+ff*xr
           fy(i)=fy(i)+ff*yr
           fz(i)=fz(i)+ff*zr
           fx(j)=fx(j)-ff*xr
           fy(j)=fy(j)-ff*yr
           fz(j)=fz(j)-ff*zr
           en=en+LORENTZ(r,sigmalj,rclor)
        END IF
     END DO
  END DO

END SUBROUTINE FORCE

SUBROUTINE THERMAL(g,delg,ngr,lbox,en,vir,etot,tequil,ps,gg,s,&
     temp1,ipr,list,npart,nhis,dt,fx,fy,fz,rclor,sigmalj,temp,&
     vx,vy,vz,x,y,z,q)
!*********************************************************************!
! Nesta estapa apenas termalizo o sistema, levo ele a condição de     !
! equilíbrio integrando as equações de movimento                      !
!*********************************************************************!
  IMPLICIT NONE

  INTEGER,INTENT(IN)::ipr,nhis,npart
  DOUBLE PRECISION,DIMENSION(nhis)::g
  DOUBLE PRECISION,DIMENSION(nhis)::gfp
  INTEGER,INTENT(INOUT),DIMENSION(npart,npart)::list
  DOUBLE PRECISION::en,vir,etot,temp1
  DOUBLE PRECISION::ps,gg,s,delg
  DOUBLE PRECISION:: tequil,t,lbox
  DOUBLE PRECISION,INTENT(IN)::dt,rclor,sigmalj,temp,q
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(npart)::x,y,z
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(npart)::vx,vy,vz
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(npart)::fx,fy,fz
  INTEGER::ngr,inn,i,j

  t=0.d0

  DO WHILE (t.lt.tequil)

     CALL INTEGRATENH(1,ps,gg,en,s,temp1,lbox,x,y,z,vx,vy,vz,&
     fx,fy,fz,npart,dt,temp,q)

     CALL FORCE(g,delg,0,ngr,lbox,en,vir,ipr,list,inn,npart, &
          x,y,z,fx,fy,fz,rclor,nhis,sigmalj)

     CALL INTEGRATENH(2,ps,gg,en,s,temp1,lbox,x,y,z,vx,vy,vz,&
     fx,fy,fz,npart,dt,temp,q)

     WRITE(3,*)t,temp1

     t=t+dt

  END DO

END SUBROUTINE THERMAL

SUBROUTINE INTEGRATENH(switch,ps,gg,en,s,temp1,lbox,x,y,z,vx,vy,vz,&
     fx,fy,fz,npart,dt,temp,q)
!*********************************************************************!
! Integração das equações de movimento usando o thermostato de        !
! nosé-hoover. Fiz as contas e mostrei como as equações de movimento  !
! tem que ser alterada para que as médias dos observávies sejam       !
! condizentes com a distribuição de probabilidades do ensemble        !
! canônico                                                            !
!                                                                     !
! 1) Após ter ter obtido as novas posições aplico as condições de     !
! contorno se uma dada partícula teve uma posição tal que x ou y ou z !
! é > lbox ou  < 0 então tenho que entrar com a partícular pelo lado  !
! oposto da caixa                                                     !
!                                                                     !
! 1.1) Se por ventura uma partícula cair em um posição maior que 2*l  !
! as condições de contorno usuais não dão conta de colocar a partícula!
! dentro da caixa novamente. Sendo assim, modifiquei as condições de  !
! contorno para que qualquer que seja a posição da partícula, as      !
! condições de contorno sempre funcionam.                             !
!                                                                     !
! 2) Hamiltoniano de nosé-hoover                                      !
!                                                                     !
! 3) Temperatura cinética                                             !
!                                                                     !
!*********************************************************************!
  IMPLICIT NONE

  INTEGER::i,iter,res
  INTEGER,INTENT(IN)::npart
  INTEGER,INTENT(IN)::switch

  DOUBLE PRECISION,INTENT(IN)::gg,temp,dt,q
  DOUBLE PRECISION,INTENT(INOUT)::ps,s
  DOUBLE PRECISION,INTENT(INOUT)::temp1,en
  DOUBLE PRECISION::sumv2,delt2,delth,erro
  DOUBLE PRECISION::lbox,pso,delps,ri,di,psn,H

  DOUBLE PRECISION,DIMENSION(npart)::bx,by,bz
  DOUBLE PRECISION,DIMENSION(npart)::vxn,vyn,vzn
  DOUBLE PRECISION,DIMENSION(npart)::vxo,vyo,vzo
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::x,y,z
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::vx,vy,vz
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::fx,fy,fz

  LOGICAL::ready

  IF(switch.EQ.1)THEN
     sumv2=0.d0
     delt2=dt*dt/2.d0
     delth=dt/2.d0

     DO  i=1,npart
        x(i)=x(i)+dt*vx(i)+(fx(i)-ps*vx(i))*delt2
        y(i)=y(i)+dt*vy(i)+(fy(i)-ps*vy(i))*delt2
        z(i)=z(i)+dt*vz(i)+(fz(i)-ps*vz(i))*delt2

        IF(x(i).GT.lbox)x(i)=MOD(x(i),lbox)                           !1.1)
        IF(y(i).GT.lbox)y(i)=MOD(y(i),lbox)
        IF(z(i).GT.lbox)z(i)=MOD(z(i),lbox)
        IF(x(i).LT.0)x(i)=lbox+MOD(x(i),lbox)
        IF(y(i).LT.0)y(i)=lbox+MOD(y(i),lbox)
        IF(z(i).LT.0)z(i)=lbox+MOD(z(i),lbox)

        sumv2=sumv2+vx(i)**2+vy(i)**2+vz(i)**2

        vx(i)=vx(i)+(fx(i)-ps*vx(i))*delth
        vy(i)=vy(i)+(fy(i)-ps*vy(i))*delth
        vz(i)=vz(i)+(fz(i)-ps*vz(i))*delth
     END DO
     s=s+ps*dt+(sumv2-gg*temp)*delt2/q
     ps=ps+(sumv2-gg*temp)*delth/q
  END IF

  IF(switch.EQ.2)THEN
     delth=dt/2.d0
     erro = 1.d-10
     sumv2=0.d0
     DO i=1,npart
        vxn(i)=vx(i)
        vyn(i)=vy(i)
        vzn(i)=vz(i)
        sumv2=sumv2+vxn(i)**2+vyn(i)**2+vzn(i)**2
     END DO
     psn=ps
     ready=.FALSE.
     iter=0
     DO WHILE(.NOT.ready.AND.iter.LT.100)
        iter=iter+1
        pso=psn
        delps=0.d0
        DO i=1,npart
           vxo(i)=vxn(i)
           vyo(i)=vyn(i)
           vzo(i)=vzn(i)
           bx(i)=-delth*(fx(i)-pso*vxo(i))-(vx(i)-vxo(i))
           ri=vxo(i)*dt/q
           delps=delps+ri*bx(i)
           by(i)=-delth*(fy(i)-pso*vyo(i))-(vy(i)-vyo(i))
           ri=vyo(i)*dt/q
           delps=delps+ri*by(i)
           bz(i)=-delth*(fz(i)-pso*vzo(i))-(vz(i)-vzo(i))
           ri=vzo(i)*dt/q
           delps=delps+ri*bz(i)
        END DO
        di=-(pso*delth+1.d0)
        delps=delps-di*((-sumv2+gg*temp)*delth/q-(ps-pso))
        delps=delps/(-dt*delth*sumv2/q+di)
        sumv2=0.d0
        DO i=1,npart
           vxn(i)=vxn(i)+(bx(i)+vxo(i)*delps*delth)/di
           vyn(i)=vyn(i)+(by(i)+vyo(i)*delps*delth)/di
           vzn(i)=vzn(i)+(bz(i)+vzo(i)*delps*delth)/di
           sumv2=sumv2+vxn(i)**2+vyn(i)**2+vzn(i)**2
        END DO
        psn=pso+delps
        ready=.true.
        i=0
        DO WHILE(i.LE.npart.AND.ready)
          i=i+1
          IF(i.LE.npart)THEN
             IF(ABS((vxn(i)-vxo(i))/vxn(i)).GT.erro)ready=.FALSE.
             IF(ABS((vyn(i)-vyo(i))/vyn(i)).GT.erro)ready=.FALSE.
             IF(ABS((vzn(i)-vzo(i))/vzn(i)).GT.erro)ready=.FALSE.
          ELSE
             IF(ABS((psn-pso)/psn).GT.erro)ready=.FALSE.
          END IF
       END DO
     END DO
     DO i=1,npart
        vx(i)=vxn(i)
        vy(i)=vyn(i)
        vz(i)=vzn(i)
     END DO
     ps=psn

     H=sumv2/(2.d0*npart)+en/npart+(ps**2*q)/2.d0+gg*temp*s           !2)
     temp1=sumv2/(3.d0*npart)                                         !3)

  END IF

END SUBROUTINE INTEGRATENH

SUBROUTINE SAMPLE(g,delg,ngr,en,etot,vir,lbox,tmax,press,kin,ntel,&
     dtime,delt,ps,gg,s,temp1,ntel2,tel,ibl,iblm,ipr,list,x,y,z,vx&
     ,vy,vz,fx,fy,fz,rho,nhis,npart,dt,ibmax,nblock,nsamp,rclor,  &
     temp,sigmalj,q,delr2,vxb,vyb,vzb)
!*********************************************************************!
! Nesta etapa, após termalizado o sistema, começo a calcular as médias!
!*********************************************************************!
  IMPLICIT NONE

  INTEGER::inn,i,j,ih
  INTEGER,INTENT(INOUT)::ngr,ntel2,ntel,iblm
  INTEGER,INTENT(IN)::ipr,nhis,ibmax,npart,nblock,nsamp
  INTEGER,DIMENSION(npart,ipr),INTENT(INOUT)::list
  INTEGER,DIMENSION(ibmax),INTENT(IN)::ibl
  INTEGER,DIMENSION(ibmax,nblock),INTENT(INOUT)::tel

  DOUBLE PRECISION::sumtemp,sumen,sumetot,sumvir,sumQl,psi2
  DOUBLE PRECISION::sumpsil2,sumpsil,t,etot,vir,pi,psi,Q6
  DOUBLE PRECISION,INTENT(IN)::lbox,gg,temp,dt,rclor,tmax
  DOUBLE PRECISION,INTENT(IN)::rho,delg,delt,sigmalj,q
  DOUBLE PRECISION,INTENT(INOUT)::s,ps,temp1,en,press,kin,dtime

  DOUBLE PRECISION,DIMENSION(nhis),INTENT(INOUT)::g
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::x,y,z
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::vx,vy,vz
  DOUBLE PRECISION,DIMENSION(npart),INTENT(INOUT)::fx,fy,fz
  DOUBLE PRECISION,DIMENSION(ibmax,nblock)::delr2
  DOUBLE PRECISION,DIMENSION(ibmax,nblock,npart)::vxb,vyb,vzb

  t=0.d0
  inn=0
  sumtemp=0.d0
  sumen=0.d0
  sumetot=0.d0
  sumvir=0.d0
  sumQl=0.d0
  sumpsil=0.d0
  ih=0

  DO WHILE(t.LE.tmax)
     inn=inn+1

     IF(MOD(inn,nsamp).EQ.0)THEN
        CALL DIFFUSION(1,ntel2,dtime,tel,ibl,iblm,ibmax,nblock&
             ,dt,vx,vy,vz,npart,nsamp,delr2,vxb,vyb,vzb)
     END IF

     CALL INTEGRATENH(1,ps,gg,en,s,temp1,lbox,x,y,z,vx,vy,vz,&
          fx,fy,fz,npart,dt,temp,q)

     CALL FORCE(g,delg,1,ngr,lbox,en,vir,ipr,list,inn,npart, &
          x,y,z,fx,fy,fz,rclor,nhis,sigmalj)

!     IF(MOD(inn,nsamp).EQ.0)THEN
!        WRITE(4,'(f12.6,1x,f12.6,1x,f12.6,1x,f12.6)')t,en,   &
!             (3/2)*temp1,en+(3/2)*temp1
!     END IF

!     IF(MOD(inn,4).EQ.0)THEN
!        ih=ih+1
!        CALL PARAM_Q6(lbox,list,ipr,Ql,psil,psil2,npart,x,y,z)
!        sumQl=sumQl+Ql
!        sumpsil=sumpsil+psil
!        sumpsil2=sumpsil2+psil2
!     END IF

     CALL INTEGRATENH(2,ps,gg,en,s,temp1,lbox,x,y,z,vx,vy,vz,&
          fx,fy,fz,npart,dt,temp,q)

     t=t+dt
     sumtemp=sumtemp+temp1
     sumen=sumen+en
     sumetot=sumetot+etot
     sumvir=sumvir+vir

!     IF(MOD(inn,4*nsamp).EQ.0)THEN
!        WRITE(8,'(i4)') npart
!        WRITE(8,'(A,3f6.2)')'Atom'
!        DO i=1,npart
!           WRITE(8,'("N",f12.6,1x,f12.6,1x,f12.6)')x(i),y(i),z(i)
!        END DO
!     END IF

  END DO

  psi=sumpsil/ih
  psi2=sumpsil2/ih
  Q6=sumQl/ih
  temp1=sumtemp*dt/t
  en=sumen*dt/(npart*t)
  etot=sumetot*dt/t
  vir=sumvir*dt/t
  Kin=etot-en
  pi=3.141592653589793d0
  press=rho*temp1+vir/(3.d0*lbox**3)

  WRITE(7, '(5(f10.5,2x))')temp1,press,rho,etot,Kin
  WRITE(12,'(5(f10.5,2x))')rho,Q6,temp1,psi,psi2

END SUBROUTINE SAMPLE

SUBROUTINE DIFFUSION(swit,ntel2,dtime,tel,ibl,iblm,ibmax,nblock,dt,&
     vx,vy,vz,npart,nsamp,delr2,vxb,vyb,vzb)

  IMPLICIT NONE

  INTEGER::j,tdifmax,iblock,i,inmax,in,iblm,ii,ihbmax,inp,ib
  INTEGER,INTENT(IN)::swit,ibmax,nblock,npart,nsamp
  INTEGER,INTENT(INOUT)::ntel2
  INTEGER,DIMENSION(ibmax)::ibl
  INTEGER,DIMENSION(ibmax,nblock),INTENT(INOUT)::tel

  DOUBLE PRECISION,INTENT(IN)::dt
  DOUBLE PRECISION,INTENT(INOUT)::dtime
  DOUBLE PRECISION,DIMENSION(npart),INTENT(IN)::vx,vy,vz
  DOUBLE PRECISION,DIMENSION(ibmax,nblock,npart)::vxb,vyb,vzb
  DOUBLE PRECISION,DIMENSION(ibmax,nblock)::delr2
  DOUBLE PRECISION::time,r2,delx,dely,delz,thmax

  tdifmax=1000

  IF(swit.EQ.0) THEN
     ntel2=0
     delr2=0.d0
     dtime=dt*nsamp

     DO ib=1,ibmax
        ibl(ib)=0
        DO j=1,nblock
           tel(ib,j)=0
           delr2(ib,j)=0.d0
           DO i=1,npart
              vxb(ib,j,i)=0.d0
              vyb(ib,j,i)=0.d0
              vzb(ib,j,i)=0.d0
           END DO
        END DO
     END DO
  END IF

  IF(swit.eq.1) then
     ntel2=ntel2+1
     iblm=1
     ii=ntel2/nblock

     DO WHILE(ii.NE.0)
        iblm=iblm+1
        ii=ii/nblock
        time=dtime*(nbLock**(iblm))
        IF(time.GT.tdifmax) ii=0
     END DO

     iblm=MIN(iblm,ibmax)

     DO ib=1,iblm
        iblock=nblock**(ib-1)
        IF (MOD(ntel2,iblock).EQ.0) THEN
           ibl(ib)=ibl(ib)+1
           inmax = MIN(ibl(ib),nblock)
           DO i=1,npart
              IF (ib.EQ.1) THEN
                 delx=vx(i)
                 dely=vy(i)
                 delz=vz(i)
              ELSE
                 delx=vxb(ib-1,1,i)
                 dely=vyb(ib-1,1,i)
                 delz=vzb(ib-1,1,i)
              END IF
              DO in=1,inmax
                 inp=in
                 IF (inmax.EQ.nblock) inp=in+1
                 IF (in.LT.inmax) THEN
                    vxb(ib,in,i)=vxb(ib,inp,i)+delx
                    vyb(ib,in,i)=vyb(ib,inp,i)+dely
                    vzb(ib,in,i)=vzb(ib,inp,i)+delz
!                    print*,vxb(ib,in,i)
                 ELSE
                    vxb(ib,in,i)=delx
                    vyb(ib,in,i)=dely
                    vzb(ib,in,i)=delz
                 END IF
              END DO
              DO in=1,inmax
                 tel(ib,in)=tel(ib,in)+1
                 delr2(ib,in)=delr2(ib,in)+&
                      (vxb(ib,inmax-in+1,i)*dtime)**2.d0+&
                      (vyb(ib,inmax-in+1,i)*dtime)**2.d0+&
                      (vzb(ib,inmax-in+1,i)*dtime)**2.d0
              END DO
           END DO
        END IF
     END DO
  END IF

  IF(swit.EQ.2)THEN
     thmax=0
     ihbmax=0
     DO ib = 1,MIN(ibmax,iblm)
        DO j = 2,MIN(ibl(ib),nblock)
           IF (tel(ib,j).NE.0)THEN
              WRITE(9,*)j*dtime*(nblock**(ib-1)),delr2(ib,j)/tel(ib,j)&
                   ,tel(ib,j)
              IF(j*dtime*(nblock**(ib-1)).GT.thmax) THEN
                 ihbmax=tel(ib,j)
                 thmax=j*dtime*(nbLock**(ib-1))
              END IF
           END IF
        END DO
     END DO
  END IF

END SUBROUTINE DIFFUSION

SUBROUTINE NEIGHBORS(icont,RR,list,ipr,r,i,j,npart)
!********************************************************************!
! Nesta subrotina vou determinar quem são as particulas mais proximas!
! de uma dada particula i.                                           !
!                                                                    !
! 1) Neste ponto preencho a matriz RR com o valores de r da          !
! partícula i com a j. Na matriz list guardo qual o índice das       !
! partículas que estao interagindo                                   !
!                                                                    !
! 2) Observe que só conto até o número de vizinhas (ipr) que vou     !
! querer avaliar o parametro de ordem,                               !
!                                                                    !
! 3) A matriz list tem a forma (exemplo npart=10, ipr=8):            !
!                                                                    !
!         /  2   3   4   5   6   7   8   9 \                         !
!         |  1   3   4   5   6   7   8   9 |                         !
!         |  1   2   4   5   6   7   8   9 |                         !
!         |  1   2   3   5   6   7   8   9 |                         !
! list =  |  1   2   3   4   6   7   8   9 |                         !
!         |  1   2   3   4   5   7   8   9 |                         !
!         |  1   2   3   4   5   6   8   9 |                         !
!         |  1   2   3   4   5   6   7   9 |                         !
!         |  1   2   3   4   5   6   7   8 |                         !
!         \  1   2   3   4   5   6   7   8 /                         !
!                                                                    !
! onde le-se: a linha 1, ou seja particula 1 interage com as         !
! 2,3,...,9 a particula 2 interage com 1,3,4,....9 e assim por diante!
!                                                                    !
! 4) Agora que tenho uma matriz com o índice das particulas que      !
! interagem, e outra com as distâncias entre elas, !vou começar a    !
! ordenar. Faco isso toda vez que chegar ao final de uma linha       !
!                                                                    !
! 5) Se a distancia r_ij entre uma dada particula i e sua vizinha j é!
! menor que de outra partícula j-1 temos que inverter a ordem dos    !
! r's na matriz RR e seus indices na matriz list.                    !
!                                                                    !
! 6) O mesmo raciocinio vale aqui, porem faço isso quando chegar no  !
! final de uma coluna. Note que o procedimento acima mais este abaixo!
! faz com que as matrizes tenham suas coluas colocas em ordem        !
! crescente, e após isso suas linhas                                 !
!                                                                    !
! 7) Quando icont for maior que ipr, tenho que ver ser se as posições!
! ainda não colocadas na mtatriz são menores do que aquela que já    !
! foram odenadas. Note que icont vai ate npart                       !
!                                                                    !
! 8) Aqui estou olhando para posições que não foram colocadas na     !
! matriz e que eventualmente podem ser menores que os elementis já   !
! alocados                                                           !
!                                                                    !
! 9) Comparando com os valores da matriz, se ele for menor eu  troco !
! os elementos                                                       !
!                                                                    !
! 10)Tanto em RR quanto em list                                      !
!                                                                    !
! 11) O mesmo raciocinio vale aqui, porem comparo aqui entre         !
! diferentes linhas da matrix                                        !
!                                                                    !
! 12) Com os valores das prosições mais proximas de uma dada         !
! particula tenho que calcular os valor dos angulos para calcular os !
! harmonicos esfericos. Vamos inicialmente nos preocupar como obter  !
! os harmonicos esfericos a !partir dos polinomios associados de     !
! Legendre                                                           !
!********************************************************************!

  IMPLICIT NONE

  INTEGER::i,j,l,ll
  INTEGER,INTENT(IN)::ipr,npart
  INTEGER,INTENT(INOUT),DIMENSION(npart)::icont
  INTEGER,INTENT(INOUT),DIMENSION(npart,ipr)::list
  DOUBLE PRECISION,DIMENSION(npart,ipr)::RR
  DOUBLE PRECISION::r,aux,iaux,jaux

  icont(i)=icont(i)+1                                                 !1)
  IF(icont(i).LE.ipr)THEN                                             !2)
     RR(i,icont(i))=r
     list(i,icont(i))=j                                               !3)
  END IF

  icont(j)=icont(j)+1
  IF(icont(j).LE.ipr)THEN
     RR(j,icont(j))=r
     list(j,icont(j))=i
  END IF
                                                                      !4)
  IF(icont(i).EQ.ipr)THEN
     DO ll=1,ipr
        DO l=2,ipr
           IF(RR(i,l).LT.RR(i,l-1))THEN                               !5)
              aux=RR(i,l)
              RR(i,l)=RR(i,l-1)
              RR(i,l-1)=aux
              iaux=list(i,l)
              list(i,l)=list(i,l-1)
              list(i,l-1)=iaux
           END IF
        END DO
     END DO
  END IF

  IF(icont(j).EQ.ipr)THEN                                             !6)
     DO ll=1,ipr
        DO l=2,ipr
           IF(RR(j,l).LT.RR(j,l-1))THEN
              aux=RR(j,l)
              RR(j,l)=RR(j,l-1)
              RR(j,l-1)=aux
              jaux=list(j,l)
              list(j,l)=list(j,l-1)
              list(j,l-1)=jaux
           END IF
        END DO
     END DO
  END IF

  IF(icont(i).GT.ipr)THEN                                             !7)
     IF(r.LT.RR(i,ipr))THEN                                           !8)
        RR(i,ipr)=r
        list(i,ipr)=j
        DO ll=1,ipr
           DO l=2,ipr
              IF(RR(i,l).LT.RR(i,l-1))THEN                            !9)
                 aux=RR(i,l)                                          !10)
                 RR(i,l)=RR(i,l-1)
                 RR(i,l-1)=aux
                 iaux=list(i,l)
                 list(i,l)=list(i,l-1)
                 list(i,l-1)=iaux
              ENDIF
           END DO
        END DO
     ENDIF
  ENDIF

  IF(icont(j).GT.ipr)THEN                                             !11)
     IF(r.LT.RR(j,ipr))THEN
        RR(j,ipr)=r
        list(j,ipr)=i
        DO ll=1,ipr
           DO l=2,ipr
              IF(RR(j,l).LT.RR(j,l-1))THEN
                 aux=RR(j,l)
                 RR(j,l)=RR(j,l-1)
                 RR(j,l-1)=aux
                 jaux=list(j,l)
                 list(j,l)=list(j,l-1)
                 list(j,l-1)=jaux
              ENDIF
           END DO
        END DO
     ENDIF
  ENDIF
                                                                      !12)
END SUBROUTINE NEIGHBORS

FUNCTION PALEGENDRE(l,m,x)
!*********************************************************************!
! Nesta etapa vamos implementar (by Num. Recipes) uma subrotina para  !
! gerar os polinomios associados de Legendre, que serão usados para   !
! o calculo dos harmonicos esfericos.                                 !
!                                                                     !
! 1) Aqui vou obter P_m^m(x), visto que pela relação de recorrência   !
! dos polinômios associados de legendre os demais polinômios podem ser!
! obtidos deste.                                                      !
!                                                                     !
! 2) Observe que 2*l-1 para l=m é equivalente a somar 2 a cada parcela!
!                                                                     !
!       P_m^m(x)=(-1)^l(2*l-1)!!*sqrt{(1.d0-x)(1.d0+x)}P_m^m(x)       !
!                                                                     !
!                                                                     !
! 3) Vamos agora dar conta dos casos particulares, que aqui são 2 :   !
! l=m e  l > m                                                        !
!                                                                     !
! 4) Para l=m é justamente o que foi calculado acima                  !
!                                                                     !
! 5) Para l.NE. m a relação de recorrência muda e temos que calcular o!
! termo 'geral' novamente, que é dado por:                            !
!                                                                     !
!                        P_{m+1}^{m}=x(2m+1)                          !
!                                                                     !
! 6) Caso contrario, calculamos P_l^{m} se l>m+1 onde :               !
!                                                                     !
!            P_l^{m}=A_l^m*\frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l          !
!*********************************************************************!

  IMPLICIT NONE

  INTEGER::i,ll
  INTEGER,INTENT(IN)::l,m
  DOUBLE PRECISION,INTENT(IN)::x
  DOUBLE PRECISION::pmm,fact,square,pll,pmmp1
  DOUBLE PRECISION::palegendre

  IF(m.LT.0.OR.m.GT.l.OR.ABS(x).GT.1.)THEN
     PRINT*,'valor nao podem ser implementado'
  END IF

  pmm=1.d0
  IF(m.GT.0)THEN                                                      !1)
     square=SQRT((1.d0-x)*(1.d0+x))
     fact=1.d0
     DO i=1,m                                                         !2)
        pmm=-pmm*fact*square
        fact=fact+2.d0
     END DO
  END IF
  IF(l.EQ.m)THEN                                                      !3)
     palegendre=pmm                                                   !4)
  ELSE
     pmmp1=x*(2*m+1)*pmm
     IF(l.EQ.m+1)THEN
        palegendre=pmmp1                                              !5)
     ELSE
        DO ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/DBLE(ll-m)             !6)
           pmm=pmmp1
           pmmp1=pll
        END DO
        palegendre=pll
     END IF
  END IF

END FUNCTION PALEGENDRE

SUBROUTINE PARAM_Q6(lbox,list,ipr,Ql,psil,psil2,npart,x,y,z)
!*********************************************************************!
! Nesta etapa vamos calcular o valor média dos harmonicos esféricos   !
! e por fim calcular o calcular o parametro de ordem orientacional.   !
! Para isto vamos usar os polinomios associados de Legendre obtidos   !
! e então os harmonicos esféricos serão obtidos por meio de :         !
!                                                                     !
!           Y_l^m = A_lm*P_l^m*exp(im\phi)               onde         !
!                                                                     !
! 1)      A_l^m = SQRT((2l+1)/(4pi)*((l-m)!/(l+m)!))                  !
!                                                                     !
!                Y_l^{-m} = (-1)^m*(Y_l^m)*                           !
!                                                                     !
! 2) Os polinomios assoc. de Legend. são funções de \cos\theta        !
!                                                                     !
! 3) Lembrando que os harmonicos esféricos é um coeficiente vezes os  !
! polinomios vezes uma exponencial complexa. Neste passo estou somando!
! todas as componento do harmonico esférico.                          !
!                                                                     !
! 4) A quantidade Ylm é real e não complexa, isso justifica o fato de !
! na linha abaixo nao precisar tomar a parte real do numero. Lembre : !
!                                                                     !
!                       z + z* = 2Re(z)                               !
!                                                                     !
! 5) Temos que calcular Y(l,m=0), lembrando que ele não tem parte     !
! complexa logo não segue a mesma forma de transformação dos demais   !
! harmonicos.                                                         !
!                                                                     !
! 6) Guardando o valor do quadrado de Y(l,m=0), note que como esta    !
! é real não preciso tomar o complexo conjugado.                      !
!                                                                     !
! 7) Lembrando que no parametro de ordem temos a soma do valor medio  !
! sobre as m componentes de Ylm e então toma-se a raiz quadrada, e    !
! soma-se sobre todas as particulas                                   !
!                                                                     !
! 8) Feito o passo sete, multiplicamos este resultado pela inverso da !
! da raiz quadrada que aparece em Y(l,m=0)                            !
!                                                                     !
! 9) Finaliza-se dividindo pelo número total de particulas            !
!                                                                     !
! 10) Agora vou começar a calcular o parametro de ordem local \psi_6, !
! note preciso apenas manipular o valor médio do harmonico esferico,  !
! fiz de duas maneiras diferentes, por isso psil e psil2.             !
!                                                                     !
! 11) Lembrando que preciso calcular:                                 !
!                                                                     !
! \psi_6=\sum_{m=-l}^l q_{\ell,m}^iq_{\ell,m}^j \times                !
!     \times \left[ \sum_{m=-l}^l q_{\ell,m}^iq_{\ell,m}^i \times     !
!     \time \sum_{m=-l}^l q_{\ell,m}^jq_{\ell,m}^j \right]^{1/2}      !
!                                                                     !
!*********************************************************************!
  IMPLICIT NONE

  INTEGER::i,l,m,jj,j
  INTEGER,INTENT(IN)::ipr,npart
  INTEGER,INTENT(IN),DIMENSION(npart,ipr)::list
  DOUBLE PRECISION,INTENT(IN)::lbox
  DOUBLE PRECISION,INTENT(INOUT)::Ql
  DOUBLE PRECISION::pi,A_lm, sumYl0,sumYlm,sumYl,bij
  DOUBLE PRECISION::c1_lm,xr,yr,zr,cos0,phi,Qli,sumYli
  DOUBLE PRECISION::r2,r,Ylm,Yl0,psil,psi,den1,den2,num
  DOUBLE PRECISION::den11,den22,num1,bij2,psil2
  DOUBLE PRECISION,DIMENSION(npart),INTENT(IN)::x,y,z
  COMPLEX,DIMENSION(npart)::qlm
  DOUBLE PRECISION,EXTERNAL::FACTORIAL,PALEGENDRE
  COMPLEX::zlm,sumzlm
  COMPLEX,DIMENSION(:,:),ALLOCATABLE::qi

  l=6

  ALLOCATE(qi(l+1,npart))

  pi=3.141592653589793d0
  c1_lm=(2.d0*l+1.d0)/(4.d0*pi)
  sumYli=0.d0
  psil=0.d0
  qlm=0.d0

  DO i=1,npart
     sumYl0=0.d0
     sumYlm=0.d0

     DO jj=1,ipr                                                      !5)
        j=list(i,jj)
        xr=x(i)-x(j)
        yr=y(i)-y(j)
        zr=z(i)-z(j)
        xr=xr-lbox*NINT(xr/lbox)
        yr=yr-lbox*NINT(yr/lbox)
        zr=zr-lbox*NINT(zr/lbox)
        r2=xr*xr+yr*yr+zr*zr
        r=SQRT(r2)
        cos0=zr/r
        Yl0=SQRT(c1_lm)*PALEGENDRE(l,0,cos0)
        sumYl0=sumYl0+Yl0
     END DO
     sumYl0=(sumYl0/DBLE(ipr))**2.                                    !6)
     qi(l+1,i)=SQRT(sumYl0)

     DO m=1,l
        sumYl=0.d0
        A_lm=c1_lm*FACTORIAL(l-m)/(FACTORIAL(l+m))                    !1)
        sumzlm=CMPLX(0.d0,0.d0)

        DO jj=1,ipr
           j=list(i,jj)
           xr=x(i)-x(j)
           yr=y(i)-y(j)
           zr=z(i)-z(j)
           xr=xr-lbox*NINT(xr/lbox)
           yr=yr-lbox*NINT(yr/lbox)
           zr=zr-lbox*NINT(zr/lbox)
           r2=xr*xr+yr*yr+zr*zr
           r=SQRT(r2)
           cos0=zr/r                                                  !2)
           phi=ATAN(yr/xr)

           Ylm=SQRT(A_lm)*PALEGENDRE(l,m,cos0)
           Zlm=CMPLX(COS(m*phi),SIN(m*phi))

           sumzlm=sumzlm+Ylm*zlm                                      !3)
        END DO

        sumzlm=sumzlm/DBLE(ipr)
        Ylm=sumzlm*CONJG(sumzlm)                                      !4)
        sumYlm=sumYlm+2.d0*Ylm

     END DO

     Ylm=sumYlm+sumYl0
     sumYli=sumYli+SQRT(Ylm)                                          !7)

     qlm(i)=sumzlm+SQRT(sumYl0)
     qi(m,i)=sumzlm
  END DO

  DO i=1,npart                                                        !10)
     den1=0.d0
     den2=0.d0
     num=0.d0
     den11=0.d0
     den22=0.d0
     num1=0.d0

     den1=den1+qlm(i)*CONJG(qlm(i))
     den2=den2+qlm(j)*CONJG(qlm(j))
     num=num+REAL(qlm(i)*CONJG(qlm(j)))


     DO jj=1,ipr
        j=list(i,jj)                                                  !11)
        DO m=1,l+1
           den11=den11+SQRT(qi(m,i)*CONJG(qi(m,i)))
           den22=den22+SQRT(qi(m,j)*CONJG(qi(m,j)))
           num1=num1+REAL(qi(m,i)*CONJG(qi(m,j)))
        END DO
     ENDDO

     bij=num/(SQRT(den1)*SQRT(den2))
     bij2=num1/(den11*den22)
     psil=psil+bij/DBLE(ipr)
     psil2=psil2+bij2/DBLE(ipr)
  END DO

  psil=psil/DBLE(npart)
  psil2=psil2/DBLE(npart)
  Qli=sumYli/SQRT(c1_lm)                                              !8)
  Ql=Qli/DBLE(npart)                                                  !9)

END SUBROUTINE PARAM_Q6

FUNCTION FACTORIAL(n)

  IMPLICIT NONE

  INTEGER::i
  INTEGER,INTENT(IN)::n
  DOUBLE PRECISION::factorial

  factorial=1.d0

  DO i=2,n
     factorial=factorial*i
  END DO

END FUNCTION FACTORIAL

SUBROUTINE TIMESTAMP

  IMPLICIT NONE

  CHARACTER (LEN  = 8):: AMPM, DATE
  INTEGER   (KIND = 4):: D, H, M, MM, N, S, VALUES(8), Y
  CHARACTER (LEN  = 9), PARAMETER, DIMENSION (12) :: MONTH = (/ &
&   'January  ', 'February ', 'March    ', 'April    ', &
&   'May      ', 'June     ', 'July     ', 'August   ', &
&   'September', 'October  ', 'November ', 'December ' /)
  CHARACTER (LEN  = 10):: TIME
  CHARACTER (LEN  = 5 ):: ZONE

  CALL DATE_AND_TIME (DATE, TIME, ZONE, VALUES)

  Y  = VALUES(1)
  M  = VALUES(2)
  D  = VALUES(3)
  H  = VALUES(5)
  N  = VALUES(6)
  S  = VALUES(7)
  MM = VALUES(8)


  WRITE(*, *) ' '
  WRITE ( *, '(1x,a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,1x,a3)' ) &
       TRIM(MONTH(M)), D, Y, H, ':', N, ':', S, 'hrs'
  WRITE(*, *) ' '

END SUBROUTINE TIMESTAMP


!*************************************************************************|
!****************************** RANDOM NUMBERS ***************************|
!*************************************************************************!

FUNCTION RAN_GAUSS(seed)
!Gera números aleatórios com uma distribuição gaussiana.
  IMPLICIT NONE

  REAL,PARAMETER:: A1=3.949846138
  REAL,PARAMETER:: A3=0.252408784
  REAL,PARAMETER:: A5=0.076542912
  REAL,PARAMETER:: A7=0.008355968
  REAL,PARAMETER:: A9=0.029899776
  REAL:: gsum,R,R2
  REAL::ran_gauss  !ran2 is a random number
  REAL,EXTERNAL::ran2
  INTEGER:: seed
  INTEGER:: i

  gsum = 0.
  DO i=1,12
     gsum=gsum+ran2(seed)
  END DO
  r=(gsum-6.0)/4.0
  r2=r*r
  ran_gauss=((((A9*R2+A7 )*R2+A5)*R2+A3)*R2+A1)*R
  RETURN

END FUNCTION RAN_GAUSS

FUNCTION RAN2(idum)

  IMPLICIT NONE

  INTEGER::idum,idum2
  INTEGER:: IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL::am,eps,rnmx,ran2
  PARAMETER (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,&
       ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,&
       ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
  INTEGER:: j,k,iv(ntab),iy
  SAVE:: iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  IF (idum.LE.0) THEN
     idum=MAX(-idum,1)
     idum2=idum
     DO j=ntab+8,1,-1
        k=idum/IQ1
        idum=ia1*(idum-k*iq1)-k*ir1
        IF (idum.LT.0) idum=idum+im1
        IF (j.LE.ntab) iv(j)=idum
     END DO
     iy=iv(1)
  END IF
  k=idum/iq1
  idum=ia1*(idum-k*iq1)-k*ir1
  IF (idum.LT.0) idum=idum+im1
  k=idum2/iq2
  idum2=ia2*(idum2-k*iq2)-k*ir2
  IF (idum2.LT.0) idum2=idum2+im2
  j=1+iy/ndiv
  iy=iv(j)-idum2
  iv(j)=idum
  IF(iy.LT.1)iy=iy+imm1
  ran2=MIN(am*iy,rnmx)
  RETURN
END FUNCTION RAN2

SUBROUTINE INIT_RANDOM_SEED()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE INIT_RANDOM_SEED


!**********************************************************************|
!*************************** POTENTIALS & VIRIAL **********************|
!**********************************************************************|
FUNCTION LORENTZ(r,sigmalj,rclor)

  IMPLICIT NONE

  DOUBLE PRECISION,PARAMETER::epslor=1.0d0
  DOUBLE PRECISION,PARAMETER::al =0.30d0
  DOUBLE PRECISION,PARAMETER::dl =1.0d0
  DOUBLE PRECISION,PARAMETER::al2=-1.2d0
  DOUBLE PRECISION,PARAMETER::dl2=1.8d0
  DOUBLE PRECISION,PARAMETER::al3=2.0d0
  DOUBLE PRECISION,PARAMETER::dl3=3.0d0

  DOUBLE PRECISION::lorentzC,lorentz
  DOUBLE PRECISION,INTENT(IN)::r,sigmalj,rclor

  lorentzC=epslor*((sigmalj/rclor)**12.d0-(sigmalj/rclor)**6.d0)+&
       (al/(al**2+(rclor-dl)**2.d0))+&
       (al2/(al2**2+(rclor-dl2)**2.d0))+&
       (al3/(al3**2+(rclor-dl3)**2.d0))

  lorentz=epslor*((sigmalj/r)**12.d0-(sigmalj/r)**6.d0)+&
       (al/(al**2+(r-dl)**2.d0))+&
       (al2/(al2**2+(r-dl2)**2.d0))+&
       (al3/(al3**2+(r-dl3)**2.d0))-lorentzC

  RETURN
END FUNCTION LORENTZ

FUNCTION VIRLORENTZ(r,sigmalj)

  IMPLICIT NONE

  DOUBLE PRECISION,PARAMETER::epslor=1.0d0
  DOUBLE PRECISION,PARAMETER::al =0.30d0
  DOUBLE PRECISION,PARAMETER::dl =1.0d0
  DOUBLE PRECISION,PARAMETER::al2=-1.20d0
  DOUBLE PRECISION,PARAMETER::dl2=1.8d0
  DOUBLE PRECISION,PARAMETER::al3=2.0d0
  DOUBLE PRECISION,PARAMETER::dl3=3.0d0
  DOUBLE PRECISION::virlorentz

  DOUBLE PRECISION,INTENT(IN)::r,sigmalj

  virlorentz=12.d0*epslor*((sigmalj/r)**14.-0.5d0*(sigmalj/r)**8.)+&
       (2/r)*al*(r-dl)/(al**2+(r-dl)**2)**2+    &
       (2/r)*al2*(r-dl2)/(al2**2+(r-dl2)**2)**2+&
       (2/r)*al3*(r-dl3)/(al3**2+(r-dl3)**2)**2

  RETURN
END FUNCTION VIRLORENTZ


!**********************************************************************|
!**********************************************************************|
