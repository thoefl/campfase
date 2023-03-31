! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_oc

implicit none

integer, parameter :: dp = selected_real_kind(15)
integer, parameter :: maxc = 100
real(dp), parameter :: one = 1.0_dp
real(dp), parameter :: zero = 0.0_dp

integer, parameter :: &
       GSBEG=0,       GSOCC=1,        GSADV=2,      GSNOGLOB=3,  &
       GSNOMERGE=4,   GSNODATA=5,     GSNOPHASE=6,  GSNOACS=7,   &
       GSNOREMCS=8,   GSNOSAVE=9,     GSVERBOSE=10, GSSETVERB=11,&
       GSSILENT=12,   GSNOAFTEREQ=13, GSXGRID=14,   GSNOPAR=15,  &
       GSNOSMGLOB=16, GSNOTELCOMP=17, GSTGRID=18,   GSOGRID=19,  &
       GSNORECALC=20, GSOLDMAP=21,    GSNOAUTOSP=22,GSYGRID=23,  &
       GSVIRTUAL=24

TYPE gtp_parerr
    INTEGER :: bmperr
END TYPE gtp_parerr

TYPE(gtp_parerr) :: gx

TYPE tpfun_parres
    integer forcenewcalc
    double precision, dimension(2) :: tpused
    double precision, dimension(6) :: results
END TYPE tpfun_parres

TYPE gtp_components
    integer :: splink,phlink,status
    character(16) :: refstate
    integer, dimension(:), allocatable :: endmember
    double precision, dimension(2) :: tpref
    double precision, dimension(2) :: chempot
    double precision mass,molat
END TYPE gtp_components

TYPE gtp_state_variable
    integer statevarid,norm,unit,phref,argtyp
    integer phase,compset,component,constituent
    double precision coeff
    integer oldstv
end TYPE gtp_state_variable

TYPE gtp_condition
    integer :: noofterms,statev,active,iunit,nid,iref,seqz,experimenttype
    integer symlink1,symlink2
    integer, dimension(:,:), allocatable :: indices
    double precision, dimension(:), allocatable :: condcoeff
    double precision prescribed, current, uncertainty
    TYPE(gtp_state_variable), dimension(:), allocatable :: statvar
    TYPE(gtp_condition), pointer :: next, previous
end TYPE gtp_condition

TYPE gtp_mqmqa_var 
     integer nquad,npair,ns1,ns2
     double precision, allocatable :: yy1(:),yy2(:),dyy1(:,:),dyy2(:,:)
     double precision, allocatable :: d2yy1(:,:),d2yy2(:,:)
     double precision, allocatable :: ceqf1(:),ceqf2(:),dceqf1(:,:),dceqf2(:,:)
     double precision, allocatable :: pair(:),dpair(:,:)
     double precision, allocatable :: eqf1(:),deqf1(:,:),d2eqf1(:,:)
     double precision, allocatable :: eqf2(:),deqf2(:,:),d2eqf2(:,:)
  end type gtp_mqmqa_var

TYPE gtp_fraction_set
     integer latd,ndd,tnoofxfr,tnoofyfr,varreslink,totdis
     character(1) id
     double precision, dimension(:), allocatable :: dsites
     integer, dimension(:), allocatable :: nooffr
     integer, dimension(:), allocatable :: splink
     integer, dimension(:), allocatable :: y2x
     double precision, dimension(:), allocatable :: dxidyj
     double precision fsites
  END TYPE gtp_fraction_set

TYPE gtp_phase_varres
     integer nextfree,phlink,status2,phstate,phtupx
     double precision, dimension(3) :: abnorm
     character(4) prefix,suffix
     integer, dimension(:), allocatable :: constat
     double precision, dimension(:), allocatable :: yfr
     real, dimension(:), allocatable :: mmyfr
     double precision, dimension(:), allocatable :: sites
     double precision, dimension(:), allocatable :: dpqdy
     double precision, dimension(:), allocatable :: d2pqdvay
     type(gtp_fraction_set) :: disfra
     type(gtp_mqmqa_var) :: mqmqaf
     double precision amfu,netcharge,dgm,qcbonds
     double precision, allocatable, dimension(:) :: qcsro
     integer nprop
     integer, dimension(:), allocatable :: listprop
     double precision, dimension(:,:), allocatable :: gval
     double precision, dimension(:,:,:), allocatable :: dgval
     double precision, dimension(:,:), allocatable :: d2gval
     double precision, dimension(3,3) :: curlat
     double precision, dimension(:,:), allocatable :: cinvy
     double precision, dimension(:), allocatable :: cxmol
     double precision, dimension(:,:), allocatable :: cdxmol
     double precision, dimension(:), allocatable :: addg
     integer invsavediter
     double precision, dimension(:,:), allocatable ::invsaved
  END TYPE gtp_phase_varres

type gtp_equilibrium_data
    integer status,multiuse,eqno,nexteq
    character eqname*24,comment*72
    double precision tpval(2),rtn
    double precision :: weight=one
    double precision, dimension(:), allocatable :: svfunres
    TYPE(gtp_condition), pointer :: lastcondition,lastexperiment
    TYPE(gtp_components), dimension(:), allocatable :: complist
    double precision, dimension(:,:), allocatable :: compstoi
    double precision, dimension(:,:), allocatable :: invcompstoi
    TYPE(gtp_phase_varres), dimension(:), allocatable :: phase_varres
    TYPE(tpfun_parres), dimension(:), allocatable :: eq_tpres
    double precision, dimension(:), allocatable :: cmuval
    double precision xconv,gdconv(2)
    double precision :: gmindif
    integer :: maxiter
    integer ::  type_change_phase_amount
    double precision :: scale_change_phase_amount
    integer :: precondsolver
    integer :: splitsolver
    integer :: conv_iter
    character (len=80), dimension(:), allocatable :: eqextra
    integer :: sysmatdim=0,nfixmu=0,nfixph=0
    integer, allocatable :: fixmu(:)
    integer, allocatable :: fixph(:,:)
    double precision, allocatable :: savesysmat(:,:)
    integer eecliq
    double precision eecliqs
    integer, dimension(:,:), allocatable :: phaseremoved
end type gtp_equilibrium_data

interface

    subroutine tqgdmat_th(phase_tuple_idx, temp, pressure, x_known, mu, dmu, current_equi, silent_opt)
        import
        integer, intent(in) :: phase_tuple_idx
        real(dp), intent(in) :: temp, pressure
        real(dp), dimension(:), intent(in) :: x_known
        real(dp), dimension(:), intent(out) :: dmu, mu
        type(gtp_equilibrium_data), pointer, intent(inout) :: current_equi
        logical, intent(in), optional :: silent_opt
    end subroutine tqgdmat_th
    
    subroutine tqsetc(stavar,n1,n2,value,cnam,ceq)
        import
        integer, intent(in) :: n1
        integer, intent(in) :: n2
        character, intent(in) :: stavar*(*)
        character, dimension(maxc), intent(in) :: cnam*24
        real(dp) value
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqsetc
    
    subroutine tqlc(lut,ceq)
        import
        integer lut
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqlc
    
    subroutine tqce(target,n1,n2,ntup,value,ceq,ysave)
        import
        character target*(*)
        integer, intent(in) :: n1,n2
        integer, intent(out) :: ntup
        real(dp) :: value
        type(gtp_equilibrium_data), pointer :: ceq
        real(dp), allocatable, dimension(:,:), intent(inout) :: ysave
    end subroutine tqce
    
    subroutine tqlr(lut,ceq)
        import
        integer lut
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqlr

    
    subroutine tqgetv(stavar,n1,n2,n3,values,ceq,cnam)
        import
        character stavar*(*)
        integer n1, n2, n3
        real(dp) values(*)
        type(gtp_equilibrium_data), pointer, intent(in) :: ceq
        character, dimension(maxc), intent(in) :: cnam*24
    end subroutine tqgetv
    
    subroutine tqgpi(phtupx,phasename,ceq)
        import
        integer phtupx
        character phasename*(*)
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqgpi

    subroutine tqgpi2(iph,ics,phasename,ceq)
        import
        integer iph, ics
        character phasename*(*)
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqgpi2
    
    subroutine tqgphc1(n1,nsub,cinsub,spix,yfrac,sites,extra,ceq)
        import
        integer n1,nsub,cinsub(*),spix(*)
        real(dp) sites(*),yfrac(*),extra(*)
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqgphc1
    
    subroutine get_phase_variance(iph,nv)
        integer iph,nv
    end subroutine get_phase_variance
    
    subroutine get_state_var_value(statevar,value,encoded,ceq)
        import
        character statevar*(*),encoded*(*)
        double precision value
        TYPE(gtp_equilibrium_data), pointer :: ceq
    end subroutine get_state_var_value
    
    subroutine get_sublattice_number(iph,nsl,ceq)
        import
        integer iph,nsl
        TYPE(gtp_equilibrium_data), pointer :: ceq
    end subroutine get_sublattice_number
    
    subroutine tqini(n, ceq)
        import
        integer n
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqini
    
    subroutine tqquiet(yes)
        logical yes
    end subroutine tqquiet
    
    subroutine tqtgsw(i)
        integer i
    end subroutine tqtgsw
    
    subroutine tqrpfil(filename,nel,ntup,selel,ceq,cnam)
        import
        character(*), intent(in) :: filename
        integer, intent(out) :: nel,ntup
        character(len=24), dimension(maxc), intent(out) :: cnam
        character(len=2), dimension(:) :: selel
        type(gtp_equilibrium_data), pointer :: ceq
    end subroutine tqrpfil
    
    subroutine tqphsts(phtupx,newstat,ntup,val,ceq)
        import
        integer, intent(in) :: phtupx,newstat,ntup
        real(dp) val
        type(gtp_equilibrium_data), pointer :: ceq 
    end subroutine tqphsts
    
    subroutine enter_composition_set(iph,prefix,suffix,icsno)
        integer iph,icsno
        character(*) prefix,suffix
    end subroutine enter_composition_set
    
    subroutine ask_default_constitution(cline,last,iph,ics,ceq)
        import
        character cline*(*)
        integer last,iph,ics
        TYPE(gtp_equilibrium_data), pointer :: ceq
    end subroutine ask_default_constitution
    
end interface

end module mod_oc