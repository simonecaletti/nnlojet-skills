* What makefortRR emits for the S,c skeleton next door
* (sc_block_skeleton.map) — shape extracted from real generator output
* and genericized: N = full kin-set index, M = N-1, L = N-2;
* i* = momenta on kin(N), j* on kin(M) (first mapping), k* on kin(L).
* ISOFT is the soft parton k of the skeleton's legend.
*
*     ipass=1
*     ipass_block=1
*     call unset_map()
*     first mapping: the inner cluster (RA,ISOFT,RB) of the iterated
*     product, N -> M; fillsonMFF fills the s{N}on{M} cross-set
*     invariant commons that SS1 reads:
*     call set_map(N,M, (/iRA,iSOFT,iRB/), (/../), ipass)
*     call fillson<M>FF(..., N, ipass)
*     second mapping M -> L (the reduced-ME kinematics):
*     call set_map(M,L, (/../), (/../), ipass)
*     call fillson<L>FF(..., N, ipass)
*     call <MYJET>(0,L,ipass)
*     call set_flav_perm(L, (/ .. /))
*     if(ipass.eq.1)then
*       jpass(NN)=1
*       call getqcdnorm(ix,partons,facnorm, L, ipass_block, IB)
*       if(ipass_block/=0) then
*       wtsoft=0d0
*       one SS1 per SFF of the .map, same signs; arguments are
*       (radiator, soft, radiator, jpset, ipset): radiators on the
*       MAPPED set jpset (j* here), the soft leg on the ORIGINAL set
*       ipset = N — the SSfortset substitution and findipsetSS in the
*       generator pick the set pair per term:
*       wtsoft=wtsoft+1d0*SS1(jW ,iSOFT,jRB,M,N)
*       wtsoft=wtsoft+1d0*SS1(jRA,iSOFT,jV ,M,N)
*       wtsoft=wtsoft-1d0*SS1(jW ,iSOFT,jRA,M,N)
*       wtsoft=wtsoft-1d0*SS1(jRB,iSOFT,jV ,M,N)
*       wt(NN)=-1d0*wtsoft*X30NAME(jW,jRA,jRB,M)*REDME(k..,k..,k..,k1,k2)
*       wt(NN)=bino(ix,partons,-wt(NN)*facnorm)
*       end if
*     endif
*     call unset_flav_perm()
*
* Facts this shape depends on (verified in source):
* - SS1(j1,i3,j2,jpset,ipset) [src/X30/SS1.f] takes radiators j1,j2 on
*   kin(jpset) and the soft leg i3 on kin(ipset); it reads the
*   invariants from the cross-set commons s{ipset}on{jpset} (s7on6,
*   s6on5, ...), which are filled by the fillson* routines emitted
*   with the set_map chain [src/map/libmap.f]. The caller must
*   guarantee that chain ran for the set pair — the generator emits
*   it; a hand-written probe must mimic it.
* - SS(i1,i3,i2,ipset) [src/X30/SS.f] is the unmapped-radiators
*   variant: radiators and the notational soft label live on ONE set,
*   and the actual soft momentum is read from common /soft/psoft(4),
*   filled by makesoft(i1,ipset) [src/map/libmap.f]. makesoft has no
*   caller in the current tree — a probe using SS must call it itself.
* - makefortRR routes a term through this path when its antenna
*   content is SS-set functions (times at most one X30): the branches
*   commented "sum SS * MM0" and "X30 * sum SS * ML0", which toggle
*   the insoft flag and extract the soft leg from the SS middle
*   argument.
