# Square root

$M_0=0$, $D_0=1$, $a_0 = \lfloor \sqrt{d} \rfloor$

Define

$$M_{k+1} = D_ka_k - M_k,$$

$$D_{k+1} = \dfrac{d - M_{k+1}^2}{D_k},$$

$$a_{k+1} = \left\lfloor \dfrac{M_{k+1} + a_0}{D_{k+1}} \right\rfloor.$$
