# References

https://realpython.com/python-fractions/

https://realpython.com/python-rounding/

https://realpython.com/python-square-root-function/

$M_0=0$, $D_0=1$, $a_0 = \lfloor \sqrt{d} \rfloor$

For positive real numbers $x$ and $y$, it is true that

$$
x - y < \left\lfloor \frac{x}{y} \right\rfloor y \leq x
$$ \tag{1}

Let $p$ and $q$ be positive integers with $p>q$. Let $v_0 = p$ and $v_1 = q$. By hypothesis, $v_0 > v_1$.

Define $a_0 = \lfloor \frac{v_0}{v_1} \rfloor$ and define $v_2 = v_0 - a_0 v_1$.

By [(1)](#eq-floor) we have $v_0 - v_1 < \lfloor \frac{v_0}{v_1} \rfloor v_1 \leq v_0$, i.e. $v_0 - v_1 < a_0 v_1 \leq v_0$,
so $v_0 - a_0 v_1 - v_1  < 0 \leq v_0 - a_0v_1$,
so $v_2 - v_1 < 0 \leq v_2$. From $v_2-v_1<0$ we get $v_1 > v_2$. From $0 \leq v_2$ we get $v_2 \geq 0$. So $v_1 > v_2 \geq 0$.


For $m > 1$, if $v_m > 0$ define

$$a_{m-1} = \left\lfloor \frac{v_{m-1}}{v_m} \right\rfloor$$

and 

$$\qquad v_{m+1} = v_{m-1} - a_{m-1}v_m.$$

Since $v_m$ is a strictly decreasing sequence of nonnegative integers, there is some $M \geq 2$ such that $v_M >0$ and $v_{M+1}=0$.

For $0 \leq m \leq M-1$,

$$a_m = \left\lfloor \frac{v_m}{v_{m+1}} \right\rfloor\tag{2}$$

and

$$v_m = a_m v_{m+1} + v_{m+2}\tag{3}$$

and

$$v_0 > v_1 > \cdots > v_M > v_{M+1}$$

where $v_0=p$, $v_1=q$, and $v_{M+1}=0$.

Let $d$ be a positive integer such that $d \vert \gcd(p,q)$, so $d \vert v_0$ and $d \vert v_1$. Suppose for $0 \leq m \leq M-1$ that $d \vert v_m$
and $d \vert v_{m+1}$. Using (3) we get $d \vert v_{m+2}$. By induction, for each $0 \leq m \leq M-1$ it holds that $d \vert v_m$ and $d \vert v_{m+1}$. In particular $d \vert v_M$. We have established that $d \vert \gcd(p,q)$ implies $d \vert v_M$.

Let $d$ be a positive integer such that $d \vert v_M$. (3) tells us $v_{M-1} = a_{M-1}v_M + v_{M+1}$. Using $v_{M+1}=0$ we have $v_{M-1} = a_{M-1}v_{M}$ and so $d \vert v_{M-1}$. For $0 \leq m \leq M-1$ suppose that $d \vert v_{M-m}$ and $d \vert v_{M-m-1}$. From (3) it follows that $d \vert v_{M-m+1}$. By induction, for each $0 \leq m \leq M-1$ it holds that $d \vert v_{M-m}$ and $d \vert v_{M-m-1}$. In particular for $m=M-1$ we get that $d \vert v_1$ and
$d \vert v_0$, namely $d \vert p$ and $d \vert q$ and thus $d \vert \gcd(p, q)$. We have established that $d \vert v_M$ implies $d \vert \gcd(p,q)$.

Therefore $v_M = \gcd(p,q)$.

For example let $p=42$ and $q=12$. $v_0=42$ and $v_1=12$. Then $a_0 =  \lfloor \frac{v_0}{v_1} \rfloor = \lfloor \frac{42}{12} \rfloor = 3$ and $v_2 = v_0-a_0v_1 = 42 - (3)(12) = 42 - 36 = 6$. Finally, $a_1 =  \lfloor \frac{v_1}{v_2} \rfloor = \lfloor \frac{12}{6} \rfloor = 2$ and $v_3 = v_1 - a_1 v_2 = 12 - (2)(6) = 0$.

