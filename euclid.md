# Examples

https://jordanbell.info/LaTeX/euclideanalgorithm/

Aristarchus, On the Sizes and Distances of the Sun and Moon, Proposition 13 [28, p. 397] asserts that 7921 : 4050 > 88 : 45.

```
x=7921 y=4050
q_list=[None, None, 1, 1, 21, 1, 1, 1, 2, 22]
r_list=[7921, 4050, 3871, 179, 112, 67, 45, 22, 1, 0]
s_list=[1, 0, 1, -1, 22, -23, 45, -68, 181, -4050]
t_list=[0, 1, -1, 2, -43, 45, -88, 133, -354, 7921]

gcd of 7921 and 4050:    1

s * x + t * y = 181 * 7921 + -354 * 4050 = 1

Continued fraction for 7921 / 4050:      [1, 1, 21, 1, 1, 1, 2, 22]

Convergents of continued fraction:
Decimal expansion for 1 / 1:     1
Decimal expansion for 2 / 1:     2
Decimal expansion for 43 / 22:   1.9545454545454545455
Decimal expansion for 45 / 23:   1.9565217391304347826
Decimal expansion for 88 / 45:   1.9555555555555555556
Decimal expansion for 133 / 68:  1.9558823529411764706
Decimal expansion for 354 / 181:         1.9558011049723756906
Decimal expansion for 7921 / 4050:       1.9558024691358024691

Decimal expansion for 7921 / 4050:       1.955802469135802469
```

Proposition 15 [28, p. 407], asserts that 71755875 : 61735500 > 43 : 37.

```
x=71755875 y=61735500
q_list=[None, None, 1, 6, 6, 4, 1, 2, 1, 2, 1, 6]
r_list=[71755875, 61735500, 10020375, 1613250, 340875, 249750, 91125, 67500, 23625, 20250, 3375, 0]
s_list=[1, 0, 1, -6, 37, -154, 191, -536, 727, -1990, 2717, -18292]
t_list=[0, 1, -1, 7, -43, 179, -222, 623, -845, 2313, -3158, 21261]

gcd of 71755875 and 61735500:    3375

s * x + t * y = 2717 * 71755875 + -3158 * 61735500 = 3375

Continued fraction for 71755875 / 61735500:      [1, 6, 6, 4, 1, 2, 1, 2, 1, 6]

Convergents of continued fraction:
Decimal expansion for 1 / 1:     1
Decimal expansion for 7 / 6:     1.1666666666666666667
Decimal expansion for 43 / 37:   1.1621621621621621622
Decimal expansion for 179 / 154:         1.1623376623376623377
Decimal expansion for 222 / 191:         1.1623036649214659686
Decimal expansion for 623 / 536:         1.1623134328358208955
Decimal expansion for 845 / 727:         1.1623108665749656121
Decimal expansion for 2313 / 1990:       1.1623115577889447236
Decimal expansion for 3158 / 2717:       1.1623113728376886272
Decimal expansion for 21261 / 18292:     1.1623113929586704570

Decimal expansion for 71755875 / 61735500:       1.1623113929586704570
```

# Quotations

Boethius De institutione Arithmetica libri duo, Book I, Chapter 18 [26, pp. 21–22, §2], “On finding those numbers that are secondary and composite with respect to each other [and numbers that are] prime and incomposite relative to others”:

> The method by which we can find such numbers, if someone proposes them to us and declares that it is not known whether they are commensurable in any measure or [whether] the unit alone measures each, is this. Should two unequal numbers be given, it will be necessary to subtract the smaller from the greater, and if what remains is greater, subtract the smaller from the greater again; [but] if it is smaller, subtract it from the greater [number] that remains, and this should be done until [either] unity finally prevents any further diminution, or, if each of the numbers proposed is odd, some number [is reached that is] necessarily odd; but you will see that the number which is left is equal to that [odd] number. And so it is that if this subtraction should, in turn, reach one, the numbers are said to be prime to each other necessarily and they are conjoined by no other measure except unity alone. If, however, the end of the subtraction arrives at some [odd] number as was said above, it will be a number that measures each sum, and we call the same number that remains the common measure of each.
