# Examples

https://jordanbell.info/LaTeX/euclideanalgorithm/

Aristarchus, On the Sizes and Distances of the Sun and Moon, Proposition 13 [28, p. 397] asserts that 7921 : 4050 > 88 : 45.

T. L. Heath (1913) Aristarchus of Samos: The Ancient Copernicus. Clarendon Press, Oxford.

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

Ptolemy, *Almagest* I.67.22 says the following [17, p. 91, Fragment 41]:

D. R. Dicks (1960) The geographical fragments of Hipparchus. University of London, Athlone Press.

> I have taken the arc from the northernmost limit to the most southerly, that is the arc between the tropics, as being always 47° and more than two-thirds but less than three-quarters of a degree, which is nearly the same estimate as that of Eratosthenes and which Hipparchus also used; for the arc between the tropics amounts to almost exactly 11 of the units of which the meridian contains 83.

Theon of Alexandria, *Commentary on Ptolemy’s Almagest*, writes the following about this passage:

> This ratio is nearly the same as that of Eratosthenes, which Hipparchus also used because it had been accurately measured; for Eratosthenes determined the whole circle as being 83 units, and found that part of it which lies between the tropics to be 11 units; and the ratio 360° : 47°42′40′′ is the same as 83 : 11.

Indeed, we calculate 360° : 47°42′40′′ = 16200 : 2147. The continued fraction of 16200 : 2147 is [7, 1, 1, 5, 195] and the convergents are 7, 3, 15 : 2, 83: 11, 16200 : 2147.

Nicomachus of Gerasa, *Introduction to Arithmetic* I.XIII.10–13 [14, pp. 206–207]:

M. L. D’Ooge, F. E. Robbins, and L. C. Karpinski (1926) Nicomachus of Gerasa: Introduction to arithmetic. University of Michigan Studies, Humanistic Series, Vol. XVI, Macmillan, New York.

> We shall now investigate how we may have a method of discerning whether numbers are prime and incomposite, or secondary and composite, relatively to each other, since of the former unity is the common measure, but of the latter some other number also besides unity; and what this number is.
>
> Suppose there be given us two odd numbers and some one sets the problem and directs us to determine whether they are prime and incomposite relatively to each other or secondary and composite, and if they are secondary and composite what number is their common measure. We must compare the given numbers and subtract the smaller from the larger as many times as possible; then after this subtraction subtract in turn from the other, as many times as possible; for this changing about and subtraction from one and the other in turn will necessarily end either in unity or in some one and the same number, which will necessarily be odd. Now when the subtractions terminate in unity they show that the numbers are prime and incomposite relatively to each other; and when they end in some other number, odd in quantity and twice produced,’ then say that they are secondary and composite relatively to each other, and that their common measure is that very number which twice appears.
>
> For example, if the given numbers were 23 and 45, subtract 23 from 45, and 22 will be the remainder; subtracting this from 23, the remainder is 1, subtracting this from 22 as many times as possible you will end with unity. Hence they are prime and incomposite to one another, and unity, which is the remainder, is their common measure.
>
> But if one should propose other numbers, 21 and 49, I subtract the smaller from the larger and 28 is the remainder. Then again I subtract the same 21 from this, for it can be done, and the remainder is 7. This I subtract in turn from 21 and 14 remains; from which I subtract 7 again, for it is possible, and 7 will remain. But it is not possible to subtract 7 from 7; hence the termination of the process with a repeated 7 has been brought about, and you may declare the original numbers 21 and 49 secondary and composite relatively to each other, and 7 their common measure in addition to the universal unit.

Martianus Capella, *The Marriage of Philology and Mercury*, VII, 785 [61, p. 306]:

W. H. Stahl, R. Johnson, and E. L. Burge (1977) Martianus Capella and the Seven Liberal Arts, volume II. The Marriage of Philology and Mercury. Records of Civilization: Sources and Studies, Columbia University Press.

> If two numbers are composite to one another, a greater and a smaller, how can their largest and their smallest common measure be found? From the larger number let the smaller be subtracted as often as possible; then let whatever amount is left from the former [larger] number be subtracted from the smaller number as often as possible. The amount of the difference will be the greatest measure of these numbers. Take the numbers 350 and 100. Let one hundred be subtracted as often as possible from 350, which is three times. The remainder is 50. From the other number of the pair, one hundred, let 50 be subtracted; the remainder is 50. This number is the greatest common measure of 350 and 100; for fifty times two is one hundred, and fifty times seven is 350. From this calculation it becomes clear how one finds, of all the numbers which measure two numbers, their greatest common measure.

Boethius, *De institutione Arithmetica libri duo*, Book I, Chapter 18 [26, pp. 21–22, §2], “On finding those numbers that are secondary and composite with respect to each other [and numbers that are] prime and incomposite relative to others”:

> The method by which we can find such numbers, if someone proposes them to us and declares that it is not known whether they are commensurable in any measure or [whether] the unit alone measures each, is this. Should two unequal numbers be given, it will be necessary to subtract the smaller from the greater, and if what remains is greater, subtract the smaller from the greater again; [but] if it is smaller, subtract it from the greater [number] that remains, and this should be done until [either] unity finally prevents any further diminution, or, if each of the numbers proposed is odd, some number [is reached that is] necessarily odd; but you will see that the number which is left is equal to that [odd] number. And so it is that if this subtraction should, in turn, reach one, the numbers are said to be prime to each other necessarily and they are conjoined by no other measure except unity alone. If, however, the end of the subtraction arrives at some [odd] number as was said above, it will be a number that measures each sum, and we call the same number that remains the common measure of each.
>
> Take two proposed numbers with respect to which we do not know whether some common measure measures them; let these be 9 and 29. Now we make an alternate subtraction. Let the smaller be subtracted from the greater, that is, 9 from 29, and 20 is left; let us now again subtract the smaller, that is, 9 from 20, and 11 is left; I again subtract from the remainder [i.e., 11] and 2 remains. If I subtract this from 9, 7 is left, and if I again take 2 from 7, 5 remains; and from this another 2 and 3 remains, which after it is diminished by 2 another leaves only unity. Again, if I subtract one from two, the end of the subtraction is fixed at one, which shows that there is no other common measure of these two numbers, namely 9 and 29. Therefore, we will call these numbers prime to each other.
>
> But should other numbers be proposed in the same situation, that is 21 and 9, they could be investigated since they would be mutually related. Again I subtract the quantity of the smaller number from the greater, that is, 9 from 21, and 12 remains. From 12 I take 9 and 3 remains, which if subtracted from 9 leaves 6; and if were taken from 6, 3 would be left, from which 3 cannot be subtracted for it is equal to it. For 3, which was reached by continually subtracting, cannot be subtracted from 3, since they are equal. Therefore, we shall pronounce them commensurable, and 3, which is the remainder, is their common measure.

Itard [38, p. 73] writes:

J. Itard (1961) Les livres arithmétiques d’Euclide. Histoire de la Pensée, Vol. X, Hermann, Paris.

> Cependant on peut trouver quelques démonstrations par récurrence ou induction complète. On ne trouvera jamais le leitmotiv moderne, un peu pédant: «  nous avons vérifié la propriété pour 2, nous avons montré que si elle est vraie pour un nombre, elle est vraie pour son suivant, donc elle est générale »et ceux qui ne voient l’induction complète qu’accompagnée de sa rengaine auront le droit de dire qu’on ne la trouve par dans les Eléments.
>
> Pour nous, nous la voyons dans les prop. 3, 27 et 36, VII, 2, 4 et 13, VIII, 8 et 9, IX.

Netz [49, pp. 268–269]:

R. Netz (1999) The shaping of deduction in Greek mathematics: a study in cognitive history. Ideas in Context, Vol. 51, Cambridge University Press.

> The Greeks cannot speak of ‘A1, A2, ..., An’. What they must do is to use, effectively, something like a dot-representation: the general set of numbers is represented by a diagram consisting of a definite number of lines. Here the generalisation procedure becomes very problematic, and I think the Greeks realised this. This is shown by their tendency to prove such propositions with a number of numbers above the required minimum. This is an odd redundancy, untypical of Greek mathematical economy, and must represent what is after all a justified concern that the minimal case, being also a limiting case, might turn out to be unrepresentative. The fear is justified, but the case of n = 3 is only quantitatively different from the case of n = 2. The truth is that in these propositions Greek actually prove for particular cases, the generalisation being no more than a guess; arithmeticians are prone to guess.
>
> To sum up: in arithmetic, the generalisation is from a particular case to an infinite multitude of mathematically distinguishable cases. This must have exercised the Greeks. They came up with something of a solution for the case of a single infinity. The double infinity of sets of numbers left them defenceless. I suspect Euclid was aware of this, and thus did not consider his particular proofs as rigorous proofs for the general statement, hence the absence of the sumperasma. It is not that he had any doubt about the truth of the general conclusion, but he did feel the invalidity of the move to that conclusion.
>
> The issue of mathematical induction belongs here.
>
> Mathematical induction is a procedure similar to the one described in this chapter concerning Greek geometry. It is a procedure in which generality is sustained by repeatability. Here the similarity stops. The repeatability, in mathematical induction, is not left to be seen by the well-educated mathematical reader, but is proved. Nothing in the practices of Greek geometry would suggest that a proof of repeatability is either possible or necessary. Everything in the practices of Greek geometry prepares one to accept the intuition of repeatability as a substitute for its proof. It is true that the result of this is that arithmetic is less tightly logically principled than geometry – reflecting the difference in their subject matters. Given the paradigmatic role of geometry in this mathematics, this need not surprise us.

Szabó [62] assembles a philological argument that the Euclidean algorithm was created as part of the Pythagorean theory of music. Szabó [62, p. 136, Chapter 2.8] summarizes, “More precisely, this method was developed in the course of experiments with the monochord and was used originally to ascertain the ratio between the lengths of two sections on the monochord. In other words, successive subtraction was first developed in the musical theory of proportions.” Earlier in this work Szabó [62, pp. 28–29] says, “Euclidean arithmetic is predominantly of musical origin not just because, following a tradition developed in the theory of music, it uses straight lines (originally ‘sections of a string’) to symbolize numbers, but also because it uses the method of successive subtraction which was developed originally in the theory of music. However, the theory of odd and even clearly derives from an ‘arithmetic of counting stones’ (ψῆφοι), which did not originally contain the method of successive subtraction.”

Á. Szabó (1978) The beginnings of Greek mathematics. Sythese Historical Library, Vol. 17, D. Reidel Publishing Company, Dordrecht. Translated from the German by A. M. Ungar.


