from sage.algebras.shuffle_algebra import ShuffleAlgebra
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word

def calclaurent(base, cartan, idx):
    """

    """
    power = 0
    for i, char in enumerate(base):
        if char[1] == 2:
            for j, passing in enumerate(base[i:]):
                if passing[1] == 1:
                    power -= cartan[idx[char[0]], idx[passing[0]]]
    return power


def specinterleave(str1, str2, min_idx=0):
    """

    """
    mylist = []

    if(len(str2) < 2):
        n1 = len(str1)
        n2 = len(str2)
        if n2 > 0: 
            for i in range(min_idx, n1+1):
                mylist.append(str1[0:i] + str2 + str1[i:n1])
        else: 
            mylist = ([str1])
        return mylist

    else:
        mylist = []
        n1 = len(str1)
        n2 = len(str2)
        minvec = range(min_idx, n1+1)
        mychar = str2[0]
        newlist = [str1[0:i] + [mychar] + str1[i:n1] for i in minvec]

        for i, st in enumerate(newlist):
            ret_val = specinterleave(st, str2[1:n2], minvec[i]+1)
            mylist.extend(ret_val)
        return mylist

class QuantumShuffleAlgebra(ShuffleAlgebra):
    """

    """
    def __init__(self, R=LaurentPolynomialRing(QQ,'q'), names='ab', cartan=0):
        """

        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")

        names = Alphabet(names)

        if cartan != 0:
            self._cartan = cartan
            if names.cardinality() != cartan.nrows():
                names = Alphabet([str(i) for i in range(0, cartan.nrows())])
        else:
            self._cartan = CartanMatrix(['A', names.cardinality()])

        ShuffleAlgebra.__init__(self, R, names)

        self._idx = dict(zip(self.variable_names(),range(0, names.cardinality())))
        self._alphabet = Alphabet(names)
        self.__ngens = self._alphabet.cardinality()



    def product_on_basis(self, w1, w2):
        """

        """
        base_multiple = self.base_ring().gen()
        mylist = list()
        name1 = [[x, 1] for x in str(w1)]
        name2 = [[x, 2] for x in str(w2)]
        name3 = specinterleave(name2, name1)
        for name in name3:
            word = ''.join([y[0] for y in name])
            nameval = calclaurent(name, self._cartan, self._idx)
            mylist.append(list((word, base_multiple**nameval)))
        return sum(u[1] * self.basis()[u[0]] for u in mylist)

