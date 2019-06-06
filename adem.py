# from sys import argv
import math
from functools import reduce

def choose(n, k):
    """
    Standard combinatorial choose function
    """

    if n<k:
        return 0
    # there will be k multiplications for both numerator and denominator
    # since choose(n,k) = choose(n,n-k), take the smaller of k and n-k
    k = min(k, n-k)
    numerator = reduce(lambda x,y: x*y, range(n-k+1, n+1), 1)
    denominator = reduce(lambda x,y: x*y, range(1,k+1), 1)
    return int(numerator / denominator)


def sum_splice(mono, splice_list, index1, index2):
    """
    mono: monomial given as List[int]
    splice_list: a List[List[int]] thought of as a sum of monomials
    index1: starting index for splice
    index2: ending index for splice

    INPUT:
        a monomial of Steenrod squares, and a sum of monomials

    OUPUT:
        a sum of monomials, created in the following way:
            for each monomial in splice_list, called new_mono, replace mono between index1 and index2 with new_mono
            add all these spliced monos to a list

    EXAMPLE:
        mono = [1,2,3,4]
        splice_list = [[6,7], [8]]
        index1 = 1
        index2 = 3
        output = [[1,6,7,4], [1,8,4]]
    """

    output = []
    for summand in splice_list:
        new_mono = mono[:]
        new_mono[index1:index2] = summand
        output.append(new_mono)
    return output


def adem(mono):
    """
    mono: length 2 monomial given as [i,j]

    INPUT:
        a length two monomial, Sq^i Sq^j

    OUTPUT:
        a sum of monomials given by applying the Adem relations
        an empty list signifies ZERO
    """

    output = []
    i = mono[0]
    j = mono[1]
    for k in range(int(math.floor(i/2))):
        if (choose(j-k-1, i-2*k) % 2) != 0:
            new_mono = [i+j-k, k]
            output.append(new_mono)
    return output


def write_as_basis(input):
    """
    input: sum of monomials given as List[List[int]]

    INPUT:
        a sum of monomials of Steenrod squares

    OUTPUT:
        a sum of monomials that is the input written in the Serre-Cartan basis, i.e., `admissible form'
        that is, apply the Adem relations to the original sum

    EXAMPLE:
         input = [[2,4,2], [1,2]] representing Sq^2 Sq^4 Sq^2 + Sq^1 Sq^2
         output = [[6,2], [3]] representing Sq^6 Sq^2 + Sq^3

    NOTE:
        the input (and hence the output) need not be homogeneous (i.e. a sum of same-degree elements)
    """

    if len(input) == 0:
        return input

    # if the input is a sum of more than one monomial, process each part separately
    if len(input) > 1:
        output = []
        for mono in input:
            output.extend(write_as_basis([mono]))
        return output

    mono = input[0]

    # if the only monomial is of the form Sq^i, we're already done
    if len(mono) == 1:
        return input

    # reading from left to right in pairs of two, find the first pair of `inadmissible' squares
    # apply the Adem relations to this pair and splice it back in to the sum
    for index, power in enumerate(mono):
        if index > 0 and (prev_power < 2*power):
            adem_sum = adem([prev_power, power])
            if not adem_sum:
                return []
            return write_as_basis(sum_splice(mono, adem_sum, index-1, index+1))
        prev_power = power

    # if we made it this far, we are already in admissible form
    return input

#print(write_as_basis([[2,4,2]]))





