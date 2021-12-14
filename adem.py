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
        an empty list signifies the identity, Sq^0
        None signifies zero
    """

    output = []
    i = mono[0]
    j = mono[1]

    if i == 0 and j == 0:
        return []

    for k in range(int(math.floor(i/2))+1):
        if (choose(j-k-1, i-2*k) % 2) != 0:
            # Sq^0 = 1
            if k == 0:
                new_mono = [i+j-k]
            elif i+j-k == 0:
                new_mono = [k]
            else:
                new_mono = [i+j-k, k]
            output.append(new_mono)
    if not output:
        return None
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

    if not input:
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


# reduce mod 2 applied to a list of monomials i.e. a list of lists
# the output is a list of tuples, because there is no need to convert back to lists
def reduce_mod_2(input):
    counts = {}
    for mono in input:
        mono_t = tuple(mono)
        if mono_t in counts:
            counts[mono_t] += 1
        else:
            counts[mono_t] = 1
    reduced_input = []
    for mono_t in counts:
        if counts[mono_t] % 2 != 0:
            reduced_input.append(mono_t)
    
    return reduced_input


# INPUT: sum of Steenrod squares, e.g. 4 2 4 + 1 2 corresponding to Sq^4 Sq^2 Sq^4 + Sq^1 Sq^2
# OUTPUT: the sum as a list of lists e.g. [[4,2,4], [1,2]]
#         return None if the input is not the valid format
def parse_sum_from_string(input):
    string_list = input.split('+')
    mono_list = [mono.split(' ') for mono in string_list] # should contain ints as strings and possibly ""

    mono_list_as_int = []
    for mono_arr in mono_list:
        new_mono = []
        for num in mono_arr:
            if num == "":
                continue
            num = try_parse_int(num)
            if num is None:
                return None
            new_mono.append(num)
        if new_mono:
            mono_list_as_int.append(new_mono)
    
    return mono_list_as_int


def try_parse_int(val):
    try:
        return int(val)
    except ValueError:
        return None


# convert a list of monomials back into the string form 
# ex: [(6,2), (3,)] -> 6 2 + 3
def tuple_list_to_string(input):
    output = ""
    for mono_t in input:
        mono_to_string = ""
        for square in mono_t:
            mono_to_string += str(square) + " "
        output += mono_to_string.strip()
        output += " + "
    
    if output == "":
        return "Zero"
    return output[:-3] # trailing " + "


# apply adem relations to a polynomial of steenrod squares
def apply_adem_to_string(input):
    mono_list = parse_sum_from_string(input)
    if not mono_list or mono_list is None:
        return "invalid format"
    else:
        return tuple_list_to_string(reduce_mod_2(write_as_basis(mono_list)))


# wrapper around apply_adem_to_string in order to print to stdout
def print_adem(input):
    print(apply_adem_to_string(input))
