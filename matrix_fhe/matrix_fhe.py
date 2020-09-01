# Fully homomorphic encryption based on https://eprint.iacr.org/2013/340.pdf and https://eccc.weizmann.ac.il/report/2018/125/

import random

# Error is sampled from range(-ERROR_MAGNITUDE, ERROR_MAGNITUDE+1)
ERROR_MAGNITUDE = 16

# Generates a key with the given length and bits of precision (modulus = 2**precision)
def generate_key(length, precision):
    return [1] + generate_random_vector(length-1, 2**(precision-2), 3 * 2**(precision-2))

# Generates a random vector with the given length, sampled from (range_low, range_high)
def generate_random_vector(length, range_low, range_high):
    return [random.randrange(range_low, range_high) for _ in range(length)]

# Generates a random matrix with the given dimensions, sampled from (range_low, range_high)
def generate_random_matrix(rows, cols, range_low, range_high):
    return [[random.randrange(range_low, range_high) for _ in range(cols)] for _ in range(rows)]

# Adds two or more vectors
def add_vectors(*args):
    inputs, precision = args[:-1], args[-1]
    assert isinstance(precision, int)
    length = len(inputs[0])
    for inp in inputs[1:]:
        assert len(inp) == length
    return [sum(inp[i] for inp in inputs) & (2**precision-1) for i in range(length)]

# M -> kM, self-explanatory
def mul_by_const(M, factor, precision):
    mask = 2**precision-1
    return [[(x*factor) & mask for x in row] for row in M]

# Matrix addition; self-explanatory
def matrix_add(*args):
    inputs, precision = args[:-1], args[-1]
    assert isinstance(precision, int)
    rows, cols = len(inputs[0]), len(inputs[0][0])
    for inp in inputs[1:]:
        assert len(inp) == rows and len(inp[0]) == cols
    return [[sum(inp[i][j] for inp in inputs) & (2**precision-1) for j in range(cols)] for i in range(rows)]

# Matrix multiplication; self-explanatory
def matrix_multiply(A, B, precision):
    assert len(A[0]) == len(B)
    rows = len(A)
    cols = len(B[0])
    C = [[0 for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            C[i][j] = sum([A[i][k] * B[k][j] for k in range(len(B))]) & (2**precision-1)
    return C

# Generates an identity matrix (or with the parameter set a constant-factor scaling matrix)
def identity(rows, mul_factor=1):
    return [[mul_factor if x==y else 0 for x in range(rows)] for y in range(rows)]

# Multiply a vector by a matrix
def vector_matrix_multiply(v, M, precision):
    assert len(M) == len(v)
    return [sum([M[i][j] * v[i] for i in range(len(M))]) & (2**precision-1) for j in range(len(M[0]))]

# Compute the inner product of a vector times a vector
def inner_product(v1, v2, precision):
    assert len(v1) == len(v2)
    return sum([v1[i] * v2[i] for i in range(len(v1))]) & (2**precision-1)

# Splits up each value in a matrix into bits.
# Example: [[1, 2], [3, 4]], precision=3
# Output: [[1, 0], [0, 1], [0, 0], [1, 0], [1, 0], [0, 1]]
def bitify(matrix, precision):
    o = []
    for row in matrix:
        for bit in range(precision):
            o.append([(value >> bit) % 2 for value in row])
    return o

# Generates a matrix of the form eg. [[1, 2, 4, 0, 0, 0], [0, 0, 0, 1, 2, 4]]
# The inverse of bitify (though applicable on inputs that are not possible outputs of bitify)
def generate_powers_matrix(dimension, precision):
    return [[1<<(y%precision) if x==y//precision else 0 for y in range(dimension*precision)] for x in range(dimension)]

# Encrypts a value
def encrypt(key, value, precision):
    dimension = len(key)
    # Step 1: generate a dim x dim matrix A such that key * A = (small error)
    random_matrix = generate_random_matrix(dimension-1, dimension, 0, 2**precision)
    cancelling_value = vector_matrix_multiply([2**precision-x for x in key[1:]], random_matrix, precision)
    error = generate_random_vector(dimension, -ERROR_MAGNITUDE, ERROR_MAGNITUDE+1)
    augmented_matrix = [add_vectors(cancelling_value, error, precision)] + random_matrix
    # Step 2: multiply by random fuzz and get a dim x dim*precision matrix, M = AR
    fuzz = generate_random_matrix(dimension, dimension * precision, -1, 2)
    fuzzed_matrix = matrix_multiply(augmented_matrix, fuzz, precision)
    # If we are encrypting zero, the output is M, where key * M = key * A * R = (key*A) * R = small error
    # If we are encrypting one, the output is M+[powers], where
    # key * M = key * (A * R + powers) = (key * A) * R + key * powers = small error + key * powers
    if value:
        return matrix_add(fuzzed_matrix, generate_powers_matrix(dimension, precision), precision)
    else:
        return fuzzed_matrix

# Decrypts a value, by distinguishing the two above scenarios
def decrypt(key, ct, precision):
    key_powers = vector_matrix_multiply(key, generate_powers_matrix(len(key), precision), precision)
    # The structure of the key ensures that key_powers[col] = 2**precision/2, making it maximally easy
    # to distinguish encryptions of 1 from 0
    col = precision-1
    # If plaintext=0, prod = small error (at that column)
    # If plaintext=1, prod = key*powers + small error (at that column)
    prod = inner_product([row[col] for row in ct], key, precision)
    return 0 if prod >> (precision-2) in (0, 3) else 1

# Returns the log2 of the error in the ciphertext (for debugging purposes)
def get_error(key, ct, precision):
    key_powers = vector_matrix_multiply(key, generate_powers_matrix(len(key), precision), precision)
    col = precision-1
    prod = inner_product([row[col] for row in ct], key, precision)
    return len(bin(min(prod, abs(2**(precision-1) - prod), 2**precision-prod)))-2

# Multiply ciphertexts
def multiply_ciphertexts(A, B, precision):
    o = matrix_multiply(A, bitify(B, precision), precision)
    assert len(o) == len(A) == len(B) and len(o[0]) == len(A[0]) == len(B[0])
    return o

# Encode an integer into a binary representation (least significant bits first)
def binary_encode(integer, length, encoded_zero, encoded_one):
    return [encoded_one if integer & (1 << i) else encoded_zero for i in range(length)]

def binary_encrypt(key, integer, length, precision):
    return [encrypt(key, (integer >> i) % 2, precision) for i in range(length)]

# Decrypt a series of ciphertexts that represent an integer in binary representation
def binary_decrypt(key, output, precision):
    return sum([decrypt(key, o, precision) << i for i,o in enumerate(output)])

# Logical operators
# Note that for all of these operators, you economize on error by putting the
# highest-error argument last
def _and(ct1, ct2, precision):
    return multiply_ciphertexts(ct1, ct2, precision)

def _or(ct1, ct2, precision):
    return matrix_add(ct1, ct2, mul_by_const(multiply_ciphertexts(ct1, ct2, precision), -1, precision), precision)

def _xor(ct1, ct2, precision):
    return matrix_add(ct1, ct2, mul_by_const(multiply_ciphertexts(ct1, ct2, precision), -2, precision), precision)

# 1 if at least two of the inputs are 1, else 0
def two_of_three(ct1, ct2, ct3, precision):
    ab_plus_ac = multiply_ciphertexts(ct1, matrix_add(ct2, ct3, precision), precision)
    bc = multiply_ciphertexts(ct2, ct3, precision)
    minus_two_abc = mul_by_const(multiply_ciphertexts(ct1, bc, precision), -2, precision)
    return matrix_add(ab_plus_ac, bc, minus_two_abc, precision)

# Adds the binary encodings of a and b
def encoded_add(a, b, precision):
    o = []
    carry = mul_by_const(a[0], 0, precision)
    for i in range(len(a)):
        two_of_three_abc = two_of_three(a[i], b[i], carry, precision)
        odd_of_three_abc = matrix_add(a[i], b[i], carry, mul_by_const(two_of_three_abc, -2, precision), precision)
        o.append(odd_of_three_abc)
        carry = two_of_three_abc
    return o + [carry]

# Converts a+b+c into v+w such that a+b+c = v+w. Multiplicative depth 1.
def three_to_two(a, b, c, precision):
    zero = mul_by_const(a[0], 0, precision)
    two_of_three_abc = [two_of_three(ai, bi, ci, precision) for ai, bi, ci in zip(a,b,c)]
    odd_of_three_abc = [matrix_add(ai, bi, ci, mul_by_const(tti, -2, precision), precision) for ai, bi, ci, tti in zip(a, b, c, two_of_three_abc)]
    return (
        odd_of_three_abc + [zero],
        [zero] + two_of_three_abc
    )

# Add together many numbers. Use the 3->2 adder in a tree structure (ok fine it's a DAG),
# then finish off with a 3-to-1 or 2-to-1 as needed
def multi_add(values, precision, bits=999999999999999):
    while len(values) > 2:
        print("Multi adding {} values".format(len(values)))
        o = []
        for i in range(0, len(values)-2, 3):
            x, y = three_to_two(values[i], values[i+1], values[i+2], precision)
            o.extend([x[:bits], y[:bits]])
        o.extend(values[len(values) - len(values) % 3:])
        values = o
    return encoded_add(values[0], values[1], precision)[:bits]
