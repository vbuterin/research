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
# Output: [[1, 0, 0, 0, 1, 0], [1, 1, 0, 0, 0, 1]]
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
    # Pick a column where we can easily distinguish a small error from an element of key*powers
    col = random.randrange(1, len(key_powers))
    while key_powers[col] >> (precision-2) in (0, 3):
        col = random.randrange(1, len(key_powers))
    # If plaintext=0, prod = small error (at that column)
    # If plaintext=1, prod = key*powers + small error (at that column)
    prod = inner_product([row[col] for row in ct], key, precision)
    return 0 if prod <= key_powers[col]//2 or prod > (key_powers[col]+2**precision)//2 else 1

# Returns the magnitude of the error in the ciphertext (for debugging purposes)
def get_error(key, ct, precision):
    key_powers = vector_matrix_multiply(key, generate_powers_matrix(len(key), precision), precision)
    col = random.randrange(1, len(key_powers))
    while key_powers[col] >> (precision-2) in (0, 3):
        col = random.randrange(1, len(key_powers))
    prod = inner_product([row[col] for row in ct], key, precision)
    target = 0 if prod <= key_powers[col]//2 or prod > (key_powers[col]+2**precision)//2 else key_powers[col]
    return min(abs(prod-target), 2**precision-abs(prod-target))

