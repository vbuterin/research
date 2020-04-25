from matrix_fhe import generate_key, encrypt, matrix_multiply, vector_matrix_multiply, \
    decrypt, get_error, matrix_add, bitify, identity, generate_powers_matrix
import random

def test():
    DIMENSION = 5
    PRECISION = 32
    print("Starting basic encryption/decryption tests")
    for i in range(10):
        k = generate_key(DIMENSION, PRECISION)
        value = random.randrange(2)
        ct = encrypt(k, value, PRECISION)
        d = decrypt(k, ct, PRECISION)
        assert d == value
    print("Basic tests passed")
    print("Starting addition/multiplication tests")
    for i in range(10):
        k = generate_key(DIMENSION, PRECISION)
        v1 = random.randrange(2)
        v2 = random.randrange(2)
        ct1 = encrypt(k, v1, PRECISION)
        ct2 = encrypt(k, v2, PRECISION)
        ct3 = matrix_add(ct1, ct2, PRECISION)
        ct4 = matrix_multiply(ct1, bitify(ct2, PRECISION), PRECISION)
        ct5 = ct4
        for i in range(3):
            ct5 = matrix_multiply(ct5, bitify(ct5, PRECISION), PRECISION)
        if not (v1 and v2):
            assert decrypt(k, ct3, PRECISION) == v1 ^ v2
        #print('ooo', decrypt(k, ct4, PRECISION))
        assert decrypt(k, ct4, PRECISION) == v1 * v2
        assert decrypt(k, ct5, PRECISION) == v1 * v1 * v2
        print(list(get_error(k, ct, PRECISION) for ct in (ct3, ct4, ct5)))

if __name__ == '__main__':
    test()

