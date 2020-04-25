from matrix_fhe import generate_key, encrypt, matrix_multiply, vector_matrix_multiply, \
    decrypt, get_error, matrix_add, bitify, identity, generate_powers_matrix, \
    binary_encrypt, binary_decrypt, encoded_add, multi_add, _xor, _and, _or, \
    two_of_three, multiply_ciphertexts, three_to_two
import random

def test():
    DIMENSION = 5
    PRECISION = 96
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
        #ct3 = matrix_add(ct1, ct2, PRECISION)
        ct3 = _xor(ct1, ct2, PRECISION)
        ct4 = _and(ct1, ct2, PRECISION)
        ct5 = ct4
        for i in range(5):
            ct5 = _and(ct5, ct5, PRECISION)
        print('ct5', get_error(k, ct5, PRECISION))
        assert decrypt(k, ct3, PRECISION) == v1 ^ v2
        #print('ooo', decrypt(k, ct4, PRECISION))
        assert decrypt(k, ct4, PRECISION) == v1 * v2
        assert decrypt(k, ct5, PRECISION) == v1 * v1 * v2
    for i in range(8):
        a, b, c = i%2, (i//2)%2, i//4
        cta, ctb, ctc = (encrypt(k, x, PRECISION) for x in (a, b, c))
        cto = two_of_three(cta, ctb, ctc, PRECISION)
        assert decrypt(k, cto, PRECISION) == (1 if a+b+c >= 2 else 0)
    print("Addition/multiplication tests passed")
    k = generate_key(DIMENSION, PRECISION)
    forty_two = binary_encrypt(k, 42, 8, PRECISION)
    assert binary_decrypt(k, forty_two, PRECISION) == 42
    sixty_nine = binary_encrypt(k, 69, 8, PRECISION)
    one_one_one = encoded_add(forty_two, sixty_nine, PRECISION)
    print([get_error(k, ct, PRECISION) for ct in one_one_one])
    assert binary_decrypt(k, one_one_one[:8], PRECISION) == 111
    print("Simple addition test passed")
    x, y = three_to_two(forty_two, sixty_nine, one_one_one, PRECISION)
    print([get_error(k, ct, PRECISION) for ct in x+y])
    assert binary_decrypt(k, x, PRECISION) + binary_decrypt(k, y, PRECISION) == 222
    two_two_two = multi_add([forty_two, sixty_nine, one_one_one], PRECISION, bits=8)
    print([get_error(k, ct, PRECISION) for ct in two_two_two])
    assert binary_decrypt(k, two_two_two[:8], PRECISION) == 222
    print("Three-item addition test passed")
    values = [random.randrange(1000) for i in range(8)]
    ciphertexts = [binary_encrypt(k, v, 10, PRECISION) for v in values]
    sum_ciphertext = multi_add(ciphertexts, PRECISION, bits=13)
    print([get_error(k, ct, PRECISION) for ct in sum_ciphertext])
    print(binary_decrypt(k, sum_ciphertext, PRECISION), values)
    assert binary_decrypt(k, sum_ciphertext, PRECISION) == sum(values)
    print("Somewhat less simple addition test passed")

if __name__ == '__main__':
    test()
