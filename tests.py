from homomorphic_encryption import next_prime, random_vector, encrypt, decrypt, \
    binary_encode, encoded_add, mk_transit_key, error_digits, binary_decrypt, \
    change_ciphertext_modulus, mk_modulus_change_key, multiply_ciphertexts, \
    three_to_two, encoded_add3, multi_add
import random

short_modulus = 997
medium_modulus = 99991
large_modulus = next_prime(10**100)
huge_modulus = next_prime(10**420)

def test_add(x, y):
    print("Testing addition circuit for {} and {}".format(x,y))
    S = random_vector(5, large_modulus)
    zero = encrypt(S, 0, large_modulus)
    one = encrypt(S, 1, large_modulus)
    encx = binary_encode(x, 8, zero, one)
    ency = binary_encode(y, 8, zero, one)
    encz = encoded_add(encx, ency, mk_transit_key(S, S, large_modulus))
    print('error digits', [error_digits(S, x) for x in encz])
    assert binary_decrypt(S, encz[:8]) == x + y

def test_add3(x, y, z):
    print("Testing addition circuit for {} {} {}".format(x,y, z))
    S = random_vector(5, large_modulus)
    zero = encrypt(S, 0, large_modulus)
    one = encrypt(S, 1, large_modulus)
    encx = binary_encode(x, 8, zero, one)
    ency = binary_encode(y, 8, zero, one)
    encz = binary_encode(z, 8, zero, one)
    enco = encoded_add3(encx, ency, encz, mk_transit_key(S, S, large_modulus))
    print('error digits', [error_digits(S, x) for x in enco])
    assert binary_decrypt(S, enco[:8]) == x + y + z

def test_three_to_two(x, y, z):
    print("Testing three-to-two for {} {} {}".format(x, y, z))
    S = random_vector(5, large_modulus)
    zero = encrypt(S, 0, large_modulus)
    one = encrypt(S, 1, large_modulus)
    encx = binary_encode(x, 8, zero, one)
    ency = binary_encode(y, 8, zero, one)
    encz = binary_encode(z, 8, zero, one)
    v, w = three_to_two(encx, ency, encz, zero, mk_transit_key(S, S, large_modulus))
    assert binary_decrypt(S, v) + binary_decrypt(S, w) == x + y + z

def test_multiadd(values):
    print("Testing addition circuit for {}".format(values))
    assert sum(values) < 256
    S = random_vector(5, huge_modulus)
    zero = encrypt(S, 0, huge_modulus)
    one = encrypt(S, 1, huge_modulus)
    encoded = [binary_encode(v, 8, zero, one) for v in values]
    z = multi_add(encoded, zero, mk_transit_key(S, S, huge_modulus))
    print([error_digits(S, x) for x in z])
    assert binary_decrypt(S, z[:8]) == sum(values)

def test():
    s = random_vector(5, medium_modulus)
    t = random_vector(5, medium_modulus)
    tk = mk_transit_key(s, t, medium_modulus)
    c1 = encrypt(s, 0, medium_modulus)
    c2 = encrypt(s, 0, medium_modulus)
    c3 = encrypt(s, 1, medium_modulus)
    c4 = encrypt(s, 1, medium_modulus)
    assert decrypt(s, c1) == decrypt(s, c2) == 0
    assert decrypt(s, c3) == decrypt(s, c4) == 1
    assert decrypt(s, c1 + c2) == 0
    assert decrypt(s, c1 + c4) == 1
    assert decrypt(s, c3 + c4) == 0
    assert decrypt(t, multiply_ciphertexts(c1, c2, tk)) == 0
    assert decrypt(t, multiply_ciphertexts(c1, c4, tk)) == 0
    assert decrypt(t, multiply_ciphertexts(c3, c4, tk)) == 1
    print("Basic tests passed")
    u = random_vector(3, short_modulus)
    uk = mk_modulus_change_key(s, u, medium_modulus, short_modulus)
    for c in (c1, c2, c3, c4):
        assert decrypt(u, change_ciphertext_modulus(c, uk)) == decrypt(s, c)
    for cl in (c1, c2, c3, c4):
        for cr in (c1, c2, c3, c4):
            cm = multiply_ciphertexts(cl, cr, tk)
            tuk = mk_modulus_change_key(t, u, medium_modulus, short_modulus)
            assert decrypt(u, change_ciphertext_modulus(cm, tuk)) == decrypt(t, cm) == decrypt(s, cl) * decrypt(s, cr)
    print("Modulus change tests passed")
    test_add(42, 69)
    test_add3(13, 37, 42)
    test_three_to_two(random.randrange(85), random.randrange(85), random.randrange(85))
    test_multiadd([random.randrange(32) for _ in range(8)])
    print("Multiadd tests passed")
    
if __name__ == '__main__':
    test()
