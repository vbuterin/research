from homomorphic_encryption import next_prime, random_vector, encrypt, decrypt, \
    binary_encode, encoded_add, mk_transit_key, error_digits, binary_decrypt, \
    multiply_ciphertexts, generate_key, bit_length, \
    three_to_two, encoded_add3, multi_add, add_and_return_remainder_parity_and_rangecheck, \
    partial_decrypt, mk_bootstrapping_key, bootstrap, flatten_ciphertext, flatten_key
import random

short_modulus = 167
medium_modulus = 10**8
large_modulus = 10**30 + 57
huge_modulus = 10**100 + 267

def generate_keys(modulus):
    s = generate_key(6, modulus)
    zero = encrypt(s, 0, modulus)
    one = encrypt(s, 1, modulus)
    tk = mk_transit_key(s, s, modulus)
    return s, zero, one, tk

def test_add(x, y, keys):
    print("Testing addition circuit for {} and {}".format(x,y))
    S, zero, one, tk = keys
    bitcount = bit_length(x+y)
    encx = binary_encode(x, bitcount, zero, one)
    ency = binary_encode(y, bitcount, zero, one)
    encz = encoded_add(encx, ency, tk)
    print('error digits', [error_digits(S, x) for x in encz])
    assert binary_decrypt(S, encz[:bitcount]) == x + y

def test_add3(x, y, z, keys):
    print("Testing addition circuit for {} {} {}".format(x,y, z))
    S, zero, one, tk = keys
    encx = binary_encode(x, 8, zero, one)
    ency = binary_encode(y, 8, zero, one)
    encz = binary_encode(z, 8, zero, one)
    enco = encoded_add3(encx, ency, encz, tk)
    print('error digits', [error_digits(S, x) for x in enco])
    assert binary_decrypt(S, enco[:8]) == x + y + z

def test_three_to_two(x, y, z, keys):
    print("Testing three-to-two for {} {} {}".format(x, y, z))
    S, zero, one, tk = keys
    encx = binary_encode(x, 8, zero, one)
    ency = binary_encode(y, 8, zero, one)
    encz = binary_encode(z, 8, zero, one)
    v, w = three_to_two(encx, ency, encz, zero, tk)
    assert binary_decrypt(S, v) + binary_decrypt(S, w) == x + y + z

def test_multiadd(values, keys):
    print("Testing addition circuit for {}".format(values))
    S, zero, one, tk = keys
    assert sum(values) < 256
    encoded = [binary_encode(v, 8, zero, one) for v in values]
    z = multi_add(encoded, zero, tk)
    print('error digits', [error_digits(S, x) for x in z])
    assert binary_decrypt(S, z[:8]) == sum(values)

def test_add_and_return_remainder_parity_and_rangecheck(values, start, end, keys):
    S, zero, one, tk = keys
    encoded = [binary_encode(v, 6, zero, one) for v in values]
    in_range, parity = add_and_return_remainder_parity_and_rangecheck(encoded, start, end, zero, one, tk)
    print(error_digits(S, in_range), error_digits(S, parity))
    assert decrypt(S, in_range) == (start <= sum(values) < end)
    assert decrypt(S, parity) == (sum(values)-start) % 2

def test():
    s = generate_key(6, medium_modulus)
    sk = mk_transit_key(s, s, medium_modulus)
    c1 = encrypt(s, 0, medium_modulus)
    c2 = encrypt(s, 0, medium_modulus)
    c3 = encrypt(s, 1, medium_modulus)
    c4 = encrypt(s, 1, medium_modulus)
    assert decrypt(s, c1) == decrypt(s, c2) == 0
    assert decrypt(s, c3) == decrypt(s, c4) == 1
    assert decrypt(flatten_key(s, medium_modulus), flatten_ciphertext(c1)) == 0
    assert decrypt(flatten_key(s, medium_modulus), flatten_ciphertext(c2)) == 0
    assert decrypt(flatten_key(s, medium_modulus), flatten_ciphertext(c3)) == 1
    assert decrypt(flatten_key(s, medium_modulus), flatten_ciphertext(c4)) == 1
    assert decrypt(s, c1 + c2) == 0
    assert decrypt(s, c1 + c4) == 1
    assert decrypt(s, c3 + c4) == 0
    assert decrypt(s, c1 + c3) == 1
    assert decrypt(s, c4 + c4) == 0
    assert decrypt(s, multiply_ciphertexts(c1, c2, sk)) == 0
    assert decrypt(s, multiply_ciphertexts(c1, c3, sk)) == 0
    assert decrypt(s, multiply_ciphertexts(c1, c4, sk)) == 0
    assert decrypt(s, multiply_ciphertexts(c3, c3, sk)) == 1
    assert decrypt(s, multiply_ciphertexts(c3, c4, sk)) == 1
    print('ed', error_digits(s, multiply_ciphertexts(c3, c4, sk)))
    print("Basic tests passed")
    x = encrypt(s, 0, medium_modulus)
    for i in range(5):
        print(i, error_digits(s, x), partial_decrypt(s, x))
        x = multiply_ciphertexts(x, x, sk)
    keys = {
        short_modulus: generate_keys(short_modulus),
        medium_modulus: generate_keys(medium_modulus),
        large_modulus: generate_keys(large_modulus),
        huge_modulus: generate_keys(huge_modulus),
    }
    print("Generated keys")
    test_add(42, 69, keys[large_modulus])
    test_add(123456789, 987654321, keys[huge_modulus])
    test_add3(13, 37, 42, keys[large_modulus])
    test_three_to_two(random.randrange(85), random.randrange(85), random.randrange(85), keys[large_modulus])
    test_multiadd([random.randrange(32) for _ in range(8)], keys[large_modulus])
    print("Multiadd tests passed")
    print("Testing add/parity/rangecheck")
    print("42 vs [40, 50]")
    test_add_and_return_remainder_parity_and_rangecheck([7, 11, 24], 40, 50, keys[large_modulus])
    print("42 vs [41, 50]")
    test_add_and_return_remainder_parity_and_rangecheck([7, 11, 24], 41, 50, keys[large_modulus])
    print("42 vs [42, 50]")
    test_add_and_return_remainder_parity_and_rangecheck([7, 11, 24], 42, 50, keys[large_modulus])
    print("42 vs [43, 50]")
    test_add_and_return_remainder_parity_and_rangecheck([7, 11, 24], 43, 50, keys[large_modulus])
    print("42 vs [40, 43]")
    test_add_and_return_remainder_parity_and_rangecheck([7, 11, 24], 40, 43, keys[large_modulus])
    print("42 vs [40, 42]")
    test_add_and_return_remainder_parity_and_rangecheck([7, 11, 24], 40, 42, keys[large_modulus])
    print("42 vs [40, 41]")
    test_add_and_return_remainder_parity_and_rangecheck([7, 11, 24], 40, 41, keys[large_modulus])
    print("Add/parity/rangecheck tests passed")
    print("Starting bootstrap test")
    t, bit0, bit1, tk = keys[huge_modulus]
    m, _, _, _ = keys[medium_modulus]
    s, smallzero, smallone, _ = keys[short_modulus]
    ENCRYPTING_BIT = 1
    bit = bit1 if ENCRYPTING_BIT == 1 else bit0
    for i in range(20):
        bit = multiply_ciphertexts(bit, bit, tk)
        print(error_digits(t, bit))
    assert decrypt(t, bit) == ENCRYPTING_BIT
    tmk = mk_modulus_change_key(t, m, huge_modulus, medium_modulus)
    mediumbit = change_ciphertext_modulus(bit, tmk)
    assert decrypt(m, mediumbit) == ENCRYPTING_BIT
    msk = mk_modulus_change_key(m, s, medium_modulus, short_modulus)
    smallbit = change_ciphertext_modulus(mediumbit, msk)
    assert decrypt(s, smallbit) == ENCRYPTING_BIT
    print("Generating bootstrap key")
    bk = mk_bootstrapping_key(s, t, short_modulus, huge_modulus)
    _, _, _, tk = keys[huge_modulus]
    print(smallbit, s)
    o = bootstrap(smallbit, bk, tk)
    assert decrypt(t, o) == ENCRYPTING_BIT
    print("Bootstrap successful")
    
if __name__ == '__main__':
    test()
