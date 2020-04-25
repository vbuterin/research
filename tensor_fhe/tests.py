from homomorphic_encryption import random_vector, encrypt, decrypt, \
    binary_encode, encoded_add, mk_transit_key, error_bits, binary_decrypt, \
    multiply_ciphertexts, generate_key, bit_length, adjust_ciphertext_precision, \
    three_to_two, encoded_add3, multi_add, partial_decrypt, mk_bootstrapping_key, \
    bootstrap, flatten_ciphertext, flatten_key
import random, sys, pickle, os

short_precision = 12
medium_precision = 48
large_precision = 112

def generate_keys(precision):
    print("Generating keys for {} bit precision".format(precision))
    s = generate_key(17, precision)
    zero = encrypt(s, 0, precision)
    one = encrypt(s, 1, precision)
    tk = mk_transit_key(s, s, precision)
    return s, zero, one, tk

def generate_all_keys():
    return {
        short_precision: generate_keys(short_precision),
        medium_precision: generate_keys(medium_precision),
        large_precision: generate_keys(large_precision),
        # huge_precision: generate_keys(huge_precision),
    }

def test_add(x, y, keys):
    print("Testing addition circuit for {} and {}".format(x,y))
    S, zero, one, tk = keys
    bitcount = bit_length(x+y)
    encx = binary_encode(x, bitcount, zero, one)
    ency = binary_encode(y, bitcount, zero, one)
    encz = encoded_add(encx, ency, tk)
    print('error bits', [error_bits(S, x) for x in encz])
    assert binary_decrypt(S, encz[:bitcount]) == x + y

def test_add3(x, y, z, keys):
    print("Testing addition circuit for {} {} {}".format(x,y, z))
    S, zero, one, tk = keys
    bitcount = bit_length(x+y+z)
    encx = binary_encode(x, bitcount, zero, one)
    ency = binary_encode(y, bitcount, zero, one)
    encz = binary_encode(z, bitcount, zero, one)
    enco = encoded_add3(encx, ency, encz, tk)
    print('error bits', [error_bits(S, x) for x in enco])
    assert binary_decrypt(S, enco[:bitcount]) == x + y + z

def test_three_to_two(x, y, z, keys):
    print("Testing three-to-two for {} {} {}".format(x, y, z))
    S, zero, one, tk = keys
    bitcount = bit_length(x+y+z)
    encx = binary_encode(x, bitcount, zero, one)
    ency = binary_encode(y, bitcount, zero, one)
    encz = binary_encode(z, bitcount, zero, one)
    v, w = three_to_two(encx, ency, encz, zero, tk)
    assert binary_decrypt(S, v) + binary_decrypt(S, w) == x + y + z

def test_multiadd(values, keys):
    print("Testing addition circuit for {}".format(values))
    S, zero, one, tk = keys
    bitcount = bit_length(sum(values))
    encoded = [binary_encode(v, bitcount, zero, one) for v in values]
    z = multi_add(encoded, zero, tk)
    print('error bits', [error_bits(S, x) for x in z])
    assert binary_decrypt(S, z[:bitcount]) == sum(values)

def test():
    print("Starting basic tests")
    s = generate_key(5, medium_precision)
    sk = mk_transit_key(s, s, medium_precision)
    print("Key generated")
    c1 = encrypt(s, 0, medium_precision)
    c2 = encrypt(s, 0, medium_precision)
    c3 = encrypt(s, 1, medium_precision)
    c4 = encrypt(s, 1, medium_precision)
    assert decrypt(s, c1) == decrypt(s, c2) == 0
    assert decrypt(s, c3) == decrypt(s, c4) == 1
    assert decrypt(flatten_key(s), flatten_ciphertext(c1)) == 0
    assert decrypt(flatten_key(s), flatten_ciphertext(c2)) == 0
    assert decrypt(flatten_key(s), flatten_ciphertext(c3)) == 1
    assert decrypt(flatten_key(s), flatten_ciphertext(c4)) == 1
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
    x = encrypt(s, 0, medium_precision)
    for i in range(5):
        print("Ciphertext after {} rounds of squaring: {}, with {} error bits".format(i, partial_decrypt(s, x), error_bits(s, x)))
        x = multiply_ciphertexts(x, x, sk)
    assert decrypt(s, x) == 0
    assert decrypt(s, adjust_ciphertext_precision(c1, short_precision)) == 0
    assert decrypt(s, adjust_ciphertext_precision(c3, short_precision)) == 1
    assert decrypt(s, adjust_ciphertext_precision(x, short_precision)) == 0
    assert decrypt(s, adjust_ciphertext_precision(c1, large_precision)) == 0
    assert decrypt(s, adjust_ciphertext_precision(c3, large_precision)) == 1
    assert decrypt(s, adjust_ciphertext_precision(x, large_precision)) == 0
    print("Basic tests passed")
    print("Generating more keys")
    keys = generate_all_keys()
    print("Generated keys")
    test_add(42, 69, keys[large_precision])
    test_add3(13, 37, 42, keys[large_precision])
    test_three_to_two(random.randrange(1000), random.randrange(1000), random.randrange(1000), keys[large_precision])
    test_multiadd([random.randrange(1000) for _ in range(8)], keys[large_precision])
    print("Multiadd tests passed")
    print("Starting bootstrap test")
    s, zero, one, tk = keys[large_precision]
    ENCRYPTING_BIT = 1
    bit = one if ENCRYPTING_BIT == 1 else zero
    for i in range(10):
        bit = multiply_ciphertexts(bit, bit, tk)
        print('{} error bits after {} rounds of squaring'.format(error_bits(s, bit), i+1))
    assert decrypt(s, bit) == ENCRYPTING_BIT
    print("Generating bootstrap key")
    bk = mk_bootstrapping_key(s, s, large_precision, short_precision)
    o = bootstrap(bit, bk, tk)
    print("Error bits in bootstrap output: {}".format(error_bits(s, o)))
    assert decrypt(s, o) == ENCRYPTING_BIT
    print("Bootstrap successful")
    
if __name__ == '__main__':
    test()
