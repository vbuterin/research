# A toy homomorphic encryption implementation, based on
# https://eprint.iacr.org/2012/078.pdf

from dataclasses import dataclass
import random

ERROR_BITS = 6

# Returns the bit length of n, eg. 6 -> 110 -> 3, 31 -> 11111 -> 5, 32 -> 100000 -> 6
def bit_length(n):
    return len(bin(n))-2

# Random length-l vector in a field with the given modulus
def random_vector(l, precision):
    return [random.randrange(2**precision) for _ in range(l)]

# Generates a noise value within some small range
def noise():
    return random.randrange(-2**ERROR_BITS, ERROR_BITS)

# Inner product of two vectors, that is compute sum(vec1[i] * vec2[i]).
# Note that it is done modulo 2**precision.
def prod(vec1, vec2, precision):
    return sum([i*j for i,j in zip(vec1, vec2)]) & (2**precision-1)

# Generates a private key. The first value must always be 1.
def generate_key(length, precision):
    return [1] + [(2**ERROR_BITS + noise()) & (2**precision-1) for _ in range(1, length)]

# A ciphertext representing some value m in {0,1}. A ciphertext is a
# vector v where v . k = m * q/2 + e where . is the inner product,
# q is the modulus (2**precision) and e is the error
@dataclass
class Ciphertext():
    values: list # [int]
    precision: int

    # Add together two ciphertexts into one, linearly adding together the values
    # Note that this combines the magnitudes of the error, so if you add waaaaay
    # too many times the error may overflow
    def __add__(self, other):
        assert self.precision == other.precision and len(self.values) == len(other.values)
        return Ciphertext(
            values = [(x+y) & (2**self.precision-1) for x,y in zip(self.values, other.values)],
            precision = self.precision,
        )

    # Convert 0 to 1 or 1 to 0
    def flip(self):
        new_first_value = (self.values[0] ^ (2**(self.precision-1)))
        return Ciphertext(values=new_first_value+self.values[1:], precision=self.precision)

# Add together more than 2 ciphertexts 
def sum_ciphertexts(ciphertexts):
    assert len(ciphertexts) >= 1
    L, p = len(ciphertexts[0].values), ciphertexts[0].precision
    for c in ciphertexts[1:]:
        assert c.precision == p and len(c.values) == L
    mask = (2**p)-1
    return Ciphertext(
        values = [sum([c.values[i] for c in ciphertexts]) & mask for i in range(L)],
        precision = p
    )

# Generate public key (bunch of encryptions of zero)
@dataclass
class PublicKey():
    samples: list #[Ciphertext]

def mk_public_key(key, q):
    o = []
    for i in range(100):
        rv = random_vector(len(key)-1, q)
        canceling_value = -prod(rv, key[1:], q) % q
        o.append(Ciphertext([(noise() + canceling_value) % q] + rv, q))
    return PublicKey(samples=o)

# Encrypt a value
def encrypt(key, message, precision):
    return partial_encrypt(key, message * 2**(precision-1), precision)

# Partially encrypt a value, assuming it's already in rescaled form
# (ie. 0 for a normal 0 and q/2 for a normal 1). Used directly when
# encrypting key shares for relinerization
def partial_encrypt(key, message, precision):
    assert key[0] == 1
    for k in key:
        assert k <= 2**(ERROR_BITS+1)
    rv = random_vector(len(key)-1, precision)
    # -k[1:].v[1:]. Adding this value in as v[0] ensures that
    # v[0] + k[1:].v[1:] = k[0]*v[0] + k[1:].v[1:] = k.v equals the desired message
    canceling_value = -prod(rv, key[1:], precision)
    return Ciphertext([(noise() + canceling_value + message) & (2**precision-1)] + rv, precision)

# Partially decrypt a ciphertext, providing the output (0 or 1) plus noise (e)
# Useful mainly for debugging
def partial_decrypt(key, ciphertext):
    return prod(ciphertext.values, key, ciphertext.precision)

# Another debugging method: get magnitude of error in bits
def error_bits(key, ciphertext):
    o = partial_decrypt(key, ciphertext)
    return bit_length(min(
        o,
        2**ciphertext.precision-o,
        abs(2**(ciphertext.precision-1)-o)
    ))

# Decrypt a ciphertext, outputting the message 0 or 1
def decrypt(key, ciphertext):
    values, p = ciphertext.values, ciphertext.precision
    return 1 if (2**(p-2) < prod(values, key, p) <= 3 * 2**(p-2)) else 0

# A transit key component contains 
@dataclass
class TransitKeyComponent():
    digits: list # [[ciphertext]]

    def get_combination_for(self, index):
        if index:
            return sum_ciphertexts([self.digits[power] for power in range(len(self.digits)) if (index >> power) % 2 == 1])
        else:
            return Ciphertext(values=[0] * len(self.digits[0].values), precision=self.digits[0].precision)

@dataclass
class TransitKey():
    pairs: list # [list][TransitKeyComponent]

# Converts a key [s1, s2, s3...] into a key [first bit of s1, second bit of s1 ..., first bit of s2, second bit of s2...]
def flatten_key(s):
    return sum([[(x >> i) % 2 for i in range(ERROR_BITS+1)] for x in s], [])

# Converts a ciphertext [c1, c2, c3...] into a ct [c1, 2c1, 4c1...., c2, 2c2, 4c2....]
# The goal of the above transformation and this one is that we achieve:
# (i) s . ct = flatten_key(s) . flatten_ct(ct)
# (ii) |flatten_key(s)[i]| in {0,1} for all i
# This prevents error blowup when multiplying
def flatten_ciphertext(ct):
    return Ciphertext(
        values=sum([[(x << i) & (2**ct.precision-1) for i in range(ERROR_BITS+1)] for x in ct.values], []),
        precision=ct.precision
    )

def mk_transit_key(s, t, precision):
    # For each v=s[i], and for each product v = s[i] * s[j], encrypt v, v*2, v*4, v*8.....
    # under the key `t`.
    # This allows us to compute s[i] * b or s[i] * s[j] * b for any b with a logarithmic number
    # of additions (and hence logarithmic-sized error blowup)
    s = flatten_key(s)
    pairs = []
    for i in range(len(s)):
        pairs.append([TransitKeyComponent(digits=[partial_encrypt(t, (s[i]*s[j])<<power , precision) for power in range(precision)]) for j in range(i+1)])
    return TransitKey(pairs=pairs)

def multiply_ciphertexts(c1, c2, transit_key):
    # The idea here is that we take the equation
    #
    # dec(v1) * dec(v2) = v1.s * v2.s = (v1⁰v2).(s⁰s)
    #
    # (where ⁰ is the outer product, ie. the set of all vec1[i]*vec2[j])
    #
    # except instead of evaluating it directly (we can't because we don't have the key), we
    # use the transit key to evaluate the equation as a linear combination of s[i]*s[j],
    # giving us the decrypted output, encrypted under `t` (the target key of the transit key)
    c1, c2 = flatten_ciphertext(c1), flatten_ciphertext(c2)
    v1, v2 = c1.values, c2.values
    dim, precision, mask = len(v1), c1.precision, (2**c1.precision)-1
    return sum_ciphertexts([transit_key.pairs[max(i,j)][min(i,j)].get_combination_for(((v1[i] * v2[j]) >> (precision-1)) & mask) for i in range(dim) for j in range(dim)])

# Encode an integer into a binary representation (least significant bits first)
def binary_encode(integer, length, encoded_zero, encoded_one):
    return [encoded_one if integer & (1 << i) else encoded_zero for i in range(length)]

# Decrypt a series of ciphertexts that represent an integer in binary representation
def binary_decrypt(key, output):
    return sum([decrypt(key, o) << i for i,o in enumerate(output)])

# Logical OR
def _or(a, b, tk):
    return a + b + multiply_ciphertexts(a, b, tk)

# Logical AND
_and = multiply_ciphertexts

# Kogge-Stone adder, see https://upload.wikimedia.org/wikipedia/commons/1/1c/4_bit_Kogge_Stone_Adder_Example_new.png
def encoded_add(a, b, tk):
    return kogge_stone_propagate(
        [ai + bi for ai, bi in zip(a,b)],
        [_and(ai, bi, tk) for ai, bi in zip(a,b)],
        tk
    )

def encoded_add3(a, b, c, tk):
    return kogge_stone_propagate(
        [ai + bi + ci for ai, bi, ci in zip(a,b,c)],
        [_and(ai + bi, ai + ci, tk) + ai for ai, bi, ci in zip(a,b,c)], # (a+b)(a+c)-a = (a+b)(a+c)-a2 = a2+ac+ab+bc-a2 = ab+ac+bc
        tk
    )

def kogge_stone_propagate(p, g, tk):
    origp = p[::]
    offset = 1
    while 2**offset <= len(p):
        newg = g[::]
        newp = p[::]
        for i in range(0, len(p)-offset):
            newg[i+offset] = _or(g[i+offset], _and(p[i+offset], g[i], tk), tk)
            newp[i+offset] = _and(p[i], p[i+offset], tk)
        g, p = newg, newp
        offset *= 2
    return [origp[0]] + [origp[i] + g[i-1] for i in range(1, len(p))] + [g[-1]]

# Converts a+b+c into v+w such that a+b+c = v+w. Multiplicative depth 1.
def three_to_two(a, b, c, zero, tk):
    return (
        [ai + bi + ci for ai, bi, ci in zip(a,b,c)] + [zero],
        [zero] + [_and(ai + bi, ai + ci, tk) + ai for ai, bi, ci in zip(a,b,c)]
    )

# Add together many numbers. Use the 3->2 adder in a tree structure (ok fine it's a DAG),
# then finish off with a 3-to-1 or 2-to-1 as needed
def multi_add(values, zero, tk, bits=999999999999999):
    while len(values) > 2:
        print("Multi adding {} values".format(len(values)))
        o = []
        for i in range(0, len(values)-2, 3):
            x, y = three_to_two(values[i], values[i+1], values[i+2], zero, tk)
            o.extend([x[:bits], y[:bits]])
        o.extend(values[len(values) - len(values) % 3:])
        values = o
    return encoded_add(values[0], values[1], tk)[:bits] if len(values) == 2 else encoded_add3(values[0], values[1], values[2], tk)[:bits]

# Adjusts a ciphertext's precision, chopping off lower-order bits. This does not
# magnify the error by more than a small amount!
def adjust_ciphertext_precision(ct, new_precision):
    if new_precision < ct.precision:
        new_values = [x >> (ct.precision - new_precision) for x in ct.values]
    else:
        new_values = [x << (new_precision - ct.precision) for x in ct.values]
    return Ciphertext(values=new_values, precision=new_precision)

# A key used for the bootstrapping procedure. This involves running a decryption circuit for
# scheme key `s` (and short modulus `q`) homomorphically encrypted under key `t` (and long modulus `p`)
# The bootstrapping key provides `s` encrypted under `t` to allow this computation to take place
@dataclass
class BootstrappingKey():
    values: list # [list[ciphertext]]
    zero: Ciphertext
    one: Ciphertext
    short_precision: int
    long_precision: int

def mk_bootstrapping_key(s, t, long_precision, short_precision):
    # Basically an encryption of every bit of every s[i] under t
    values = []
    for x in s:
        bits = []
        for j in range(ERROR_BITS+1):
            bits.append(encrypt(t, (x >> j) % 2, long_precision))
        values.append(bits)
    return BootstrappingKey(
        values = values,
        zero = encrypt(t, 0, long_precision),
        one = encrypt(t, 1, long_precision),
        short_precision = short_precision,
        long_precision = long_precision
    )

def bootstrap(ct, bk, tk):
    print("Bootstrapping")
    # Start by squashing the ciphertext to a shorter precision for easier calculation
    squashed_ct = adjust_ciphertext_precision(ct, bk.short_precision)
    # The i'th bin represents bits with place value 2**i
    inner_product_bits = [[] for _ in range(bk.short_precision)]
    # Compute s.ct = s[0] * ct[0] + ... + s[k-1] * ct[k-1]
    # We do this by walking through every bit in every s[i] and adding it to
    # every bin where the corresponding bit of ct[i] is 1. This is basically the
    # same as textbook addition and multiplication, except we can't add or carry
    # (yet) because the s[i] bits are all encrypted
    for ct_value, bk_value in zip(squashed_ct.values, bk.values):
        for ct_bit in range(bk.short_precision):
            if (ct_value >> ct_bit) % 2:
                for key_bit in range(min(ERROR_BITS+1, bk.short_precision-ct_bit)):
                    inner_product_bits[ct_bit + key_bit].append(bk_value[key_bit])
    print('Packed bits')
    # To combine all the bins, we'll just keep grabbing one bit from each bin, pretend
    # that's an integer, and add up all the integers
    max_inner_product_bit_count = max(len(x) for x in inner_product_bits)
    as_integer_encodings = [[inner_product_bits[i][j] if j < len(inner_product_bits[i]) else bk.zero for i in range(bk.short_precision)] for j in range(max_inner_product_bit_count)]
    print('Adding {} integers'.format(len(as_integer_encodings)))
    # Final sum mod 2**short_precision
    total = multi_add(as_integer_encodings, bk.zero, tk, bk.short_precision)
    # 1 if the top two digits are 10 or 01, 0 if they are 00 or 11
    # This is equivalent to "1 if the value is closer to q/2, 0 if it's
    # closer to 0"
    return total[bk.short_precision-1] + total[bk.short_precision-2]
