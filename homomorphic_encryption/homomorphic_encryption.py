# A toy homomorphic encryption implementation, based on
# https://eprint.iacr.org/2011/344.pdf
# Does not include bootstrapping, so only supports circuits up to
# logarithmically high degree

from dataclasses import dataclass
import random

ERROR_DIGITS = 4

# Returns the bit length of n, eg. 6 -> 110 -> 3, 31 -> 11111 -> 5, 32 -> 100000 -> 6
def bit_length(n):
    return len(bin(n))-2

# Random length-l vector in a field with the given modulus
def random_vector(l, modulus):
    return [random.randrange(modulus) for _ in range(l)]

# Generates an (even) noise value within some small range
def noise():
    return random.randrange(-2**ERROR_DIGITS, ERROR_DIGITS)

# Inner product of two vectors, that is compute sum(vec1[i] * vec2[i]).
# Allows special case for empty vectors
def prod(vec1, vec2, modulus):
    return sum([i*j for i,j in zip(vec1, vec2)]) % modulus if (vec1 and vec2) else 0

# Generates a private key
def generate_key(length, modulus):
    return [1] + [(2**ERROR_DIGITS + noise()) % modulus for _ in range(1, length)]

# A ciphertext representing some value m in {0,1}
# Ciphertexts are of the form (aux, out) where out = aux . key + m + 2e
# where . is the inner product and 2e is an even error term
@dataclass
class Ciphertext():
    values: list # [int]
    modulus: int

    # Add together two ciphertexts into one, linearly combining the aux and out.
    # Note that this combines the magnitudes of the error, so if you add waaaaay
    # too many times the error may overflow
    def __add__(self, other):
        # Special cases for "null" ciphertexts
        if self.modulus == 0:
            return other
        if other.modulus == 0:
            return self
        assert self.modulus == other.modulus and len(self.values) == len(other.values)
        return Ciphertext(
            values = [(x+y) % self.modulus for x,y in zip(self.values, other.values)],
            modulus = self.modulus,
        )

    # Convert 0 to 1 or 1 to 0
    def flip(self):
        new_first_value = (self.values[0] + self.modulus // 2) % self.modulus
        return Ciphertext(values=new_first_value+self.values[1:], modulus=self.modulus)

# A dummy ciphertext to represent zero
ZERO_CT = Ciphertext(values=[], modulus=0)

def sum_ciphertexts(ciphertexts):
    assert len(ciphertexts) >= 1
    L, m = len(ciphertexts[0].values), ciphertexts[0].modulus
    for c in ciphertexts[1:]:
        assert c.modulus == m and len(c.values) == L
    return Ciphertext(
        values = [sum([c.values[i] for c in ciphertexts]) % m for i in range(L)],
        modulus = m
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

# Encrypt a value (remember: out = aux . key + m + 2e)
def encrypt(key, message, q):
    return partial_encrypt(key, message * q // 2, q)

# Partially encrypt a value, assuming it's already in rescaled form
# (ie. 0 for a normal 0 and q//2 for a normal 1). Used for encrypting
# key shares for relinerization
def partial_encrypt(key, message, q):
    assert key[0] == 1
    for k in key:
        assert k <= 2**(ERROR_DIGITS+1)
    rv = random_vector(len(key)-1, q)
    canceling_value = -prod(rv, key[1:], q) % q
    return Ciphertext([(noise() + canceling_value + message) % q] + rv, q)

# Partially decrypt a ciphertext, providing the output (0 or 1) plus noise (2e)
# Useful mainly for debugging
# Note that if the noise turns out to be negative (true half the time), then
# you get a negative value, which appears to be a value close to the modulus
# with a flipped parity, eg. if modulus = 97, -4 appears as 93
def partial_decrypt(key, ciphertext):
    return prod(ciphertext.values, key, ciphertext.modulus)

# Another debugging method: get magnitude of error
def error_digits(key, ciphertext):
    o = partial_decrypt(key, ciphertext)
    return len(str(min(
        o,
        ciphertext.modulus-o,
        abs(ciphertext.modulus//2-o)
    )))

# Decrypt a ciphertext, outputting the message 0 or 1
def decrypt(key, ciphertext):
    values, q = ciphertext.values, ciphertext.modulus
    return 1 if (q//4 < prod(values, key, q) <= (q*3)//4) else 0

@dataclass
class TransitKeyComponent():
    digits: list # [[ciphertext]]

    def get_combination_for(self, index):
        if index:
            return sum_ciphertexts([self.digits[power] for power in range(len(self.digits)) if (index >> power) % 2 == 1])
        else:
            return Ciphertext(values=[0] * len(self.digits[0].values), modulus=self.digits[0].modulus)

@dataclass
class TransitKey():
    pairs: list # [list][TransitKeyComponent]

@dataclass
class ModulusChangeKey():
    items: list # [TransitKeyComponent]
    old_modulus: int
    new_modulus: int

def flatten_key(s, modulus):
    return sum([[(x >> i) % 2 for i in range(ERROR_DIGITS+1)] for x in s], [])

def flatten_ciphertext(ct):
    return Ciphertext(
        values=sum([[(x << i) % ct.modulus for i in range(ERROR_DIGITS+1)] for x in ct.values], []),
        modulus=ct.modulus
    )

def mk_transit_key(s, t, q):
    # For each v=s[i], and for each product v = s[i] * s[j], encrypt v, v*2, v*3... v*base, v*2*base, v*3*base...
    # under the key `t`.
    # This allows us to compute s[i] * b or s[i] * s[j] * b for any b with a logarithmic number
    # of additions (and hence logarithmic-sized error blowup)
    s = flatten_key(s, q)
    pairs = []
    for i in range(len(s)):
        pairs.append([TransitKeyComponent(digits=[partial_encrypt(t, ((s[i]*s[j])<<power) % q, q) for power in range(bit_length(q))]) for j in range(i+1)])
    return TransitKey(pairs=pairs)

def multiply_ciphertexts(c1, c2, transit_key):
    # The idea here is that we take the equation
    #
    # dec(c1) * dec(c2) = (b1 - s.a1) * (b2 - s.a2) = b1b2 - (b1*a2).s - (b2*a1).s + s1⁰s2.a1⁰a2
    #
    # (where ⁰ is the outer product, ie. the set of all vec1[i]*vec2[j])
    #
    # except instead of evaluating it directly (we can't because we don't have the key), we
    # use the transit key to evaluate the equation as a linear combination of s[i] and s[i]*s[j],
    # giving us the decrypted output, encrypted under `t` (the target key of the transit key)
    c1, c2 = flatten_ciphertext(c1), flatten_ciphertext(c2)
    v1, v2 = c1.values, c2.values
    dim, modulus = len(v1), c1.modulus
    return sum_ciphertexts([transit_key.pairs[max(i,j)][min(i,j)].get_combination_for((v1[i] * v2[j] * 2 // modulus) % modulus) for i in range(dim) for j in range(dim)])

def binary_encode(integer, length, encoded_zero, encoded_one):
    return [encoded_one if integer & (1 << i) else encoded_zero for i in range(length)]

def binary_decrypt(key, output):
    return sum([decrypt(key, o) << i for i,o in enumerate(output)])

def _or(a, b, tk):
    return a + b + multiply_ciphertexts(a, b, tk)

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
def multi_add(values, zero, tk):
    while len(values) > 2:
        o = []
        for i in range(0, len(values)-2, 3):
            x, y = three_to_two(values[i], values[i+1], values[i+2], zero, tk)
            o.extend([x, y])
        o.extend(values[len(values) - len(values) % 3:])
        values = o
    return encoded_add(values[0], values[1], tk) if len(values) == 2 else encoded_add3(values[0], values[1], values[2], tk)

def next_prime(n):
    x = n + (n % 2) + 1
    while pow(2, x, x) != 2:
        x += 2
    return x

# Takes the sum of the values and returns two bits:
# (i) is sum(values) in (range_start....range_end-1)
# (ii) is sum(values) - range_start odd (1) or even (0)?
def add_and_return_remainder_parity_and_rangecheck(values, range_start, range_end, zero, one, tk):
    assert range_end > range_start >= 1
    bit_count = bit_length(range_end)
    # 10000000000 - range_start = 1111111111 - (range_start-1)
    encoded_complement_rstart = binary_encode(range_start - 1, bit_count, one, zero)
    # 10000000000 - range_end = 1111111111 - (range_end-1)
    encoded_complement_rend = binary_encode(range_end - 1, bit_count, one, zero)
    # Should be >= 10000000000 if in range
    sum_with_range_start = multi_add(values + [encoded_complement_rstart], zero, tk)
    # Should be < 10000000000 if in range
    sum_with_range_end = multi_add(values + [encoded_complement_rend], zero, tk)
    is_in_range = _and(sum_with_range_start[bit_count], one + sum_with_range_end[bit_count], tk)
    parity = sum_with_range_start[0]
    return is_in_range, parity

# A key used for the bootstrapping procedure. This involves running a decryption circuit for
# scheme key `s` (and short modulus `q`) homomorphically encrypted under key `t` (and long modulus `p`)
# The bootstrapping key provides `s` encrypted under `t` to allow this computation to take place
@dataclass
class BootstrappingKeyComponent():
    factors: list # [encoded ciphertext]

@dataclass
class BootstrappingKey():
    values: list # [BootstrappingKeyComponent]
    zero: Ciphertext
    one: Ciphertext
    old_modulus: int
    new_modulus: int

def mk_bootstrapping_key(s, t, q, p):
    bit_count = bit_length(q)
    # To make bootstrapping even easier, we generate *all* possible multiples of each s[i]
    # This means that evaluating an inner product s[0]*aux[0] + ... + s[k-1]*aux[k-1] is just a sum
    values = []
    for i in range(len(s)):
        factors = []
        for j in range(q):
            factors.append(binary_encode((s[i] * j)%q, bit_count, encrypt(t, 0, p), encrypt(t, 1, p)))
        values.append(BootstrappingKeyComponent(factors=factors))
    return BootstrappingKey(
        values = values,
        zero = encrypt(t, 0, p),
        one = encrypt(t, 1, p),
        old_modulus = q,
        new_modulus = p
    )

def bootstrap(ct, bk, tk):
    print("Bootstrapping")
    # Components of -(s[0] * aux[0] + ... + s[k-1] * aux[k-1])
    inner_product_components = [bk.values[i].factors[bk.old_modulus-ct.aux[i]] for i in range(len(ct.aux))]
    # Components of q/2 + o - (s[0] * aux[0] + ... + s[k-1] * aux[k-1])
    # For q/2, we take the highest even number below q/2, to avoid changing the parity of the result
    # The goal of the q/2 offset is to move the range -q/2...q/2 into the range 0.....q, so we do not
    # have to deal with negative numbers
    all_components = [binary_encode(ct.out + (bk.old_modulus // 4) * 2, len(inner_product_components[0]), bk.zero, bk.one)] + inner_product_components
    # For a series of ranges, we compute (1 if start <= sum(all_components) < end else 0) and (sum(all_components) - start) % 2
    inrange_and_parity = []
    for i in range(len(bk.values)+1):
        print("Computing inrange and parity bit {} of {}".format(i, len(bk.values)+1))
        inrange_and_parity.append(add_and_return_remainder_parity_and_rangecheck(all_components, (bk.old_modulus * i) or 2, bk.old_modulus * (i+1), bk.zero, bk.one, tk))
    # Pick out the correct range (this is basically an inefficient but depth-minimizing way of taking the
    # sum mod q) and return the parity within that range
    print("Inrange and parity bits computed, computing inner product...")
    return sum([multiply_ciphertexts(r, p, tk) for r,p in inrange_and_parity], ZERO_CT)
