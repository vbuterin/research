# A toy homomorphic encryption implementation, based on
# https://eprint.iacr.org/2011/344.pdf
# Does not include bootstrapping, so only supports circuits up to
# logarithmically high degree

from dataclasses import dataclass
import random

# Returns the bit length of n, eg. 6 -> 110 -> 3, 31 -> 11111 -> 5, 32 -> 100000 -> 6
def bit_length(n):
    return len(bin(n))-2

# Random length-l vector in a field with the given modulus
def random_vector(l, modulus):
    return [random.randrange(modulus) for _ in range(l)]

mk_key = random_vector

# Generates an (even) noise value within some small range
def noise():
    return 2 * random.randrange(-3, 4)

# Inner product of two vectors, that is compute sum(vec1[i] * vec2[i]).
# Allows special case for empty vectors
def prod(vec1, vec2, modulus):
    return sum([i*j for i,j in zip(vec1, vec2)]) % modulus if (vec1 and vec2) else 0

# A ciphertext representing some value m in {0,1}
# Ciphertexts are of the form (aux, out) where out = aux . key + m + 2e
# where . is the inner product and 2e is an even error term
@dataclass
class Ciphertext():
    aux: list # [int]
    out: int
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
        assert self.modulus == other.modulus and len(self.aux) == len(other.aux)
        return Ciphertext(
            aux = [(x+y) % self.modulus for x,y in zip(self.aux, other.aux)],
            out = (self.out + other.out) % self.modulus,
            modulus = self.modulus,
        )

    # Convert 0 to 1 or 1 to 0
    def flip(self):
        return Ciphertext(aux=self.aux, out=(self.out+1) % self.modulus, modulus=self.modulus)

# A dummy ciphertext to represent zero
ZERO_CT = Ciphertext(aux=[], out=0, modulus=0)

def sum_ciphertexts(ciphertexts):
    L, m = len(ciphertexts[0].aux), ciphertexts[0].modulus
    if len(ciphertexts) == 0:
        return Ciphertext(aux=[], out=0, modulus=0)
    for c in ciphertexts[1:]:
        assert c.modulus == m and len(c.aux) == L
    return Ciphertext(
        aux = [sum([c.aux[i] for c in ciphertexts]) % m for i in range(L)],
        out = sum(c.out for c in ciphertexts) % m,
        modulus = m
    )

# Encrypt a value (remember: out = aux . key + m + 2e)
def encrypt(key, message, q):
    aux = random_vector(len(key), q)
    return Ciphertext(aux, (prod(aux, key, q) + noise() + message) % q, q)

# Partially decrypt a ciphertext, providing the output (0 or 1) plus noise (2e)
# Useful mainly for debugging
# Note that if the noise turns out to be negative (true half the time), then
# you get a negative value, which appears to be a value close to the modulus
# with a flipped parity, eg. if modulus = 97, -4 appears as 93
def partial_decrypt(key, ciphertext):
    if ciphertext.modulus == 0:
        return ciphertext.out
    return (ciphertext.out - prod(ciphertext.aux, key, ciphertext.modulus)) % ciphertext.modulus

# Another debugging method: get magnitude of error
def error_digits(key, ciphertext):
    o = partial_decrypt(key, ciphertext)
    return min(len(str(o)), len(str(ciphertext.modulus-o)))

# Decrypt a ciphertext, outputting the message 0 or 1
def decrypt(key, ciphertext):
    # Special case for empty aux
    if ciphertext.modulus == 0:
        return ciphertext.out % 2
    aux, out, q = ciphertext.aux, ciphertext.out, ciphertext.modulus
    half_q = (q // 2) - (q // 2) % 2
    return ((out - prod(aux, key, q) + half_q) % q) % 2

# Expresses a number in a given base
def express_in_base(n, base, length):
    # Optimized special case for base 16
    if BASE == 16:
        o = []
        for byte in n.to_bytes(length//2+1, 'little'): 
            o.append(byte%16)
            o.append(byte//16)
        return o[:length]
    else:
        o = []
        while n > 0:
            o.append(n % base)
            n //= base
        return o + [0] * (length - len(o))

BITS_PER_DIGIT = 4
BASE = 2**BITS_PER_DIGIT

@dataclass
class TransitKeyComponent():
    digits: list # [[ciphertext]]

    def get_combination_for(self, index):
        representation = express_in_base(index, BASE, len(self.digits))
        return sum_ciphertexts([self.digits[power][representation[power]-1] for power in range(len(self.digits)) if representation[power] > 0])

@dataclass
class TransitKey():
    singles: list # [TransitKeyComponent]
    pairs: list # [list][TransitKeyComponent]

@dataclass
class ModulusChangeKey():
    items: list # [TransitKeyComponent]
    old_modulus: int
    new_modulus: int

def mk_transit_key(s, t, q):
    # For each v=s[i], and for each product v = s[i] * s[j], encrypt v, v*2, v*3... v*base, v*2*base, v*3*base...
    # under the key `t`.
    # This allows us to compute s[i] * b or s[i] * s[j] * b for any b with a logarithmic number
    # of additions (and hence logarithmic-sized error blowup)

    # Determine the highest power we need to compute up to
    power_count = 0
    while BASE**power_count < q:
        power_count += 1

    # Precompute s[i]*s[j] products
    s_products = [[s[i] * s[j] % q if j>=i else None for j in range(len(s))] for i in range(len(s))]

    return TransitKey(
        singles = [
            TransitKeyComponent(digits=[[encrypt(t, (s[i] * digit) << (BITS_PER_DIGIT * power), q) for digit in range(1, BASE)] for power in range(power_count)])
            for i in range(len(s))
        ],
        pairs = [
            [
                TransitKeyComponent(digits=[[encrypt(t, (s_products[min(i,j)][max(i,j)] * digit) << (BITS_PER_DIGIT * power), q) for digit in range(1, BASE)] for power in range(power_count)])
                for j in range(len(s))
            ]
            for i in range(len(s))
        ],
    )

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
    a1, b1 = c1.aux, c1.out
    a2, b2 = c2.aux, c2.out
    assert c1.modulus == c2.modulus and len(a1) == len(a2) == len(transit_key.singles)
    dim, modulus = len(a1), c1.modulus
    # All monomials in -(b1*a2).s - (b2*a1).s
    b1a2_b2a1_s = [transit_key.singles[i].get_combination_for((-b1 * a2[i] - b2 * a1[i]) % modulus) for i in range(dim)]
    # All monomials in s1⁰s2.a1⁰a2
    a1a2_ss = [transit_key.pairs[i][j].get_combination_for(a1[i] * a2[j] % modulus) for i in range(dim) for j in range(dim)]
    # Combined ciphertext
    comb_ciphertext = sum(b1a2_b2a1_s + a1a2_ss, ZERO_CT)
    # Add b1*b2
    comb_ciphertext.out = (comb_ciphertext.out + b1 * b2) % modulus
    return comb_ciphertext

# "Rescale" a value mod q into a value mod p, in a way that preserves evenness of noise
# This works by first halving x (note that "even" here means {-2e: -q/4 < e < q/4}), getting
# us to an error of {-q/4 < e < q/4}. Multiplying by p/q gets us to {-p/4 < e' < p/4}, with
# up to 1/2 extra error from imperfect division. Then we multiply back by two.
# Note that the goal here is not to preserve evenness of x itself, the goal is that for any x,y
# we want to preserve the evenness of x-y
#
# This method is used to convert a high-modulus ciphertext into a low-modulus ciphertext, making
# it easier to decrypt in a circuit for bootstrapping
def rescale_modulus(x, q, p):
    half = (q+1)//2
    half_x = (x*half) % q
    scaled_half_x = (half_x * p) // q
    scaled_x = (scaled_half_x * 2) % p
    return scaled_x

# A transfer key similar to the key above, except with no quadratic terms and with a modulus change
# Goal is to convert a ciphertext from being encrypted under key s (modulo q) to key t (modulo p)
def mk_modulus_change_key(s, t, q, p):
    power_count = 0
    while BASE**power_count < q:
        power_count += 1
    return ModulusChangeKey(
        items = [
            TransitKeyComponent(digits=[[encrypt(t, rescale_modulus(s[i] * digit * BASE**power, q, p), p) for digit in range(1, BASE)] for power in range(power_count)])
            for i in range(len(s))
        ],
        old_modulus=q,
        new_modulus=p,
    )

# Change a ciphertext using a given modulus change key
def change_ciphertext_modulus(ct, modulus_change_key):
    oldq, newq = modulus_change_key.old_modulus, modulus_change_key.new_modulus
    new_ct = sum([
        modulus_change_key.items[i].get_combination_for(oldq - ct.aux[i])
        for i in range(len(ct.aux))
    ], ZERO_CT)
    new_ct.out = (new_ct.out + rescale_modulus(ct.out, oldq, newq)) % newq
    return new_ct

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
        for i in range(0, len(p)-offset):
            g[i+offset] = _or(g[i+offset], _and(p[i+offset], g[i], tk), tk)
            p[i+offset] = _and(p[i], p[i+offset], tk)
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
