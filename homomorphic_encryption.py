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
    depth: int # Counter to track multiplicative depth, for debug purposes

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
            depth = max(self.depth, other.depth)
        )

    # Convert 0 to 1 or 1 to 0
    def flip(self):
        return Ciphertext(aux=self.aux, out=(self.out+1) % self.modulus, modulus=self.modulus, depth=self.depth)

# A dummy ciphertext to represent zero
ZERO_CT = Ciphertext(aux=[], out=0, modulus=0, depth=1)

# Encrypt a value (remember: out = aux . key + m + 2e)
def encrypt(key, message, q):
    aux = random_vector(len(key), q)
    return Ciphertext(aux, (prod(aux, key, q) + noise() + message) % q, q, 1)

# Partially decrypt a ciphertext, providing the output (0 or 1) plus noise (2e)
# Useful mainly for debugging
# Note that if the noise turns out to be negative (true half the time), then
# you get a negative value, which appears to be a value close to the modulus
# with a flipped parity, eg. if modulus = 97, -4 appears as 93
def partial_decrypt(key, ciphertext):
    if ciphertext.modulus == 0:
        return ciphertext.out
    return (ciphertext.out - prod(ciphertext.aux, key, ciphertext.modulus)) % ciphertext.modulus

# Decrypt a ciphertext, outputting the message 0 or 1
def decrypt(key, ciphertext):
    # Special case for empty aux
    if ciphertext.modulus == 0:
        return ciphertext.out % 2
    aux, out, q = ciphertext.aux, ciphertext.out, ciphertext.modulus
    half_q = (q // 2) - (q // 2) % 2
    return ((out - prod(aux, key, q) + half_q) % q) % 2

@dataclass
class TransitKeyComponent():
    bits: list

    def get_combination_for(self, index):
        return sum([self.bits[i] for i in range(len(self.bits)) if (index & (1 << i))], ZERO_CT)

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
    # For each v=s[i], and for each product v=s[i] * s[j], encrypt v, v*2, v*4... under the key `t`
    # This allows us to compute s[i] * b or s[i] * s[j] * b for any b with a logarithmic number
    # of additions (and hence logarithmic-sized error blowup)
    bit_count = bit_length(q)
    return TransitKey(
        singles = [
            TransitKeyComponent(bits=[encrypt(t, s[i] * 2**b, q) for b in range(bit_count)])
            for i in range(len(s))
        ],
        pairs = [
            [
                TransitKeyComponent(bits=[encrypt(t, s[i] * s[j] * 2**b, q) for b in range(bit_count)])
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
    comb_ciphertext.depth = max(c1.depth, c2.depth) + 1
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
    bit_count = bit_length(q)
    return ModulusChangeKey(
        items = [
            TransitKeyComponent(bits=[encrypt(t, rescale_modulus(s[i] * 2**b, q, p), p) for b in range(bit_count)])
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
    origp = [ai + bi for ai, bi in zip(a,b)]
    p = origp[::]
    g = [_and(ai, bi, tk) for ai, bi in zip(a,b)]
    offset = 1
    while 2**offset <= len(a):
        for i in range(0, len(a)-offset):
            g[i+offset] = _or(g[i+offset], _and(p[i+offset], g[i], tk), tk)
            p[i+offset] = _and(p[i], p[i+offset], tk)
        offset *= 2
    return [origp[0]] + [origp[i] + g[i-1] for i in range(1, len(a))] + [g[-1]]

# Lower-depth adder for k>2 values. For log(p) bit values, has depth log(k)+log(log(k))*log(p)
def multi_add(values, zero, tk):
    # Step 1: split out each digit, eg. for the ith digit v1[i], v2[i] ... vk[i]
    # and compute the log(k) bit sized sum(v1[i] + ...  vk[i])
    # Depth: log(k)
    bit_collections = [[[v[i]] for v in values] for i in range(len(values[0]))]
    for i in range(len(bit_collections)):
        print("Processing the %d'th collection" % i)
        b = bit_collections[i]
        while len(b) > 1:
            print(len(b))
            if len(b) % 2:
                b.append([zero] * len(b[-1]))
            b = [encoded_add(b[j], b[j+1], tk) for j in range(0, len(b), 2)]
        # Output the sums
        bit_collections[i] = b[0]
    # Perform the rearranging procedure again, converting log(p) log(k)-bit numbers with offsets
    # into log(k) log(p)-bit numbers
    summands = [[zero]*i+[b[i] for b in bit_collections]+[zero]*(len(bit_collections[0])-1-i) for i in range(len(bit_collections[0]))]
    # Add together the log(k) values of size p each. Complexity log(log(k))*log(p)
    print("Combining outputs")
    while len(summands) > 1:
        print(len(summands))
        if len(summands) % 2:
            summands.append([zero] * len(summands[-1]))
        summands = [encoded_add(summands[j], summands[j+1], tk) for j in range(0, len(summands), 2)]
    return summands[0]

def next_prime(n):
    x = n + (n % 2) + 1
    while pow(2, x, x) != 2:
        x += 2
    return x

short_modulus = 997
medium_modulus = 99991
huge_modulus = next_prime(10**150)

def test_add(x, y):
    S = random_vector(5, huge_modulus)
    zero = encrypt(S, 0, huge_modulus)
    one = encrypt(S, 1, huge_modulus)
    encx = binary_encode(x, 8, zero, one)
    ency = binary_encode(y, 8, zero, one)
    encz = encoded_add(encx, ency, mk_transit_key(S, S, huge_modulus))
    assert sum([decrypt(S, x) << i for i, x in enumerate(encz[:8])]) == x + y

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
    u = random_vector(3, short_modulus)
    uk = mk_modulus_change_key(s, u, medium_modulus, short_modulus)
    for c in (c1, c2, c3, c4):
        assert decrypt(u, change_ciphertext_modulus(c, uk)) == decrypt(s, c)
    for cl in (c1, c2, c3, c4):
        for cr in (c1, c2, c3, c4):
            cm = multiply_ciphertexts(cl, cr, tk)
            tuk = mk_modulus_change_key(t, u, medium_modulus, short_modulus)
            assert decrypt(u, change_ciphertext_modulus(cm, tuk)) == decrypt(t, cm) == decrypt(s, cl) * decrypt(s, cr)
    S = random_vector(5, huge_modulus)
    zero = encrypt(S, 0, huge_modulus)
    one = encrypt(S, 1, huge_modulus)
    forty_two = binary_encode(42, 8, zero, one)
    sixty_nine = binary_encode(69, 8, zero, one)
    one_one_one = encoded_add(forty_two, sixty_nine, mk_transit_key(S, S, huge_modulus))
    assert binary_decrypt(S, one_one_one) == 111
    two, three, four, five = (binary_encode(i, 4, zero, one) for i in (2, 3, 4, 5))
    fourteen = multi_add([two, three, four, five], zero, mk_transit_key(S, S, huge_modulus))
    assert binary_decrypt(S, fourteen[:4]) == 14
    print("Tests passed")
    
if __name__ == '__main__':
    test()
