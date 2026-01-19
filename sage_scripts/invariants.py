from sboxUv2 import (
    Sb,
    oplus,
    scal_prod,
    swap_halves,
    algebraic_degree,
    algebraic_normal_form_coordinate,
)
from sage.all import (
    GF,
    VectorSpace,
    span,
    vector,
    Permutation,
    log,
)
from sage.crypto.boolean_function import BooleanFunction
from sage.crypto.sboxes import sboxes


# ============================================================
# Auxiliary functions (from TLS19)
# ============================================================

def hamming_weight(v):
    """Returns the Hamming weight of a vector."""
    return sum(1 for x in v if x != 0)


def is_smaller_equal(a, v):
    """Checks whether vector a is component-wise smaller or equal to v."""
    return all(x <= y for x, y in zip(a, v))


def basis_invariants(S):
    """
    Computes a basis of invariant Boolean functions
    associated to the permutation S (TLS19).
    """
    all_even = True
    R = []

    cycles = Permutation([x + 1 for x in S]).to_cycles()

    for cycle in cycles:
        if len(cycle) % 2 == 1:
            all_even = False

    B = [0] * len(S)
    for cycle in cycles:
        for l in cycle:
            B[l - 1] = 1
        R.append(BooleanFunction(B))
        B = [0] * len(S)

    if all_even:
        B = [0] * len(S)
        for cycle in cycles:
            value = 0
            for i in range(len(cycle)):
                B[cycle[i] - 1] = value
                value = (value + 1) % 2
        R.append(BooleanFunction(B))

    return R


def all_from_basis(B_S):
    """Generates all invariants from a basis."""
    L = []
    bitlength = B_S[0].nvariables()
    V = VectorSpace(GF(2), len(B_S))

    for v in V:
        B = BooleanFunction([0] * (2 ** bitlength))
        for i in range(len(B_S)):
            if v[i] == 1:
                B += B_S[i]
        L.append(B)

    return L


def coefficient_vector_ANF(b):
    """Returns the ANF coefficient vector of a Boolean function."""
    V = VectorSpace(GF(2), b.nvariables())
    n = V.dimension()

    r = vector(GF(2), [0] * (2 ** n))
    vectors = list(V)
    b_values = [b(list(v)) for v in vectors]

    for i in range(len(vectors)):
        value = 0
        for j in range(len(vectors)):
            if is_smaller_equal(vectors[j], vectors[i]):
                value += b_values[j]
        r[i] = value % 2

    return r


def basis_from_anf_space(V):
    """Converts a vector space of ANF coefficient vectors into Boolean functions."""
    BV = V.basis()
    bitlength = int(log(len(BV[0]), 2))

    W = VectorSpace(GF(2), bitlength)
    W_vectors = list(W)

    R = []
    for v in BV:
        values = list(v)
        b_values = []

        for k in range(len(W_vectors)):
            value = sum(
                values[i]
                for i in range(len(values))
                if is_smaller_equal(W_vectors[i], W_vectors[k])
            )
            b_values.append(value % 2)

        R.append(b_values)

    return R


def all_invariants_up_to_degree(S, d):
    """Computes all invariants of degree at most d (TLS19)."""
    B_S = basis_invariants(S)
    bitlength = B_S[0].nvariables()

    B = [coefficient_vector_ANF(b) for b in B_S]
    VS = span(B, GF(2))

    W = VectorSpace(GF(2), bitlength)
    C = [
        vector([0] * i + [1] + [0] * ((2 ** bitlength) - 1 - i))
        for i in range(2 ** bitlength)
        if hamming_weight(W[i]) <= d
    ]

    V = span(C, GF(2))
    return basis_from_anf_space(V.intersection(VS))


# ============================================================
# Lifting invariants to Feistel structures
# ============================================================

def from_small_to_big_invariant(u, f):
    """
    Generates the Boolean function:
        g(x || y) = u · f(y) + u · x
    """
    lut = []
    for y in f.input_space():
        for x in f.output_space():
            lut.append(oplus(scal_prod(u, f[y]), scal_prod(u, x)))
    return Sb(lut)


# ============================================================
# Experiments
# ============================================================

print("Quadratic invariants of Scream's S-box\n")
print("Using the TLS19 algorithm, we obtain:\n")

S1 = list(sboxes["Scream"])
for g in all_invariants_up_to_degree(S1, 2):
    print(BooleanFunction(g).algebraic_normal_form())


f1 = Sb([0,2,0,0xB,3,0,0,0xA,1,0xE,0,6,0xA,4,5,2])
f3 = Sb([2,0,0xB,0,0,3,0xA,0,0xE,1,6,0,4,0xA,2,5])

common_components = [
    u for u in range(1, 16)
    if len(set((f1.component(u) + f3.component(u)).lut())) == 1
]

print(
    "\nValues of u such that u · f1 + u · f3 is constant:",
    common_components
)

print(
    "\nCorresponding invariants from Proposition 2:",
    [
        algebraic_normal_form_coordinate(
            from_small_to_big_invariant(u, f1) * swap_halves(8)
        )
        for u in common_components
    ],
    "\n"
)


print("Quadratic invariants of iScream's S-box\n")
print("Using the TLS19 algorithm, we obtain:\n")

S2 = list(sboxes["iScream"])
for g in all_invariants_up_to_degree(S2, 2):
    print(BooleanFunction(g).algebraic_normal_form())


f4 = Sb([0, 8, 6, 13, 5, 15, 7, 12, 4, 14, 2, 3, 9, 1, 11, 10])

quad_components = [
    u for u in range(1, 16)
    if algebraic_degree(f4.component(u)) <= 2
]

print(
    "\nValues of u such that u · f4 is quadratic:",
    quad_components
)

print(
    "\nCorresponding invariants:",
    [
        algebraic_normal_form_coordinate(
            from_small_to_big_invariant(u, f4)
        )
        for u in quad_components
    ],
    "\n"
)

print("Quadratic invariants of Skinny64's S-box\n")
print("Using the TLS19 algorithm, we obtain:\n")

S3 = list(sboxes["SKINNY_4"])
for g in all_invariants_up_to_degree(S3, 2):
    print(BooleanFunction(g).algebraic_normal_form())
