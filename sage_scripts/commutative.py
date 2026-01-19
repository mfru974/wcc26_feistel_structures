from sboxUv2 import Sb,feistel_round,swap_halves,algebraic_degree,oplus,rand_Sbox

from sage.crypto.sboxes import sboxes


# ============================================================
# Definition of the commuting function G
# ============================================================

def G(a, f):
    """
    Computes the function G such that for all x, y:
        G(x || y) = (x + D_a f(y) || y + a).

    :param a: element of F_2^m
    :param f: Sbox from F_2^m to F_2^n
    """
    lut = []
    N = f.input_space_size()
    der = f.derivative(a)

    for y in f.input_space():
        for x in f.output_space():
            lut.append(N * oplus(y, a) + oplus(x, der[y]))

    return Sb(lut)


# ============================================================
# Parameters
# ============================================================

n = 4
m = 4
nb_tries = 100


# ============================================================
# Verification of Scream and iScream decompositions
# ============================================================

print("Checking the Feistel decompositions of Scream and iScream:\n")

SW = swap_halves(8)

f1 = Sb([0, 2, 0, 0xB, 3, 0, 0, 0xA, 1, 0xE, 0, 6, 0xA, 4, 5, 2])
f2 = Sb([0, 2, 0xC, 7, 5, 0xF, 0xD, 6, 4, 0xE, 8, 9, 3, 1, 0xB, 0xA])
f3 = Sb([2, 0, 0xB, 0, 0, 3, 0xA, 0, 0xE, 1, 6, 0, 4, 0xA, 2, 5])

print("Constructing Scream's S-box as a 3-round Feistel structure:")
S_scream = feistel_round(f3) * feistel_round(f2) * feistel_round(f1) * SW
print("Succeeded\n" if S_scream == Sb(sboxes["Scream"]) else "Failed\n")

print("Constructing iScream's S-box as a 3-round Feistel structure:")
f4 = Sb([0, 8, 6, 13, 5, 15, 7, 12, 4, 14, 2, 3, 9, 1, 11, 10])
S_iscream = SW * (feistel_round(f4) ** 3)
print("Succeeded\n" if S_iscream == Sb(sboxes["iScream"]) else "Failed\n")


# ============================================================
# Experimental verification of Proposition 1
# ============================================================

print("-" * 50)
print(
    f"Testing Proposition 1 ({nb_tries} random instances) "
    f"with n={n}, m={m}:\n"
)

test = True

for _ in range(nb_tries):
    g1 = rand_Sbox(m, n)
    g2 = rand_Sbox(n, m)
    g3 = rand_Sbox(m, n)

    S = swap_halves(8) * feistel_round(g3) * feistel_round(g2) * feistel_round(g1)

    for a in range(1, 16):
        if S * G(a, g1) != G(a, g3) * S:
            test = False
            break

print("Successful experiment\n" if test else "Experiment failed\n")


# ============================================================
# Analysis of derivatives
# ============================================================

print("-" * 50)
print("Checking the derivatives of f1 and f3:\n")

for a in range(1, 16):
    if f1.derivative(a) == f3.derivative(a):
        print(
            f"D_{a} f1 = D_{a} f3 "
            f"(degree {algebraic_degree(f1.derivative(a))})"
        )

print("-" * 50)
print("Checking the derivatives of f4:\n")

degrees_f4 = [algebraic_degree(f4.derivative(a)) for a in range(1, 16)]
print("Degrees of the 15 non-trivial derivatives of f4:")
print(degrees_f4)
