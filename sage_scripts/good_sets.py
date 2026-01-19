from math import log
from sboxUv2 import (
    Sb, ddt, xddt, yddt, zddt,
    feistel_round, swap_halves,
    BinLinearBasis, oplus, F2_trans,
    algebraic_degree, is_affine,
    differential_uniformity,
)


# ============================================================
# Candidate sets given by Theorem 1
# ============================================================

def candidate_good_sets_feistel(f1, f3):
    """Returns candidate good sets of differences given by Theorem 1."""
    if not (
        f1.input_space_size() == f1.output_space_size()
        == f3.input_space_size() == f3.output_space_size()
    ):
        raise ValueError("Incompatible Sbox dimensions")

    N = f1.input_space_size()
    DDT_f1, DDT_f3 = ddt(f1), ddt(f3)
    res = []

    for b in range(1, N):
        if (
            set(DDT_f1[b]) == {0, 2}
            and set(DDT_f3[b]) == {0, 2}
            and set(f1.derivative(b)) == set(f3.derivative(b))
        ):
            res.append([x + N * b for x in set(f1.derivative(b))])
    return res


def are_translated(f, g):
    """Returns True if f(x + c) = g(x) for some c."""
    for c in f.input_space():
        if f * F2_trans(c, bit_length=f.get_input_length()) == g:
            return True
    return False


def candidate_perfect_sets_feistel(f1, f3):
    """Returns candidate perfect sets of differences given by Theorem 1."""
    N = f1.input_space_size()
    DDT_f1, DDT_f3 = ddt(f1), ddt(f3)
    res = []

    for b in range(1, N):
        if (
            set(DDT_f1[b]) == {0, 2}
            and set(DDT_f3[b]) == {0, 2}
            and set(f1.derivative(b)) == set(f3.derivative(b))
            and algebraic_degree(f1.derivative(b)) <= 1
            and are_translated(f1.derivative(b), f3.derivative(b))
        ):
            res.append([x + N * b for x in set(f1.derivative(b))])
    return res


# ============================================================
# Verification tools (Definitions 3 and 4)
# ============================================================

def get_basis(myset):
    """Returns a linear basis of myset if it is a vector space."""
    s = int(log(len(myset), 2))
    if 2**s != len(myset):
        return None
    basis = list(BinLinearBasis(myset))
    return basis if len(basis) == s else None


def get_E(a, outputs, XDDT):
    """Computes the space E appearing in Definition 3."""
    temp = []
    for b in outputs:
        temp += XDDT[a][b]
    return get_basis([oplus(temp[0], x) for x in temp])


def test_V(A, YDDT, E):
    """Checks Condition (1) of Definition 3."""
    base_E = BinLinearBasis(E)
    for a in A:
        for b in A:
            base_Y = get_basis([oplus(YDDT[a][b][0], x) for x in YDDT[a][b]])
            if base_Y is None:
                return False
            for elt in base_Y:
                if not base_E.is_in_span(elt):
                    return False
    return True


def sum_condition(A, XDDT, E):
    """
    Checks Condition (1) of Definition 4.
    """
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            for k in range(j + 1, len(A)):
                lhs = oplus(
                    oplus(XDDT[A[i]][A[i]][0], XDDT[A[j]][A[j]][0]),
                    oplus(
                        XDDT[A[k]][A[k]][0],
                        XDDT[oplus(A[i], oplus(A[j], A[k]))]
                            [oplus(A[i], oplus(A[j], A[k]))][0],
                    ),
                )
                if lhs not in E:
                    return False
    return True


def xy_condition(A, XDDT, YDDT, E):
    """
    Checks Condition (2) of Definition 4.
    """
    coset_shift = oplus(
        XDDT[A[0]][A[0]][0],
        YDDT[A[0]][A[0]][0],
    )
    coset = [oplus(t, coset_shift) for t in E]

    for a in A:
        for b in A:
            if oplus(XDDT[a][b][0], YDDT[b][a][0]) not in coset:
                return False
    return True



def is_good_set(A, S):
    """Checks whether A is a good set of differences for S."""
    DDT = ddt(S)
    XDDT, YDDT, ZDDT = xddt(S), yddt(S), zddt(S)
    w = DDT[A[0]][A[0]]

    for a in A:
        for b in A:
            if DDT[a][b] != w or not is_affine(ZDDT[a][b]):
                return False

    E = get_E(A[0], A, XDDT)
    return E is not None and test_V(A, YDDT, E)


def is_perfect_set(A, S):
    """Checks whether A is a perfect set of differences for S."""
    if not is_good_set(A, S):
        return False

    XDDT, YDDT = xddt(S), yddt(S)
    temp = []
    for b in A:
        temp += XDDT[A[0]][b]
    E = [oplus(temp[0], x) for x in temp]

    return (
        sum_condition(A, XDDT, E)
        and xy_condition(A, XDDT, YDDT, E)
    )


# ============================================================
# Experiments
# ============================================================

def main():
    f1 = Sb([0,2,0,0xb,3,0,0,0xa,1,0xe,0,6,0xa,4,5,2])
    f2 = Sb([0,2,0xc,7,5,0xf,0xd,6,4,0xe,8,9,3,1,0xb,0xa])
    f3 = Sb([2,0,0xb,0,0,3,0xa,0,0xe,1,6,0,4,0xa,2,5])

    Feistel = swap_halves(8) * feistel_round(f3) * feistel_round(f2) * feistel_round(f1)

    good_sets = candidate_good_sets_feistel(f1, f3)
    print(f"Theorem 1 gives {len(good_sets)} good sets for Scream's Sbox.")
    print("Verification:", all(is_good_set(A, Feistel) for A in good_sets))

    perfect_sets = candidate_perfect_sets_feistel(f1, f3)
    print(f"Theorem 1 gives {len(perfect_sets)} perfect sets.")
    print("Verification:", all(is_perfect_set(A, Feistel) for A in perfect_sets))

    print("\nExamining the round functions f1 and f3:\n")
    print(
        f"Both are APN: the differential uniformity of f1 is "
        f"{differential_uniformity(f1)}, and that of f3 is "
        f"{differential_uniformity(f3)}.\n"
    )

    same_image_derivatives = []
    same_affine_image_derivatives = []
    same_affine_derivatives = []

    for a in range(1, 16):
        Df1 = f1.derivative(a)
        Df3 = f3.derivative(a)

        if set(Df1) == set(Df3):
            same_image_derivatives.append(a)

            if is_affine(list(set(Df1))):
                same_affine_image_derivatives.append(a)

                if algebraic_degree(Df1) <= 1:
                    same_affine_derivatives.append(a)

    print(
        "Values of b such that Im(D_b f1) = Im(D_b f3):",
        same_image_derivatives, "\n"
    )
    print(
        "Values of b such that this image is affine:",
        same_affine_image_derivatives, "\n"
    )
    print(
        "Values of b such that D_b f1 is affine:",
        same_affine_derivatives, "\n"
    )

    print(
        "Verification that D_0x1 f1 = D_0x1 f3:",
        f1.derivative(1) == f3.derivative(1), "\n"
    )


    print("\nNow examining the function f4:\n")

    f4 = Sb([0, 8, 6, 13, 5, 15, 7, 12, 4, 14, 2, 3, 9, 1, 11, 10])

    affine_derivatives = []
    affine_image_derivatives = []

    for a in range(1, 16):
        Df4 = f4.derivative(a)

        if algebraic_degree(Df4) <= 1:
            affine_derivatives.append(a)

        if is_affine(list(set(Df4))):
            affine_image_derivatives.append(a)

    print(
        "Values of b such that D_b f4 is affine:",
        affine_derivatives, "\n"
    )
    print(
        "Values of b such that Im(D_b f4) is affine:",
        affine_image_derivatives, "\n"
    )



if __name__ == "__main__":
    main()


    







        


