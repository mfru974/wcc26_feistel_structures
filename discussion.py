from sboxUv2 import Sb, ddt, algebraic_degree, F2_trans


# ============================================================
# Testing sufficient conditions for each cryptanalytic property
# ============================================================

def has_affine_commutant(f1, f2, f3):
    """Checks Proposition 1: existence of an affine commuting function."""
    if f1.input_space_size() != f3.input_space_size():
        raise ValueError("Incompatible Sbox dimensions")

    for a in range(1, f1.input_space_size()):
        da = f1.derivative(a)
        if da == f3.derivative(a) and algebraic_degree(da) <= 1:
            return True
    return False


def has_quadratic_invariant(f1, f2, f3):
    """Checks Proposition 2: existence of a quadratic invariant."""
    if f1.output_space_size() != f3.output_space_size():
        raise ValueError("Incompatible Sbox dimensions")

    for u in range(1, f1.output_space_size()):
        if len(set(f1.component(u) + f3.component(u))) == 1:
            return True
    return False


def has_perfect_sets(f1, f2, f3):
    """Checks Theorem 1: existence of perfect sets."""
    if not f2.is_invertible():
        return False

    n = f1.input_space_size()
    DDT_f1, DDT_f3 = ddt(f1), ddt(f3)

    for b in range(1, n):
        if (
            set(DDT_f1[b]) == {0, 2}
            and set(DDT_f3[b]) == {0, 2}
            and f1.derivative(b) == f3.derivative(b)
            and algebraic_degree(f1.derivative(b)) <= 1
        ):
            for t in range(n):
                if f1.derivative(b) == f3.derivative(b) * F2_trans(
                    t, bit_length=f1.get_input_length()
                ):
                    return True
    return False


# ============================================================
# Declaring the round functions
# ============================================================

f1 = Sb([0,2,0,0xb,3,0,0,0xa,1,0xe,0,6,0xa,4,5,2])
f2 = Sb([0,2,0xc,7,5,0xf,0xd,6,4,0xe,8,9,3,1,0xb,0xa])
f3 = Sb([2,0,0xb,0,0,3,0xa,0,0xe,1,6,0,4,0xa,2,5])
f4 = Sb([0,8,6,0xd,5,0xf,7,0xc,4,0xe,2,3,9,1,0xb,0xa])
f5 = Sb([0xc,5,0xa,2,5,5,7,6,5,0xb,0,0xf,4,3,5,3])
f6 = Sb([3,1,9,0,6,7,8,2,9,6,0xd,9,9,5,9,0xe])


def main():
    print("Verifying the results of Table 1:\n")
    tests = [(f1,f2,f3),(f1,f1,f3),(f4,f4,f4),(f5,f2,f6),(f5,f1,f6),(f1,f4,f6)]

    for f, g, h in tests:
        print(
            has_affine_commutant(f, g, h),
            has_quadratic_invariant(f, g, h),
            has_perfect_sets(f, g, h),
        )


if __name__ == "__main__":
    main()
