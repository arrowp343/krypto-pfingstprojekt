from __future__ import annotations

from time import perf_counter

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    print(
        'Please install the Python library "matplotlib", otherwise Task 2c will not work.'
    )
    print("https://matplotlib.org/stable/users/getting_started/index.html")
    input("Press enter to ignore this for now > ")

########## DO NOT TOUCH ##########


class Curve(object):
    def __init__(self, p: int, a: int, b: int, x: int, y: int, q: int):
        super(Curve, self).__init__()
        # elliptic curve in short Weierstrass form
        # (you can find information about what each of the parameters does
        # in NIST SP 800-186)
        self.p = p
        self.a = a
        self.b = b
        self.x = x
        self.y = y
        self.q = q

    def oncurve(self, x1: int, y1: int):
        # verify that the given point satisfies the curve equation
        x2 = (x1**3 + self.a * x1 + self.b) % self.p
        y2 = y1**2 % self.p
        return x2 == y2


class Point(object):
    def __init__(self, x: int, y: int, curve: Curve):
        super(Point, self).__init__()
        # A point on an elliptic curve.
        self.x = x
        self.y = y
        self.curve = curve

    @classmethod
    def copy(cls, p):
        return Point(p.x, p.y, p.curve)

    @classmethod
    def inf(cls, curve: Curve):
        # the point (None, None) identifies the point at infinity
        return Point(None, None, curve)

    def oncurve(self) -> bool:
        return self.curve.oncurve(self.x, self.y)

    def pointdouble(self):
        if self.y == 0:
            return Point(None, None, self.curve)
        s = (3 * self.x**2 + self.curve.a) % self.curve.p
        s = s * inv_euklid(self.y * 2, self.curve.p)
        x_new = (s**2 - 2 * self.x) % self.curve.p
        y_new = (s * (self.x - x_new) - self.y) % self.curve.p
        return Point(x_new, y_new, self.curve)

    def pointaddition(self, p2):
        if self.y == -p2.y % self.curve.p and self.x == p2.x:
            return Point(None, None, self.curve)
        s = (p2.y - self.y) * inv_euklid(p2.x - self.x, self.curve.p) % self.curve.p
        x_new = (s**2 - self.x - p2.x) % self.curve.p
        y_new = (s * (self.x - x_new) - self.y) % self.curve.p
        return Point(x_new, y_new, self.curve)

    def get_coordinates(self) -> tuple[int, int]:
        return (self.x, self.y)

    def get_inverse(self):
        return Point(self.x, -self.y % self.curve.p, self.curve)

    def __repr__(self) -> str:
        return (
            "x: %s\ny: %s" % (hex(self.x), hex(self.y))
            if type(self.x) is int
            else "x: O\ny: O"
        )

    ## overloaded operators:

    def __eq__(self, obj):
        if (self.x == None and obj.x != None) or (self.x != None and obj.x == None):
            return False
        if self.x == None and obj.x == None:
            return True
        return self.x == obj.x % obj.curve.p and self.y == obj.y % obj.curve.p

    def __add__(self, obj):
        assert self.curve == obj.curve, "Points not on the same curve"
        if self.x == None:
            return obj
        if obj.x == None:
            return self
        if self == obj:
            return self.pointdouble()
        else:
            return self.pointaddition(obj)

    def __sub__(self, obj):
        if obj.x == None:
            return self
        inv = Point(obj.x, -obj.y % obj.curve.p, obj.curve)
        return self.__add__(inv)

    def __mul__(self, a):
        return naf_double_and_add(self, a)

    def __rmul__(self, a):
        return self.__mul__(self, a)


def inv_euklid(a, b):
    x0, x1, y0, y1 = 0, 1, 1, 0
    a = a % b
    m = b
    while a != 0:
        q = b // a
        t = a
        a = b % a
        b = t
        t = y0
        y0 = y1
        y1 = t - q * y1
        t = x0
        x0 = x1
        x1 = t - q * x1
    return x0 % m


# Helper function for task 2c
def plot_performance(
    perf_normal: list[float], perf_naf: list[float], output_filename="ecdh_speed.png"
) -> None:
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex="all")

    fig.suptitle("Execution speed for NAF and standard Double-and-Add")

    # if you're interested in how this works, you can find examples here:
    # https://matplotlib.org/stable/gallery/statistics/hist.html#sphx-glr-gallery-statistics-hist-py
    ax_top.hist(perf_normal, bins=50)
    ax_bottom.hist(perf_naf, bins=50)

    avg_standard = sum(perf_normal) / len(perf_normal)
    avg_naf = sum(perf_naf) / len(perf_naf)

    ax_top.set_title(f"Standard, {avg_standard:.2f} s average")
    ax_top.set_xlabel("time (seconds)")
    ax_top.set_ylabel("count")  # or "frequency"
    ax_bottom.set_title(f"NAF, {avg_naf:.2f} s average")
    ax_bottom.set_xlabel("time (seconds)")
    ax_bottom.set_ylabel("count")  # or "frequency"

    plt.tight_layout()
    plt.savefig(output_filename, format="png", dpi=300)
    plt.show()


########## ENTER YOUR SOLUTION BELOW ##########


# Task 1b
def double_and_add(p: Point, a: int) -> Point:
    # TODO implement
    ...


# Task 1c
def calc_naf_representation(exponent: int) -> list[int]:
    # TODO implement
    # In the returned list, the least significant bit (LSB) should be at index 0.
    ...


# Task 1d
def naf_double_and_add(p: Point, a: int) -> Point:
    naf_array = calc_naf_representation(a)
    # TODO implement
    ...


# DO NOT remove this line
if __name__ == "__main__":
    # setup
    nistp256 = Curve(
        0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF,
        0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC,
        0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B,
        0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
        0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5,
        0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
    )

    # Task 1d,e
    point_1 = Point(
        0xB0A3D57E543778709F9D74910DEB1B5B2C6405C87FFF60A175A68E866B0388B3,
        0xCF2879FD7EED0AC3000C484C488A6256C3887F39AFDD62E764A74E814C7E7EBE,
        nistp256,
    )

    k = 0xFEEDDEADBEEFBADCAB1EDECAFBADC0DEC001D00DC0DEBA5EC0CAC01AADD511FE

    # TODO multiply point and scalar using NAF
    # TODO report number of operations needed

    # Task 2a
    kprivA = 0x7A041F183A600B1D34E2D20D0FD994E003257481ACF9A107184599BF43F75AE7
    kprivB = 0x85249CDDD97CCFABC10400F3E0C03207721E7A079CDB93928C88DAD7204C2B12

    # TODO calculate shared key
    kpubA = None
    kpubB = None

    T_AB = None

    # Task 2b

    n = 10
    # TODO measure time for 1000 executions of double_and_add
    # TODO measure time for 1000 executions of naf_double_and_add

    # You can use perf_counter() to get the current time

    # Task 2c

    n = 5000

    # TODO measure time for each execution of double_and_add
    time_normal_list = []
    # TODO measure time for each execution of naf_double_and_add
    time_naf_list = []

    plot_performance(time_normal_list, time_naf_list, output_filename="ecdh_speed.png")

    # Task 2d

    # TODO calculate the AES key using T_AB
