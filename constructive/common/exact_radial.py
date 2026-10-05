"""Small exact acceptance kernel. Only built-in integer arithmetic accepts.

Coordinates arrive as actual binary64 values. A norm-squared is represented
as an integer numerator over one common power of two per configuration.
"""
import ast
import hashlib
import math
import struct


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def read_npy(data, expected_shape):
    """Decode uncompressed C-order binary64 NPY; no pickle or code loading."""
    if data[:6] != b"\x93NUMPY":
        raise ValueError("not NPY")
    version = tuple(data[6:8])
    if version == (1, 0):
        header_length = struct.unpack_from("<H", data, 8)[0]
        start = 10
    elif version in ((2, 0), (3, 0)):
        header_length = struct.unpack_from("<I", data, 8)[0]
        start = 12
    else:
        raise ValueError("unsupported NPY version")
    header = ast.literal_eval(data[start:start + header_length].decode("utf-8"))
    if set(header) != {"descr", "fortran_order", "shape"}:
        raise ValueError("unexpected NPY header fields")
    if header["fortran_order"] is not False:
        raise ValueError("only C order is accepted")
    if tuple(header["shape"]) != tuple(expected_shape):
        raise ValueError("NPY shape mismatch")
    if header["descr"] not in ("<f8", ">f8"):
        raise ValueError("only explicit-endian binary64 is accepted")
    offset = start + header_length
    size = math.prod(expected_shape) * 8
    if len(data) != offset + size:
        raise ValueError("NPY payload length mismatch")
    return memoryview(data)[offset:], header["descr"][0]


def centers_at(payload, endian, row=0):
    values = struct.unpack_from(endian + "200d", payload, row * 1600)
    if not all(math.isfinite(x) for x in values):
        raise ValueError("nonfinite coordinate")
    return [(values[2 * i], values[2 * i + 1]) for i in range(100)]


def coordinate_digest(points):
    # Canonical source chunks already contain lexicographically sorted centers.
    return sha256(b"".join(struct.pack("<2d", x, y) for x, y in points))


def exact_norm_squares(points):
    """Return N_i and q with ||p_i||^2 = N_i / 2^q exactly."""
    if len(points) != 100:
        raise ValueError("exactly all100 centers required")
    ratios = []
    exponent = 0
    for x, y in points:
        row = []
        for value in (x, y):
            if not math.isfinite(value):
                raise ValueError("nonfinite coordinate")
            numerator, denominator = float(value).as_integer_ratio()
            power = denominator.bit_length() - 1
            if denominator != 1 << power:
                raise ValueError("binary64 denominator is not dyadic")
            row.append((numerator, power))
            exponent = max(exponent, power)
        ratios.append(row)
    numerators = []
    for (x, ex), (y, ey) in ratios:
        X, Y = x << (exponent - ex), y << (exponent - ey)
        numerators.append(X * X + Y * Y)
    return numerators, 2 * exponent


def count_mask(norm_squares, threshold_numerator, threshold_denominator):
    """Count ALL100 exact ||p|| <= nonnegative rational threshold.

    N_i D^2 <= u^2 2^q is precisely ||p_i|| <= u/D for u>=0.
    Bit i records inclusion of center i, including equality.
    """
    if threshold_numerator < 0 or threshold_denominator <= 0:
        raise ValueError("invalid nonnegative rational threshold")
    numerators, exponent = norm_squares
    denominator_squared = threshold_denominator * threshold_denominator
    right = (threshold_numerator * threshold_numerator) << exponent
    mask = 0
    for index, numerator in enumerate(numerators):
        if numerator * denominator_squared <= right:
            mask |= 1 << index
    if len(numerators) != 100 or mask >> 100:
        raise ValueError("not an all100 mask")
    return mask.bit_count(), f"{mask:025x}"


def hall_accept(low, high, threshold_numerator, denominator, epsilon_numerator):
    if epsilon_numerator <= 0:
        raise ValueError("epsilon must be positive")
    low_count, low_mask = count_mask(low, threshold_numerator, denominator)
    high_count, high_mask = count_mask(
        high, threshold_numerator + epsilon_numerator, denominator)
    return {
        "low_count": low_count, "high_count": high_count,
        "low_mask_hex": low_mask, "high_mask_hex": high_mask,
        "strict_count_deficit": low_count - high_count,
        "accepted": low_count > high_count,
    }
