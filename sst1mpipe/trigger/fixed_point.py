"""Fixed-point arithmetic of the TDSCAN firmware (Vitis HLS ``ap_fixed``), in NumPy.

A number is stored as an integer code: value = code / 2**frac_bits.
Formats are written like ``SQ5.2`` (signed, 5 integer bits including the
sign, 2 fractional bits) or ``UQ4.0`` (unsigned).
"""
import re
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Format:
    signed: bool
    int_bits: int
    frac_bits: int

    @property
    def bits(self):
        return self.int_bits + self.frac_bits

    @property
    def min_code(self):
        return -(2 ** (self.bits - 1)) if self.signed else 0

    @property
    def max_code(self):
        return 2 ** (self.bits - 1) - 1 if self.signed else 2**self.bits - 1


def parse(text):
    """``"SQ5.2"`` -> Format(signed=True, int_bits=5, frac_bits=2)."""
    match = re.fullmatch(r"([SU])Q(\d+)\.(\d+)", text)
    if match is None:
        raise ValueError(f"Fixed-point format {text!r} must look like SQ5.2 or UQ4.0")
    return Format(signed=match.group(1) == "S", int_bits=int(match.group(2)), frac_bits=int(match.group(3)))


def fit(codes, fmt, overflow_mode):
    """Bring codes into the range of ``fmt``: saturate (AP_SAT) or wrap around (AP_WRAP)."""
    codes = np.asarray(codes, dtype=np.int64)
    if overflow_mode == "AP_SAT":
        return np.clip(codes, fmt.min_code, fmt.max_code)
    if overflow_mode == "AP_WRAP":
        codes = np.mod(codes, 2**fmt.bits)
        if fmt.signed:
            codes = np.where(codes > fmt.max_code, codes - 2**fmt.bits, codes)
        return codes
    raise ValueError(f"overflow_mode must be AP_SAT or AP_WRAP, got {overflow_mode!r}")


def to_codes(values, fmt, overflow_mode, quantization_mode):
    """Real values -> codes. AP_TRN truncates (floor), AP_RND rounds half up."""
    scaled = np.asarray(values, dtype=np.float64) * 2**fmt.frac_bits
    if quantization_mode == "AP_RND":
        scaled = scaled + 0.5
    elif quantization_mode != "AP_TRN":
        raise ValueError(f"quantization_mode must be AP_TRN or AP_RND, got {quantization_mode!r}")
    return fit(np.floor(scaled), fmt, overflow_mode)


def to_real(codes, fmt):
    """Codes -> real values."""
    return np.asarray(codes, dtype=np.float64) / 2**fmt.frac_bits


def requantize(codes, source, target, shift, overflow_mode):
    """Move codes from a wide register ``source`` to a narrower register ``target``.

    Keeps the most significant bits: arithmetic right shift by
    ``source.bits - target.bits - shift``, then saturate/wrap to ``target``.
    Same register on both sides: no shift, as in the firmware.
    """
    if source == target:
        return fit(codes, target, overflow_mode)
    codes = fit(codes, source, overflow_mode)
    right_shift = source.bits - target.bits - shift
    if right_shift >= 0:
        codes = np.floor_divide(codes, 2**right_shift)
    else:
        codes = codes * 2 ** (-right_shift)
    return fit(codes, target, overflow_mode)
