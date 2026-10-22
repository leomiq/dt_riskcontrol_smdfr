"""Verify the portable poly4 constants against the archived cfit MAT file.

This parser is deliberately pinned to the exact released MAT-file SHA-256.
It reads the five scalar doubles in the MCOS coeffValues cell array and
does not depend on SciPy or Curve Fitting Toolbox.
"""
from pathlib import Path
import hashlib
import math
import re
import struct
import zlib

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "models" / "fit_poly.mat"
SOURCE = ROOT / "matlab" / "fit_poly_coefficients.m"
EXPECTED_SHA256 = (
    "cda7e1973a7529aafbe5b9cf453c4f0b5c258bf3868b3223ca3ce4b34e934c6d"
)


def archived_coefficients():
    data = ARCHIVE.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError("fit_poly.mat differs from the pinned v1.0.0 archive")
    pos = 128  # MATLAB 5 header
    chunks = []
    while pos + 8 <= len(data):
        kind, length = struct.unpack_from("<II", data, pos)
        if kind != 15 or pos + 8 + length > len(data):
            raise ValueError("Unexpected MAT v5 compressed element")
        chunks.append(zlib.decompress(data[pos + 8 : pos + 8 + length]))
        pos += 8 + length
    if len(chunks) != 2:
        raise ValueError("Unexpected MAT v5 element count")
    payload = chunks[1]
    if b"coeffValues" not in payload or b"poly4" not in payload:
        raise ValueError("Expected cfit poly4 metadata is absent")
    values = []
    for pos in (1496, 1560, 1624, 1688, 1752):
        kind, length = struct.unpack_from("<II", payload, pos)
        array_class = struct.unpack_from("<I", payload, pos + 16)[0]
        rows, columns = struct.unpack_from("<II", payload, pos + 32)
        value_kind, value_length = struct.unpack_from("<II", payload, pos + 48)
        if (kind, length, array_class, rows, columns, value_kind, value_length) != (
            14, 56, 6, 1, 1, 9, 8
        ):
            raise ValueError("Unexpected coeffValues scalar layout")
        values.append(struct.unpack_from("<d", payload, pos + 56)[0])
    return values


def source_coefficients():
    source = SOURCE.read_text(encoding="utf-8")
    match = re.search(r"coefficients\s*=\s*\[(.*?)\];", source, re.DOTALL)
    if not match:
        raise ValueError("No coefficients vector in MATLAB source")
    numbers = re.findall(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        match.group(1),
    )
    return [float(number) for number in numbers]


if __name__ == "__main__":
    archived = archived_coefficients()
    current = source_coefficients()
    if len(current) != 5 or not all(
        math.isclose(a, b, rel_tol=0, abs_tol=0)
        for a, b in zip(archived, current)
    ):
        raise SystemExit("FAIL: MATLAB coefficients differ from the archive")
    print("PASS: five poly4 coefficients match the pinned cfit archive exactly")
