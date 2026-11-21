"""
Read and write the channel matrices file (JSON), the input of iondude.

The in-memory form is a dict with the same shape as the JSON file. The matrices
(counts, channel and inverse) are numpy arrays. The format is in iondude/matrices/README.md.
"""

import json
import math
import re
import uuid
from datetime import date
from pathlib import Path
from typing import IO

import numpy as np

FORMAT = "iondude-matrices"
FORMAT_VERSION = 1
SLOT_TYPES = ("position", "quality")
SOURCE_FIELDS = ("tool", "version", "date", "input")
MATRICES = ("counts", "channel", "inverse")

# version 1 has fixed symbols, iondude uses them in this order
SYMBOLS = {
    "substitutions": ["A", "C", "G", "T"],
    "homopolymers": list(range(1, 11)),
}

def write_matrices(file: IO[str], data: dict) -> None:
    """
    Write the matrices as JSON. Each matrix row is on one line, so that a diff of two files is readable.

    :param file: an open text file
    :param data: the matrices, as read_matrices returns them
    :raises ValueError: if a matrix has a NaN or an infinite value
    """
    inline = []
    # json.dumps cannot put only some lists on one line: put a placeholder, then replace it.
    # The token is not in the metadata, so no metadata string can match a placeholder
    metadata = json.dumps({name: value for name, value in data.items() if name not in SYMBOLS})
    token = uuid.uuid4().hex
    while token in metadata:
        token = uuid.uuid4().hex

    def one_line(values: list) -> str:
        inline.append(json.dumps(values, allow_nan=False))
        return f"@{token}:{len(inline) - 1}@"

    document = dict(data)
    for kind in SYMBOLS:
        section = data[kind]
        document[kind] = {
            "symbols": one_line(list(section["symbols"])),
            "slots": [
                {"min": slot["min"], "max": slot["max"],
                 **{name: [one_line(row) for row in np.asarray(slot[name]).tolist()] for name in MATRICES}}
                for slot in section["slots"]
            ],
        }

    text = json.dumps(document, indent=2)
    placeholder = re.compile(f'"@{token}:(\\d+)@"')
    file.write(placeholder.sub(lambda match: inline[int(match.group(1))], text) + "\n")


def read_matrices(path: str | Path) -> dict:
    """
    Read a matrices file and check all the rules of the format.

    :param path: the JSON file
    :return: the matrices, with numpy arrays for counts (int), channel and inverse (float)
    :raises ValueError: if the file does not follow the format, with the field and the slot
    """
    with open(path) as file:
        data = json.load(file)
    _check(data)
    for kind in SYMBOLS:
        for slot in data[kind]["slots"]:
            slot["counts"] = np.array(slot["counts"], dtype=np.int64)
            slot["channel"] = np.array(slot["channel"], dtype=float)
            slot["inverse"] = np.array(slot["inverse"], dtype=float)
    return data


def get_slot(data: dict, kind: str, value: int) -> dict:
    """
    Get the slot that holds a read position (slot_type position) or a quality score (slot_type quality).

    :param data: the matrices, as read_matrices returns them
    :param kind: "substitutions" or "homopolymers"
    :param value: the read position or the quality score
    :return: the slot, a dict with min, max, counts, channel and inverse
    :raises ValueError: if no slot holds the value
    """
    slots = data[kind]["slots"]
    index = value // data["range"]
    if value < 0 or index >= len(slots):
        raise ValueError(f"{kind}: no slot for {data['slot_type']} {value}")
    return slots[index]


def _fail(field: str, message: str) -> None:
    raise ValueError(f"{field}: {message}")


def _is_int(value) -> bool:
    # bool is a subclass of int, but true and false are not numbers in the file
    return isinstance(value, int) and not isinstance(value, bool)


def _check_matrix(field: str, matrix, dim: int, integers: bool) -> None:
    if not isinstance(matrix, list) or len(matrix) != dim or \
            any(not isinstance(row, list) or len(row) != dim for row in matrix):
        _fail(field, f"must be {dim} x {dim}")
    for row in matrix:
        for value in row:
            if integers and not _is_int(value):
                _fail(field, "must have only integers")
            if not integers and (not (_is_int(value) or isinstance(value, float)) or not math.isfinite(value)):
                _fail(field, "must have only finite numbers")


def _check(data) -> None:
    if not isinstance(data, dict):
        _fail("file", "must be a JSON object")
    if data.get("format") != FORMAT:
        _fail("format", f'must be "{FORMAT}"')
    if not _is_int(data.get("format_version")) or data["format_version"] != FORMAT_VERSION:
        _fail("format_version", f"must be {FORMAT_VERSION}, the file has {data.get('format_version')!r}")

    source = data.get("source")
    if not isinstance(source, dict):
        _fail("source", "missing")
    for name in SOURCE_FIELDS:
        if not isinstance(source.get(name), str):
            _fail(f"source.{name}", "missing")
    try:
        # fromisoformat accepts more forms than YYYY-MM-DD in Python 3.11+, iondude does not
        if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", source["date"]):
            raise ValueError
        date.fromisoformat(source["date"])
    except ValueError:
        _fail("source.date", "must be an ISO 8601 date (YYYY-MM-DD)")
    if not isinstance(data.get("estimation"), dict):
        _fail("estimation", "missing")

    slot_type = data.get("slot_type")
    if slot_type not in SLOT_TYPES:
        _fail("slot_type", f"must be one of {', '.join(SLOT_TYPES)}")
    width = data.get("range")
    if not _is_int(width) or width < 1:
        _fail("range", "must be an integer greater than 0")
    if slot_type == "quality" and width != 1:
        _fail("range", "must be 1 for slot_type quality")
    max_length = data.get("max_length")
    if not _is_int(max_length) or max_length < 1:
        _fail("max_length", "must be an integer greater than 0")
    n_slots = -(-max_length // width)

    for kind, symbols in SYMBOLS.items():
        section = data.get(kind)
        if not isinstance(section, dict):
            _fail(kind, "missing")
        if section.get("symbols") != symbols:
            _fail(f"{kind}.symbols", f"must be {json.dumps(symbols)}")
        slots = section.get("slots")
        if not isinstance(slots, list) or len(slots) != n_slots:
            _fail(f"{kind}.slots", f"must have {n_slots} slots, from 0 to {max_length}")
        for i, slot in enumerate(slots):
            field = f"{kind}.slots[{i}]"
            if not isinstance(slot, dict):
                _fail(field, "must be an object")
            if slot.get("min") != i * width or slot.get("max") != (i + 1) * width:
                _fail(field, f"min and max must be {i * width} and {(i + 1) * width}, no gaps and no overlaps")
            for name in MATRICES:
                _check_matrix(f"{field}.{name}", slot.get(name), len(symbols), name == "counts")
