import copy
import io
import json
import os
import tempfile
import unittest
import uuid
from unittest import mock

import numpy as np

from estimator.matrices import get_slot, read_matrices, write_matrices


def make_slot(low, high, dim):
    counts = np.eye(dim, dtype=np.int64) * 9
    counts[0][1] = 1
    channel = counts / counts.sum(axis=1)[:, None]
    return {"min": low, "max": high, "counts": counts, "channel": channel, "inverse": np.linalg.inv(channel)}


def make_matrices():
    """
    A valid file: 3 slots of 20 positions, the last one past max_length
    """
    return {
        "format": "iondude-matrices",
        "format_version": 1,
        "source": {"tool": "iondude-estimator", "version": "0.4", "date": "2026-11-21", "input": "test.sam"},
        "estimation": {"max_read_length": 50, "range": 20, "min_mapping_quality": 20},
        "slot_type": "position",
        "range": 20,
        "max_length": 50,
        "substitutions": {
            "symbols": ["A", "C", "G", "T"],
            "slots": [make_slot(low, low + 20, 4) for low in (0, 20, 40)],
        },
        "homopolymers": {
            "symbols": list(range(1, 11)),
            "slots": [make_slot(low, low + 20, 10) for low in (0, 20, 40)],
        },
    }


def to_json(data):
    out = io.StringIO()
    write_matrices(out, data)
    return json.loads(out.getvalue())


class MatricesTests(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "matrices.json")

    def tearDown(self):
        self.tmp.cleanup()

    def save(self, document):
        with open(self.path, "w") as file:
            json.dump(document, file)

    def assert_invalid(self, change, field):
        """
        The file with the change does not load, and the error names the field
        """
        document = to_json(make_matrices())
        change(document)
        self.save(document)
        with self.assertRaises(ValueError) as error:
            read_matrices(self.path)
        self.assertIn(field, str(error.exception))

    def test_write_and_read_back(self):
        data = make_matrices()
        with open(self.path, "w") as file:
            write_matrices(file, data)

        loaded = read_matrices(self.path)

        self.assertEqual(loaded["source"], data["source"])
        self.assertEqual(loaded["estimation"], data["estimation"])
        for kind in ("substitutions", "homopolymers"):
            for expected, slot in zip(data[kind]["slots"], loaded[kind]["slots"], strict=True):
                self.assertEqual(slot["counts"].dtype, np.int64)
                np.testing.assert_array_equal(slot["counts"], expected["counts"])
                # full precision: the values are the same, not only close
                np.testing.assert_array_equal(slot["channel"], expected["channel"])
                np.testing.assert_array_equal(slot["inverse"], expected["inverse"])

    def test_each_row_on_one_line(self):
        out = io.StringIO()
        write_matrices(out, make_matrices())

        lines = out.getvalue().splitlines()

        self.assertIn('    "symbols": ["A", "C", "G", "T"],', lines)
        self.assertIn("          [9, 1, 0, 0],", lines)

    def test_metadata_like_a_placeholder(self):
        """
        Metadata strings that look like the internal row placeholders are written as they are
        """
        data = make_matrices()
        data["source"]["input"] = "@inline0@"
        data["estimation"]["note"] = "@inline999999@"
        with open(self.path, "w") as file:
            write_matrices(file, data)

        loaded = read_matrices(self.path)

        self.assertEqual(loaded["source"]["input"], "@inline0@")
        self.assertEqual(loaded["estimation"]["note"], "@inline999999@")

    def test_token_in_metadata_is_not_used(self):
        """
        A metadata string that has the placeholder token does not change: the writer takes a new token
        """
        used = uuid.UUID(int=1)
        data = make_matrices()
        data["source"]["input"] = f"@{used.hex}:0@"
        with mock.patch("estimator.matrices.uuid.uuid4", side_effect=[used, uuid.UUID(int=2)]):
            out = io.StringIO()
            write_matrices(out, data)

        self.assertEqual(json.loads(out.getvalue())["source"]["input"], f"@{used.hex}:0@")

    def test_nan_is_not_written(self):
        data = make_matrices()
        data["substitutions"]["slots"][0]["channel"][0][0] = np.nan
        with self.assertRaises(ValueError):
            write_matrices(io.StringIO(), data)

    def test_get_slot(self):
        data = make_matrices()

        self.assertEqual(get_slot(data, "substitutions", 0)["min"], 0)
        self.assertEqual(get_slot(data, "substitutions", 25)["min"], 20)
        self.assertEqual(get_slot(data, "homopolymers", 49)["min"], 40)
        with self.assertRaises(ValueError):
            get_slot(data, "substitutions", 60)
        with self.assertRaises(ValueError):
            get_slot(data, "substitutions", -1)

    def test_get_slot_quality(self):
        data = make_matrices()
        data.update(slot_type="quality", range=1, max_length=3)
        for kind in ("substitutions", "homopolymers"):
            data[kind]["slots"] = [make_slot(q, q + 1, len(data[kind]["symbols"])) for q in range(3)]
        with open(self.path, "w") as file:
            write_matrices(file, data)

        loaded = read_matrices(self.path)

        self.assertEqual(get_slot(loaded, "substitutions", 2)["min"], 2)

    def test_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            read_matrices(os.path.join(self.tmp.name, "missing.json"))

    def test_not_json(self):
        with open(self.path, "w") as file:
            file.write("#SF 0-20 Substitutions frequencies\n")
        with self.assertRaises(ValueError):
            read_matrices(self.path)

    def test_wrong_format(self):
        self.assert_invalid(lambda d: d.update(format="other"), "format")

    def test_wrong_format_version(self):
        self.assert_invalid(lambda d: d.update(format_version=2), "format_version")

    def test_missing_source_field(self):
        self.assert_invalid(lambda d: d["source"].pop("input"), "source.input")

    def test_source_date_not_iso(self):
        self.assert_invalid(lambda d: d["source"].update(date="21/11/2026"), "source.date")

    def test_missing_estimation(self):
        self.assert_invalid(lambda d: d.pop("estimation"), "estimation")

    def test_wrong_slot_type(self):
        self.assert_invalid(lambda d: d.update(slot_type="read"), "slot_type")

    def test_wrong_range(self):
        self.assert_invalid(lambda d: d.update(range=0), "range")

    def test_quality_range_must_be_1(self):
        self.assert_invalid(lambda d: d.update(slot_type="quality"), "range")

    def test_wrong_max_length(self):
        self.assert_invalid(lambda d: d.update(max_length="50"), "max_length")

    def test_missing_section(self):
        self.assert_invalid(lambda d: d.pop("homopolymers"), "homopolymers")

    def test_wrong_symbols(self):
        self.assert_invalid(lambda d: d["substitutions"].update(symbols=["A", "C", "G"]), "substitutions.symbols")

    def test_missing_slot(self):
        self.assert_invalid(lambda d: d["homopolymers"]["slots"].pop(), "homopolymers.slots")

    def test_slot_gap(self):
        def change(document):
            document["substitutions"]["slots"][1] = copy.deepcopy(document["substitutions"]["slots"][2])
        self.assert_invalid(change, "substitutions.slots[1]")

    def test_slot_wrong_width(self):
        self.assert_invalid(lambda d: d["substitutions"]["slots"][0].update(max=10), "substitutions.slots[0]")

    def test_matrix_wrong_shape(self):
        self.assert_invalid(lambda d: d["substitutions"]["slots"][2]["channel"].pop(),
                            "substitutions.slots[2].channel")

    def test_missing_matrix(self):
        self.assert_invalid(lambda d: d["homopolymers"]["slots"][0].pop("inverse"), "homopolymers.slots[0].inverse")

    def test_counts_not_integers(self):
        self.assert_invalid(lambda d: d["substitutions"]["slots"][0]["counts"][0].__setitem__(0, 9.5),
                            "substitutions.slots[0].counts")

    def test_channel_not_numbers(self):
        self.assert_invalid(lambda d: d["substitutions"]["slots"][0]["channel"][0].__setitem__(0, "x"),
                            "substitutions.slots[0].channel")
