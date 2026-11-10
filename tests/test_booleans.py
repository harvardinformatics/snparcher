import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "lib"))

from booleans import parse_bool_column  # noqa: E402


def test_accepts_words_letters_and_1_0_in_any_case():
    values = pd.Series(["TRUE", "f", "Yes", "n", "T", "no", 1, 0, True, False, 1.0])
    assert parse_bool_column(values, "exclude").tolist() == [
        True, False, True, False, True, False, True, False, True, False, True,
    ]


@pytest.mark.parametrize("default", [False, True])
def test_missing_and_blank_values_take_the_default(default):
    values = pd.Series([None, float("nan"), "", "  "], dtype=object)
    assert parse_bool_column(values, "mark_duplicates", default=default).tolist() == [default] * 4


@pytest.mark.parametrize("value", [2, 0.5, "maybe"])
def test_other_values_raise_with_the_csv_row(value):
    with pytest.raises(ValueError, match=r"Invalid exclude value at row 3: "):
        parse_bool_column(pd.Series(["yes", value], dtype=object), "exclude")
