"""Parse true/false columns of user CSVs (sample sheet and sample metadata)."""

import pandas as pd

TRUTHY = {"true", "t", "yes", "y", "1"}
FALSY = {"false", "f", "no", "n", "0"}


def parse_bool_column(series, column_name, default=False):
    """Parse a column of true/false, t/f, yes/no, y/n or 1/0, in any case.

    Missing or blank values become `default`. Any other value raises, naming
    its CSV row (the header is row 1).
    """
    parsed = []
    for idx, value in series.items():
        text = str(value).strip().lower()
        if pd.isna(value) or not text:
            parsed.append(default)
        elif isinstance(value, bool):
            parsed.append(value)
        elif isinstance(value, (int, float)) and value in (0, 1):
            parsed.append(bool(value))
        elif text in TRUTHY:
            parsed.append(True)
        elif text in FALSY:
            parsed.append(False)
        else:
            raise ValueError(
                f"Invalid {column_name} value at row {idx + 2}: {value!r}. "
                "Use true/false (or equivalent boolean values)."
            )
    return pd.Series(parsed, index=series.index, dtype=bool)
