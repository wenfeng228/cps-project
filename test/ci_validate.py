import sys
import pandas as pd
import numpy as np

try:
    actual = pd.read_csv('actual_output.csv', header=None)
    expected = pd.read_csv('test/expected_output.csv', header=None)

    # Drop columns 3, 4, 5 (0-based index corresponds to columns 4, 5, 6)
    actual = actual.drop(columns=[3, 4, 5])
    expected = expected.drop(columns=[3, 4, 5])

    # Compare with tolerance (1e-5)
    pd.testing.assert_frame_equal(actual, expected, check_dtype=False, rtol=1e-5)
    print("SUCCESS: Output matches expected data within tolerance.")

except AssertionError as e:
    print("FAILURE: Data mismatch found.")
    print(e)
    sys.exit(1)
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)
