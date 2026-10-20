import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from canonical_observer import frozen_cohort_mask


class CohortMembershipTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'cohort.npz'
        environment = patch.dict(os.environ, DSPH_PREEXISTING_ID_PATH=str(self.path))
        environment.start()
        self.addCleanup(environment.stop)

    def save(self, ids):
        np.savez_compressed(self.path, particle_id=ids)

    def test_unordered_uint64_membership_matches_independent_set(self):
        cohort = np.array([2**64-1, 2**60+7, 0, 33], dtype=np.uint64)
        current = np.array([81, 33, 2**60+7, 2**60+8, 0, 2**64-1], dtype=np.uint64)
        self.save(cohort)
        before = current.copy()
        mask, state = frozen_cohort_mask(current, 'unused')
        members = {int(value) for value in cohort}
        expected = np.array([int(value) in members for value in current])
        np.testing.assert_array_equal(mask, expected)
        np.testing.assert_array_equal(current, before)
        self.assertEqual(state['cohort_N'], len(cohort))

    def test_missing_count_is_exact_including_empty_snapshot(self):
        self.save(np.array([1, 2, 3], dtype=np.uint64))
        for values, missing in (([1, 8], 2), ([1, 2], 1), ([], 3)):
            with self.subTest(values=values), self.assertRaisesRegex(ValueError, 'count='+str(missing)):
                frozen_cohort_mask(np.array(values, dtype=np.uint64), 'unused')

    def test_duplicate_ids_still_fail_before_completeness_check(self):
        self.save(np.array([1, 2], dtype=np.uint64))
        with self.assertRaisesRegex(ValueError, 'unique'):
            frozen_cohort_mask(np.array([1, 1, 2], dtype=np.uint64), 'unused')
        self.save(np.array([1, 1], dtype=np.uint64))
        with self.assertRaisesRegex(ValueError, 'unique'):
            frozen_cohort_mask(np.array([1, 2], dtype=np.uint64), 'unused')

    def test_empty_cohort_and_wrong_integer_dtype_still_fail(self):
        self.save(np.array([], dtype=np.uint64))
        with self.assertRaisesRegex(ValueError, 'empty'):
            frozen_cohort_mask(np.array([1], dtype=np.uint64), 'unused')
        self.save(np.array([1], dtype=np.uint64))
        with self.assertRaisesRegex(ValueError, 'uint64'):
            frozen_cohort_mask(np.array([1], dtype=np.int64), 'unused')


if __name__ == '__main__':
    unittest.main()
