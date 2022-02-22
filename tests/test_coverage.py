import numpy as np
import pytest
from mbf_sampledata import get_sample_path
from mbf_bam import calculate_coverage, calculate_coverage_sum


class TestCoverage:
    def test_simple(self):
        input = str(get_sample_path("mbf_align/chipseq_chr22.bam"))
        start = 16_000_000
        stop = start + 1_000_000
        intervals = [("chr22", start, stop, False)]
        forward = calculate_coverage(input, None, intervals)
        assert len(forward) == 1
        assert len(forward[0]) == stop - start
        assert np.sum(forward[0]) > 0
        intervals = [("chr22", start, stop, True)]
        reverse = calculate_coverage(input, None, intervals)
        assert np.sum(reverse[0]) == np.sum(forward[0])
        eq = forward[0][::-1] == reverse[0]
        assert eq

    def test_overlapping(self):
        input = str(get_sample_path("mbf_align/chipseq_chr22.bam"))
        offset = 16_000_000
        k = 1_000_000  # first 100k are empty...
        intervals = [
            ("chr22", offset, offset + k, False),
            ("chr22", offset + k // 2, offset + k // 2 + k, False),
        ]
        forward = calculate_coverage(input, None, intervals)
        assert len(forward) == 2
        assert len(forward[0]) == k
        assert len(forward[1]) == k
        assert np.sum(forward[0]) > 0
        assert np.sum(forward[1]) > 0
        eq = forward[0][500_000:] == forward[1][:500_000]
        assert eq
        assert np.sum(forward[0][500_000:]) > 0

    def test_example(self):
        input = str(get_sample_path("mbf_align/chipseq_chr22.bam"))
        intervals = [
            ("chr22", 16097540 - 1, 16097652 - 1, False),
        ]
        forward = calculate_coverage(input, None, intervals)
        assert np.sum(forward[0]) == 36

        intervals = [
            ("chr22", 16097552 - 1, 16097552 - 1 + 36, False),
        ]
        forward = calculate_coverage(input, None, intervals)
        assert np.sum(forward[0]) == 36

        intervals = [
            ("chr22", 16097552 - 1 + 3, 16097552 - 1 + 36 - 2, False),
        ]
        forward = calculate_coverage(input, None, intervals)
        assert np.sum(forward[0]) == 36 - 3 - 2
        assert np.sum(forward[0]) == len(forward[0])

        intervals = [
            ("chr22", 51038715 - 1, 51038746, False),
        ]
        forward = calculate_coverage(input, None, intervals)
        # fmt: off
        assert forward[0] == [ 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, ]
        # fmt: on

    def test_cov_sum(self):
        input = str(get_sample_path("mbf_align/chipseq_chr22.bam"))
        start = 16_000_000
        stop = start + 1_000_000
        intervals = []
        for ss in range(start, stop, 1000):
            st = ss + 1000
            intervals.append(("chr22", ss, st, False))
        forward = calculate_coverage(input, None, intervals)
        forward_sum = calculate_coverage_sum(input, None, intervals)
        npf = np.array(forward)
        assert npf.sum().sum() > 0
        assert npf.shape[0] == len(intervals)
        assert len(forward_sum) == 1000
        assert (npf.sum(axis=0) == forward_sum).all()

    def test_cov_sum_raises_on_unequal_sizes(self):

        input = str(get_sample_path("mbf_align/chipseq_chr22.bam"))
        intervals = [
            ("chr22", 0, 1000, False),
            ("chr22", 1000, 1500, False),
        ]
        with pytest.raises(ValueError):
            calculate_coverage_sum(input, None, intervals)

    def test_negative(self):
        input = str(get_sample_path("mbf_align/chipseq_chr22.bam"))
        k = 17_000_000  # first 100k are empty...
        intervals = [
            ("chr22", 0, k, False),
            ("chr22", -10, k, False),
        ]
        forward = calculate_coverage(input, None, intervals)
        assert len(forward) == 2
        assert len(forward[0]) == k
        assert len(forward[1]) == k + 10
        assert np.sum(forward[0]) > 0
        assert np.sum(forward[1]) > 0
        assert np.sum(forward[0]) == np.sum(forward[1])
