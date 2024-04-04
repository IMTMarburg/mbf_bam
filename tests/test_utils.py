from mbf_bam import (
    reheader_and_rename_chromosomes,
    job_reheader_and_rename_chromosomes,
    job_filter_and_rename,
    fix_sorting_to_be_deterministic,
)
from pathlib import Path
import pysam
import pypipegraph as ppg
import shutil
import pytest


def get_sample_path(name):
    return Path(name.replace("mbf_align", "../../../sample_data"))


class TestReheader:
    def test_rename(self, new_pipegraph):
        ppg.util.global_pipegraph.quiet = False
        input = get_sample_path("mbf_align/ex2.bam")
        output = "out.bam"
        job_reheader_and_rename_chromosomes(
            input, output, {"chr1": "shu", "chr2": "sha"}
        )
        ppg.run_pipegraph()
        assert Path("out.bam").exists()
        f = pysam.Samfile("out.bam")
        assert set(f.references) == set(["shu", "sha"])

    def test_rename_raises_on_no_replacement(self, new_pipegraph):
        ppg.util.global_pipegraph.quiet = False
        input = get_sample_path("mbf_align/ex2.bam")
        output = "out.bam"
        j = job_reheader_and_rename_chromosomes(input, output, {})
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert not Path("out.bam").exists()
        assert "No replacement happened" in str(j.exception)


class TestSubtract:
    def test_subtract_subset(self, new_pipegraph):
        from mbf_bam import subtract_bam

        input = get_sample_path("mbf_align/chipseq_chr22.bam")
        minued = get_sample_path("mbf_align/chipseq_chr22_subset_plus_unmapped.bam")
        output = "output.bam"
        print(input, input.exists())
        print(minued, minued.exists())
        subtract_bam(str(output), str(input.absolute()), str(minued.absolute()))
        f = pysam.Samfile(output)
        should = 80495
        total = sum((x.total for x in f.get_index_statistics()))
        assert should == total


class TestFilterAnd_Rename:
    def test_filter_and_rename_ommited(self, new_pipegraph):
        ppg.util.global_pipegraph.quiet = False
        input = get_sample_path("mbf_align/ex2.bam")
        output = "out.bam"
        job_filter_and_rename(input, output, {"chr2": "sha"})  # "chr1": None,
        ppg.run_pipegraph()
        assert Path("out.bam").exists()
        f = pysam.Samfile("out.bam")
        assert set(f.references) == set(["sha"])
        assert len(list(f.fetch("sha"))) == 7

    def test_filter_and_rename_None(self, new_pipegraph):
        ppg.util.global_pipegraph.quiet = False
        input = get_sample_path("mbf_align/ex2.bam")
        output = "out.bam"
        job_filter_and_rename(input, output, {"chr1": None, "chr2": "sha"})
        ppg.run_pipegraph()
        assert Path("out.bam").exists()
        f = pysam.Samfile("out.bam")
        assert set(f.references) == set(["sha"])
        assert len(list(f.fetch("sha"))) == 7


class TestFixSort:
    def test_semi_sorted(self, new_pipegraph):
        input = get_sample_path("mbf_align/subread_semi_sorted.bam")
        output = "output.bam"
        fix_sorting_to_be_deterministic(input, output)
        assert Path(output).exists()
        pysam.index(output)  # that must work...
        last = None

        f = pysam.Samfile(input)
        violations = 0
        for read in f.fetch("1"):
            if last:
                if not (
                    (last.pos < read.pos)
                    or ((last.pos == read.pos) and (last.query_name < read.query_name))
                ):
                    violations += 1
            last = read
        assert violations > 0  # so that we are actually testing something.

        f = pysam.Samfile(output)
        last = None
        for read in f.fetch("1"):
            print(read.pos)
            if last:
                assert (last.pos < read.pos) or (
                    (last.pos == read.pos) and (last.query_name < read.query_name)
                )
            last = read

    def test_semi_sorted_but_no_index(self, new_pipegraph):
        real_input = get_sample_path("mbf_align/subread_semi_sorted.bam")
        input = "input.bam"
        shutil.copy(real_input, input)
        output = "output.bam"
        fix_sorting_to_be_deterministic(input, output)
        assert Path(output).exists()
        pysam.index(output)  # that must work...

        f = pysam.Samfile(output)
        last = None
        for read in f.fetch("1"):
            print(read.pos)
            if last:
                assert (last.pos < read.pos) or (
                    (last.pos == read.pos) and (last.query_name < read.query_name)
                )
            last = read

    def test_spliced(self, new_pipegraph):
        input = get_sample_path("mbf_align/spliced_reads.bam")
        output = "output.bam"
        fix_sorting_to_be_deterministic(input, output)
        assert Path(output).exists()
        last = None

        f = pysam.Samfile(input)
        violations = 0
        for read in f.fetch(until_eof=True):
            if last:
                if last.tid != read.tid:
                    if not (last.pos < read.pos) or (
                        (last.pos == read.pos) and (last.query_name < read.query_name)
                    ):
                        violations += 1
            last = read
        # assert violations > 0  # so that we are actually testing something.

        f = pysam.Samfile(output)
        for read in f.fetch(until_eof=True):
            if last:
                if last.tid != read.tid:
                    assert (last.pos < read.pos) or (
                        (last.pos == read.pos) and (last.query_name < read.query_name)
                    )
            last = read
        pysam.index(output)  # that must work...

    def test_unsorted_raises(self):
        input = get_sample_path("mbf_align/unsorted.bam")
        output = "output.bam"
        with pytest.raises(ValueError):
            fix_sorting_to_be_deterministic(input, output)
        assert not Path(output).exists()
