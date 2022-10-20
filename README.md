# mbf_bam

Fast, multi-core read counters and the like based on the BAM fileformat and rust-htslib.


Part of the mbf_* suite from https://github.com/IMTMarburg

build wheels using
'docker run --rm -v $(pwd):/io ghcr.io/pyo3/maturin build --release -f'
then 'twine upload target/wheels/*'


To use in nix, a flake that provides a function taking
nixpkgs, mach-nix and the version to report to nix is
provided.
