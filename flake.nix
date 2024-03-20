{
  description = "Wraps mbf-bam into an mach-nix importable builder";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/23.11";
    import-cargo.url = github:edolstra/import-cargo;
  };

  outputs = {
    self,
    nixpkgs,
    import-cargo,
  }: let
    inherit (import-cargo.builders) importCargo;
    system = "x86_64-linux";
    npkgs = import nixpkgs {inherit system;};
  in let
    build_mbf_bam = pkgs: pythonpkgs: outside_version: let
      cargo_in = importCargo {
        lockFile = ./Cargo.lock;
        inherit pkgs;
      };
    in
      pythonpkgs.buildPythonPackage
      {
        src = ./.;
        version = outside_version;

        nativeBuildInputs = [
          cargo_in.cargoHome

          # Build-time dependencies
          pkgs.rustc
          pkgs.cargo
          pkgs.openssl.dev
          pkgs.perl
        ];
        requirementsExtra = ''
          setuptools-rust
        '';
      };
  in {
    # pass in nixpkgs, mach-nix and what you want it to report back as a version
    mach-nix-build-python-package = build_mbf_bam;

    devShell.x86_64-linux = npkgs.mkShell {
      # be sure to set this back in your build scripts,
      # otherwise pyo3 will get recompiled all the time
      CARGO_TARGET_DIR = "target_rust_analyzer";

      nativeBuildInputs = [
        npkgs.rustc
        npkgs.cargo
        npkgs.cargo-binutils
        npkgs.rust-analyzer
        npkgs.git
        npkgs.cargo-udeps
        npkgs.cargo-audit
        npkgs.cargo-vet
        npkgs.cargo-outdated
        npkgs.bacon
        npkgs.maturin
        (npkgs.python3.withPackages (p: [p.pytest p.pytest-cov p.pandas p.tomlkit]))
      ];
    };
  };
}
