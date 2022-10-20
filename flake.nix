{
  description = "Wraps mbf-bam into an mach-nix importable builder";

  inputs = {
    import-cargo.url = github:edolstra/import-cargo;
  };

  outputs = {
    self,
    import-cargo,
  }: let
    inherit (import-cargo.builders) importCargo;
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
  };
}
