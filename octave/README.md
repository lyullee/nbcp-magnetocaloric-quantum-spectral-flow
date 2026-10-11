# GNU Octave cross-check

This directory contains a package-free Octave/MATLAB reproduction of the
bulk polarized-state stability calculation.

Run from the project root:

```matlab
addpath('octave');
bulk_polarized_stability_octave();
```

The script writes `results/bulk_polarized_stability_octave.json`.  The
expected tin-foil values from the independently tested Python implementation
are M0 `1.7743602002 T` and M4 `1.8676712357 T`.

The calculation was run successfully with Octave 11.3.0.  Smart App Control
initially blocked the unsigned official `libcurl-4.dll`.  The project-local
runtime now uses the same-version, identical-export, validly signed Git for
Windows `libcurl-4.dll` and its four signed dependencies.  The original DLL is
preserved as `libcurl-4.dll.octave-original`; Windows security policy was not
changed.  `tools/setup_octave_signed_shim.ps1` reproduces the local setup.

The Octave output agrees with the independent Python M0--M4 critical fields to
at worst `2.22e-16 T`.
