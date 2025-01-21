## New package

This is a new package.

## Test environments

* Mac aarch64-apple-darwin20 (64-bit), R 4.4.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Tests

All tests pass.

## Previous submission (0.1.0)

Possibly misspelled words in DESCRIPTION:
  Hemodynmaic (3:8)
  hemodynamic (34:27)
  hrf (34:58)

* The misspelling has been corrected. "hemodynamic" and "hrf" have been added to the WORDLIST.

Package has a VignetteBuilder field but no prebuilt vignette index.

* A vignette has been added.

The Title field should be in title case. Current version is:
  'Hemodynmaic response function'
In title case that is:
  'Hemodynmaic Response Function'

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
Check: S3 generic/method consistency, Result: NOTE
  Mismatches for apparent methods not registered:
  is.matrix:
    function(x)
  is.matrix.or.df:
    function(q)
  See section 'Registering S3 methods' in the 'Writing R Extensions'
  manual.

* This function has been renamed to is_matrix_or_df
