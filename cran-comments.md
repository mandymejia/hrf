## New package

This is a new package.

## Test environments

* Mac aarch64-apple-darwin20 (64-bit), R 4.4.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Tests

All tests pass.

## Previous submission (0.1.1)

The Description field is intended to be a (one paragraph) description of
what the package does and why it may be useful. Please add more details
about the package functionality and implemented methods in your
Description text.
For more details:
<https://contributor.r-project.org/cran-cookbook/general_issues.html#description-length>

* A fuller Description has been added.

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'https:' and angle brackets for
auto-linking. (If you want to add a title as well please put it in
quotes: "Title")
For more details:
<https://contributor.r-project.org/cran-cookbook/description_issues.html#references>

* There are no references necessary (yet).

You write information messages to the console that cannot be easily
suppressed.
It is more R like to generate objects that can be used to extract the
information a user is interested in, and then print() that object.
Instead of print()/cat() rather use message()/warning() or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console. (except for print, summary, interactive functions)
-> R/check_data_inputs.R; R/GLM_multi.R; R/make_design.methods.R;
R/make_design.R; R/multiGLM.R; R/util.R
For more details:
<https://contributor.r-project.org/cran-cookbook/code_issues.html#using-printcat>

* `print()` and `cat()` messages not inside an `if(verbose){...}` statement have been replaced with `message` or `warning`, or wrapped in an `if(verbose){...}` statement.