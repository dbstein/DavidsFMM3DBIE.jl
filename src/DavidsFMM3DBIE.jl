# In case you want to know, why the last line of the docstring below looks like it is:
# It will show the package (local) path when help on the package is invoked like     help?> DavidsFMM3DBIE
# but it will interpolate to an empty string on CI server, 
# preventing appearing the server local path in the documentation built there.

"""
    Package DavidsFMM3DBIE v$(pkgversion(DavidsFMM3DBIE))

Experiment with wrapping FMM3DBIE and making installation easy.

$(isnothing(get(ENV, "CI", nothing)) ? ("\n" * "Package local path: " * pathof(DavidsFMM3DBIE)) : "") 
"""

module DavidsFMM3DBIE

# Write your package code here.

end
