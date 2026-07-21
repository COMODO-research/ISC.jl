"""
    iscdir()

# Description 

This function simply returns the string for the ISC package path. This is helpful for
instance to load items, such as meshes, from the `assets` folder. 
"""
function iscdir()
    pkgdir(@__MODULE__)
end
