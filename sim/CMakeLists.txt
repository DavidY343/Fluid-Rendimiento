# Add to this list all files related to sim library
add_library(sim
        grid.cpp
        grid.hpp
        block.cpp
        block.hpp
        progargs.hpp
        progargs.cpp
)
# Use this line only if you have dependencies from stim to GSL
# target_link_libraries (sim PRIVATE Microsoft.GSL::GSL)

