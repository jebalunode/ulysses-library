cmake_minimum_required(VERSION 3.15)

set(GFN2_XTB_SOURCES
    gfn2-xtb_opt_lib.cpp
    gfn2-xtb_spe_lib.cpp
)

# Static library
add_library(gfn2_xtb_static STATIC ${GFN2_XTB_SOURCES})
set_target_properties(gfn2_xtb_static PROPERTIES OUTPUT_NAME gfn2_xtb)

# Shared library
# add_library(gfn2_xtb_shared SHARED ${GFN2_XTB_SOURCES})
#set_target_properties(gfn2_xtb_shared PROPERTIES OUTPUT_NAME gfn2_xtb)

# Compiler optimizations
target_compile_options(gfn2_xtb_static PRIVATE -O3 -march=native -O0 -fopenmp)
#target_compile_options(gfn2_xtb_shared PRIVATE -O3 -march=native -O0 -fopenmp)

# OpenMP support
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(gfn2_xtb_static PUBLIC OpenMP::OpenMP_CXX)
    #    target_link_libraries(gfn2_xtb_shared PUBLIC OpenMP::OpenMP_CXX)
endif()

# === Install rules ===
install(TARGETS gfn2_xtb_static 
        ARCHIVE DESTINATION lib    # .a
        LIBRARY DESTINATION lib    # .so
)

