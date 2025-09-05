cmake_minimum_required(VERSION 3.15)

# Option: choose static or shared library
option(USE_SHARED_LIB "Link executables against shared library" OFF)

if(USE_SHARED_LIB)
    set(GFN2_XTB_LIB gfn2_xtb_shared)
else()
    set(GFN2_XTB_LIB gfn2_xtb_static)
endif()

# Executable: gfn2_xtb_spe
add_executable(gfn2_xtb_spe main1.cpp)
target_compile_options(gfn2_xtb_spe PRIVATE -O3 -march=native -fopenmp)
target_link_libraries(gfn2_xtb_spe PRIVATE ${GFN2_XTB_LIB})

# Executable: gfn2_xtb_opt
add_executable(gfn2_xtb_opt main2.cpp)
target_compile_options(gfn2_xtb_opt PRIVATE -O3 -march=native -fopenmp)
target_link_libraries(gfn2_xtb_opt PRIVATE ${GFN2_XTB_LIB})

# === Install rules ===
install(TARGETS gfn2_xtb_spe gfn2_xtb_opt
        RUNTIME DESTINATION bin
)

