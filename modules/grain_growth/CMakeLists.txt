add_library(grain_growth_model STATIC model.cpp)
target_compile_features(grain_growth_model PUBLIC cxx_std_20)
target_include_directories(grain_growth_model PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_library(grain_growth_prismspf STATIC prismspf/equation.cpp)
target_compile_features(grain_growth_prismspf PUBLIC cxx_std_20)
target_link_libraries(grain_growth_prismspf PUBLIC grain_growth_model bravais_grain_representation PRISMS-PF::prisms_pf)
