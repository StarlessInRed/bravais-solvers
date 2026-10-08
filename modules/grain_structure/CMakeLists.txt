add_library(prisms_grain_structure STATIC
  upstream/src/grains/FloodFiller.cc
  upstream/src/grains/SimplifiedGrainRepresentation.cc
  upstream/src/grains/OrderParameterRemapper.cc
  native.cpp)
target_compile_features(prisms_grain_structure PUBLIC cxx_std_20)
target_include_directories(prisms_grain_structure PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_SOURCE_DIR}/upstream/include)
target_link_libraries(prisms_grain_structure PUBLIC PRISMS-PF::prisms_pf)
