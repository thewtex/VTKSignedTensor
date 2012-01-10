cmake_minimum_required( VERSION 2.8 )

project( VTK_STRAIN_GLYPH )

find_package( VTK REQUIRED )
include( ${VTK_USE_FILE} )

add_executable( vtk-tensor-glyph vtk_tensor_glyph.cxx )
target_link_libraries( vtk-tensor-glyph ${VTK_LIBRARIES} )

add_executable( vtk-sphere vtk_sphere.cxx )
target_link_libraries( vtk-sphere ${VTK_LIBRARIES} )
