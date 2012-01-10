add_executable( vtk-tensor-glyph vtkTensorGlyphTest.cxx )
target_link_libraries( vtk-tensor-glyph ${VTK_LIBRARIES} )
add_test( vtkTensorGlyphTest vtk-tensor-glyph strain_flavors_strain.vtk )

add_executable( vtk-sphere vtkSphereTest.cxx )
target_link_libraries( vtk-sphere ${VTK_LIBRARIES} )
add_test( vtkSphereTest vtk-sphere )

add_executable( vtk-one-sheeted-hyperboloid vtkOneSheetedHyperboloidTest.cxx )
target_link_libraries( vtk-one-sheeted-hyperboloid ${VTK_LIBRARIES} vtkSignedTensor )
add_test( vtkOneSheetedHyperboloidTest vtk-one-sheeted-hyperboloid )

add_test( GenerateStrainFlavors ${CMAKE_CURRENT_SOURCE_DIR}/generate_strain_flavors.py )

set_tests_properties( vtkTensorGlyphTest PROPERTIES DEPENDS GenerateStrainFlavors )
