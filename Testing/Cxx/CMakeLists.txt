add_executable( vtk-tensor-glyph vtkTensorGlyphTest.cxx )
target_link_libraries( vtk-tensor-glyph ${VTK_LIBRARIES} )
add_test( vtkTensorGlyphTest vtk-tensor-glyph strain_flavors_strain.vtk )

add_executable( vtk-sphere vtkSphereTest.cxx )
target_link_libraries( vtk-sphere ${VTK_LIBRARIES} )
add_test( vtkSphereTest vtk-sphere )

add_executable( vtk-one-sheeted-hyperboloid vtkOneSheetedHyperboloidTest.cxx )
target_link_libraries( vtk-one-sheeted-hyperboloid ${VTK_LIBRARIES} vtkSignedTensor )
add_test( vtkOneSheetedHyperboloidTest vtk-one-sheeted-hyperboloid )

add_executable( vtk-two-sheeted-hyperboloid vtkTwoSheetedHyperboloidTest.cxx )
target_link_libraries( vtk-two-sheeted-hyperboloid ${VTK_LIBRARIES} vtkSignedTensor )
add_test( vtkTwoSheetedHyperboloidTest vtk-two-sheeted-hyperboloid )

add_executable( vtk-signed-tensor-glyph vtkSignedEigenvalueTensorGlyphTest.cxx )
target_link_libraries( vtk-signed-tensor-glyph ${VTK_LIBRARIES} vtkSignedTensor )
add_test( vtkSignedEigenvalueTensorGlyphTest vtk-signed-tensor-glyph
  strain_flavors_strain.vtk )

add_executable( vtk-uniaxial-load vtkUniAxialLoadingTest.cxx )
target_link_libraries( vtk-uniaxial-load ${VTK_LIBRARIES} vtkSignedTensor )
add_test( vtkUniAxialLoadingTensionTest vtk-uniaxial-load uniaxial_tension.vtk 0 )
add_test( vtkUniAxialLoadingCompressionTest vtk-uniaxial-load uniaxial_compression.vtk 1 )

add_test( GenerateStrainFlavors ${CMAKE_CURRENT_SOURCE_DIR}/generate_strain_flavors.py )
add_test( GenerateUniaxialLoad ${CMAKE_CURRENT_SOURCE_DIR}/generate_uniaxial_load.py )

set_tests_properties( vtkTensorGlyphTest vtkSignedEigenvalueTensorGlyphTest
  PROPERTIES DEPENDS GenerateStrainFlavors )
set_tests_properties( vtkUniAxialLoadingTensionTest vtkUniAxialLoadingCompressionTest
  PROPERTIES DEPENDS GenerateUniaxialLoad )

vtk_add_test_cxx(${vtk-module}CxxTests tests
  )
vtk_test_cxx_executable(${vtk-module}CxxTests tests)
