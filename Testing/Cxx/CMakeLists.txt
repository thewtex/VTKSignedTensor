set(TestTensorGlyph_ARGS DATA{../Data/Input/StrainFlavorsStrain.vtk}
  DATA{../Data/Input/StrainFlavorsEigenvalues.txt}
  DATA{../Data/Input/StrainFlavorsAngles.txt})
set(TestSignedEigenvalueTensorGlyph_ARGS DATA{../Data/Input/StrainFlavorsStrain.vtk}
  DATA{../Data/Input/StrainFlavorsEigenvalues.txt}
  DATA{../Data/Input/StrainFlavorsAngles.txt})
set(TestUniAxialLoading_ARGS DATA{../Data/Input/UniaxialCompression.vtk} 1)
vtk_add_test_cxx(${vtk-module}CxxTests tests
  NO_DATA
  TestTensorGlyph.cxx
  TestSignedEigenvalueTensorGlyph.cxx
  TestSphere.cxx
  TestOneSheetedHyperboloid.cxx
  TestTwoSheetedHyperboloid.cxx
  TestUniAxialLoading.cxx
  )
vtk_test_cxx_executable(${vtk-module}CxxTests tests)
