set(Module_SRCS
  vtkOneSheetedHyperboloidSource.cxx
  vtkTwoSheetedHyperboloidSource.cxx
  vtkSignedEigenvalueTensorGlyph.cxx
  )

vtk_module_library(SignedTensor ${Module_SRCS})
