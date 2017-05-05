/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestContourTriangulator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This example demonstrates how to use vtkTensorGlyph
//
// The command line arguments are:
// -I        => run in interactive mode; unless this is used, the program will
//              not allow interaction and exit
#include "vtkActor.h"
#include "vtkActor2D.h"
#include "vtkCamera.h"
#include "vtkCubeSource.h"
#include "vtkExtractVOI.h"
#include "vtkLookupTable.h"
#include "vtkNew.h"
#include "vtkInteractorStyleImage.h"
#include "vtkOneSheetedHyperboloidSource.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkSignedEigenvalueTensorGlyph.h"
#include "vtkTesting.h"
#include "vtkTwoSheetedHyperboloidSource.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

int TestSignedEigenvalueTensorGlyph( int argc, char * argv[] )
{
  if ( argc < 2 )
  {
    std::cerr << "Usage: " << argv[0]
              << " <TensorImage.vtk>"
              << " -T TemporaryDirectory [-I]" << std::endl;
    return EXIT_FAILURE;
  }
  const char * tensorImageFileName = argv[1];

  vtkSmartPointer< vtkStructuredPointsReader > reader = vtkSmartPointer< vtkStructuredPointsReader >::New();
  reader->SetFileName( tensorImageFileName );
  reader->Update();

  vtkSmartPointer< vtkStructuredPoints > structuredPoints = reader->GetOutput();
  const int * dimensions = structuredPoints->GetDimensions();

  vtkSmartPointer< vtkSphereSource > sphere = vtkSmartPointer< vtkSphereSource >::New();
  sphere->SetThetaResolution( 20 );
  sphere->SetPhiResolution( 20 );

  vtkSmartPointer< vtkOneSheetedHyperboloidSource > oneSheetedHyperboloid = vtkSmartPointer< vtkOneSheetedHyperboloidSource >::New();
  oneSheetedHyperboloid->SetThetaResolution( 20 );
  oneSheetedHyperboloid->SetZResolution( 20 );
  vtkSmartPointer< vtkPolyDataNormals > oneSheetedHyperboloidWithNormals = vtkSmartPointer< vtkPolyDataNormals >::New();
  oneSheetedHyperboloidWithNormals->SetInputConnection( oneSheetedHyperboloid->GetOutputPort() );

  vtkSmartPointer< vtkTwoSheetedHyperboloidSource > twoSheetedHyperboloid = vtkSmartPointer< vtkTwoSheetedHyperboloidSource >::New();
  twoSheetedHyperboloid->SetThetaResolution( 20 );
  twoSheetedHyperboloid->SetZResolution( 20 );
  vtkSmartPointer< vtkPolyDataNormals > twoSheetedHyperboloidWithNormals = vtkSmartPointer< vtkPolyDataNormals >::New();
  twoSheetedHyperboloidWithNormals->SetInputConnection( twoSheetedHyperboloid->GetOutputPort() );

  vtkSmartPointer< vtkCubeSource > cube = vtkSmartPointer< vtkCubeSource >::New();

  vtkSmartPointer< vtkSignedEigenvalueTensorGlyph > tensorGlyph = vtkSmartPointer< vtkSignedEigenvalueTensorGlyph >::New();
  tensorGlyph->SetInputConnection(1, sphere->GetOutputPort() );
  tensorGlyph->SetInputConnection(2, oneSheetedHyperboloidWithNormals->GetOutputPort() );
  tensorGlyph->SetInputConnection(3, twoSheetedHyperboloidWithNormals->GetOutputPort() );
  tensorGlyph->SetInputConnection(4, cube->GetOutputPort() );
  tensorGlyph->SetInputConnection( reader->GetOutputPort() );
  //tensorGlyph->SetColorModeToEigenvalues();
  tensorGlyph->SetExtractEigenvalues( true );
  tensorGlyph->SetScaleFactor( 0.7 );
  tensorGlyph->SetScaling( true );
  tensorGlyph->Print( std::cout );

  // this is needed to prevent some of the glyphs going black.
  // taken from Graphics/Testing/Python/TestTensorGlyph.py
  vtkSmartPointer< vtkPolyDataNormals > normals = vtkSmartPointer< vtkPolyDataNormals >::New();
  normals->SetInputConnection( tensorGlyph->GetOutputPort() );

  //vtkSmartPointer< vtkLookupTable > lut = vtkSmartPointer< vtkLookupTable >::New();
  // lut.SetHueRange( 0.1, 0.667 )
  //lut->Build();

  vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  mapper->SetInputConnection( normals->GetOutputPort() );
  //mapper->SetLookupTable( lut );
  //mapper->SetColorModeToMapScalars();
  //mapper->SetScalarRange( 0.0, 10.11 );
  //mapper->SetScalarModeToUsePointData();

  vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
  actor->SetMapper( mapper );
  actor->RotateZ( 2.0 );
  actor->RotateY( 10.0 );
  actor->RotateX( 10.0 );

  vtkProperty * property = actor->GetProperty();
  property->SetColor( 0.2, 0.7, 0.3 );
  property->SetAmbient( 0.1 );
  property->SetDiffuse( 0.9 );
  property->SetSpecular( 0.2 );
  property->SetSpecularPower( 5 );
  //property->EdgeVisibilityOn();
  property->SetEdgeColor( 0.1, 0.6, 0.1 );
  property->SetLineWidth( 0.5 );
  property->SetInterpolationToGouraud();
  property->ShadingOn();

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground( 0.2, 0.2, 0.2 );
  renderer->AddActor( actor );

  vtkSmartPointer< vtkRenderWindow > renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
  renderWindow->AddRenderer( renderer );
  renderWindow->SetSize( 600, 800 );

  vtkSmartPointer< vtkInteractorStyleImage > interactorStyle = vtkSmartPointer< vtkInteractorStyleImage >::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  vtkSmartPointer< vtkCamera > camera = vtkSmartPointer< vtkCamera >::New();
  camera->SetPosition( 0, 0, 15 );
  camera->SetFocalPoint( 2.5, 3.5, 0 );
  renderer->SetActiveCamera( camera );

  renderWindow->Render();

  renderWindowInteractor->Initialize();
  const int returnValue = vtkTesting::Test(argc, argv, renderWindow, 20);
  if( returnValue == vtkTesting::DO_INTERACTOR )
  {
    renderWindowInteractor->Start();
  }

  if ((returnValue == vtkTesting::PASSED) || (returnValue == vtkTesting::DO_INTERACTOR))
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}
