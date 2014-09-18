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
#include "vtkExtractVOI.h"
#include "vtkLabelPlacementMapper.h"
#include "vtkLookupTable.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPointData.h"
#include "vtkPointSetToLabelHierarchy.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStringArray.h"
#include "vtkSphereSource.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkTensorGlyph.h"
#include "vtkTesting.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkNew.h"

#include <fstream>

int TestTensorGlyph( int argc, char * argv[] )
{
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0]
              << " <TensorImage.vtk> <Eigenvalues.txt> <Angles.txt>"
              << " -T TemporaryDirectory [-I]" << std::endl;
    return EXIT_FAILURE;
    }
  const char * tensorImageFileName = argv[1];
  const char * eigenvaluesFileName = argv[2];
  const char * anglesFileName = argv[3];

  vtkNew<vtkTesting> testHelper;
  testHelper->AddArguments(argc, const_cast<const char **>(argv));

  vtkSmartPointer< vtkStructuredPointsReader > reader = vtkSmartPointer< vtkStructuredPointsReader >::New();
  reader->SetFileName( tensorImageFileName );
  reader->Update();

  vtkSmartPointer< vtkStructuredPoints > structuredPoints = reader->GetOutput();
  const int * dimensions = structuredPoints->GetDimensions();

  vtkSmartPointer< vtkSphereSource > sphere = vtkSmartPointer< vtkSphereSource >::New();
  sphere->SetThetaResolution( 50 );
  sphere->SetPhiResolution( 50 );
  vtkSmartPointer< vtkTensorGlyph > tensorGlyph = vtkSmartPointer< vtkTensorGlyph >::New();
  tensorGlyph->SetSourceConnection( sphere->GetOutputPort() );
  tensorGlyph->SetInputConnection( reader->GetOutputPort() );
  tensorGlyph->SetColorModeToEigenvalues();
  tensorGlyph->SetExtractEigenvalues( true );
  tensorGlyph->SetScaleFactor( 0.7 );
  tensorGlyph->SetScaling( true );

  // this is needed to prevent some of the glyphs going black.
  // taken from Graphics/Testing/Python/TestTensorGlyph.py
  vtkSmartPointer< vtkPolyDataNormals > normals = vtkSmartPointer< vtkPolyDataNormals >::New();
  normals->SetInputConnection( tensorGlyph->GetOutputPort() );

  vtkSmartPointer< vtkLookupTable > lut = vtkSmartPointer< vtkLookupTable >::New();
  // lut.SetHueRange( 0.1, 0.667 )
  lut->SetRange( 0.0, 0.11 );
  lut->Build();

  vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  mapper->SetInputConnection( normals->GetOutputPort() );
  mapper->SetLookupTable( lut );
  mapper->SetColorModeToMapScalars();
  mapper->SetScalarRange( 0.0, 10.11 );
  mapper->SetScalarModeToUsePointData();

  vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
  actor->SetMapper( mapper );
  //actor->RotateZ( -90.0 );

  const unsigned int numberOfTensors = 8;
  const unsigned int numberOfAngles = 6;
  vtkSmartPointer< vtkPolyData > labelPoints = vtkSmartPointer< vtkPolyData >::New();
  vtkSmartPointer< vtkPoints > labelPointsPts = vtkSmartPointer< vtkPoints >::New();
  labelPointsPts->Allocate( numberOfAngles + numberOfTensors  + 2 );
  labelPoints->SetPoints( labelPointsPts );

  vtkSmartPointer< vtkStringArray > labels = vtkSmartPointer< vtkStringArray >::New();
  labels->SetNumberOfValues( numberOfAngles + numberOfTensors + 2 );
  labels->SetName( "labels" );
  std::ifstream tensorsFile( eigenvaluesFileName, std::ifstream::in );

  if( !tensorsFile.is_open() )
    {
    std::cerr << "Could not open the eigenvalues file." << std::endl;
    return 1;
    }
  std::string line;
  double location[3];
  location[0] = -2.0;
  location[1] = 0.0;
  location[2] = 0.0;
  for( unsigned int ii = 0; ii < numberOfTensors; ++ii )
    {
    std::getline(tensorsFile, line);
    labels->SetValue( ii, line.c_str() );
    labelPointsPts->InsertNextPoint( location );
    location[1] += 1.0;
    }
  tensorsFile.close();
  std::ifstream anglesFile( anglesFileName, std::ifstream::in );
  if( !anglesFile.is_open() )
    {
    std::cerr << "Could not open the angles file." << std::endl;
    return 1;
    }
  location[0] = 0.0;
  location[1] = -1.0;
  location[2] = 0.0;
  for( unsigned int ii = numberOfTensors; ii < numberOfAngles + numberOfTensors; ++ii )
    {
    std::getline(anglesFile, line);
    labels->SetValue( ii, line.c_str() );
    labelPointsPts->InsertNextPoint( location );
    location[0] += 1.0;
    }
  anglesFile.close();
  labelPoints->GetPointData()->AddArray( labels );

  location[0] = -3.3;
  location[1] = numberOfTensors / 2.0;
  location[2] = 0.0;
  labelPointsPts->InsertNextPoint( location );
  labels->SetValue( numberOfAngles + numberOfTensors, "Eigenvalues" );
  location[0] = numberOfAngles / 2.0;
  location[1] = -2.0;
  location[2] = 0.0;
  labelPointsPts->InsertNextPoint( location );
  labels->SetValue( numberOfAngles + numberOfTensors + 1, "Angles (degrees)" );

  vtkSmartPointer< vtkPointSetToLabelHierarchy > pointSetToLabelHierarchy = vtkSmartPointer< vtkPointSetToLabelHierarchy >::New();
  pointSetToLabelHierarchy->SetInputData( labelPoints );
  pointSetToLabelHierarchy->SetLabelArrayName( "labels" );
  pointSetToLabelHierarchy->Update();

  vtkSmartPointer< vtkLabelPlacementMapper > labelMapper = vtkSmartPointer< vtkLabelPlacementMapper >::New();
  labelMapper->SetInputConnection( pointSetToLabelHierarchy->GetOutputPort() );
  vtkSmartPointer< vtkActor2D > labelActor = vtkSmartPointer< vtkActor2D >::New();
  labelActor->SetMapper( labelMapper );

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground( 0.2, 0.2, 0.2 );
  //renderer->SetActiveCamera( camera );
  renderer->AddActor( actor );
  renderer->AddActor( labelActor );

  vtkSmartPointer< vtkRenderWindow > renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
  renderWindow->AddRenderer( renderer );
  renderWindow->SetSize( 800, 800 );

  vtkSmartPointer< vtkInteractorStyleImage > interactorStyle = vtkSmartPointer< vtkInteractorStyleImage >::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

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
