#include "vtkActor.h"
#include "vtkActor2D.h"
#include "vtkCamera.h"
#include "vtkCubeSource.h"
#include "vtkExtractVOI.h"
#include "vtkLabelPlacementMapper.h"
#include "vtkLookupTable.h"
#include "vtkInteractorStyleImage.h"
#include "vtkOneSheetedHyperboloidSource.h"
#include "vtkPointData.h"
#include "vtkPointSetToLabelHierarchy.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkStringArray.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkSignedEigenvalueTensorGlyph.h"
#include "vtkTwoSheetedHyperboloidSource.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include <iostream>
#include <fstream>

int main( int argc, char * argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " <tensor_image>.vtk" << std::endl;
    return 1;
    }

  vtkSmartPointer< vtkStructuredPointsReader > reader = vtkSmartPointer< vtkStructuredPointsReader >::New();
  reader->SetFileName( argv[1] );
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
  std::ifstream tensorsFile( "strain_flavors_eigenvalues.txt", std::ifstream::in );

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
  std::ifstream anglesFile( "strain_flavors_angles.txt", std::ifstream::in );
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
  pointSetToLabelHierarchy->SetInput( labelPoints );
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
  renderWindowInteractor->Start();

  return 0;
}
