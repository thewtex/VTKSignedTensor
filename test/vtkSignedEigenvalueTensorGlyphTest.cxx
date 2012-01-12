#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCubeSource.h"
#include "vtkExtractVOI.h"
#include "vtkLookupTable.h"
#include "vtkInteractorStyleImage.h"
#include "vtkOutlineFilter.h"
#include "vtkOneSheetedHyperboloidSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkSignedEigenvalueTensorGlyph.h"
#include "vtkTwoSheetedHyperboloidSource.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include <iostream>

int main( int argc, char * argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " <tensor_image>.vtk" << std::endl;
    return 1;
    }

  //vtkSmartPointer< vtkCamera > camera = vtkSmartPointer< vtkCamera >::New();
  //camera->SetViewUp     (0.0, 1.0, 0.0);
  //camera->SetFocalPoint (404.70907910888917, -202.13628697726463, 0.0);
  //camera->SetPosition   (404.70907910888917, -202.13628697726463, 850.1861762830398);

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
  tensorGlyph->SetSourceConnection(0, sphere->GetOutputPort() );
  tensorGlyph->SetSourceConnection(1, oneSheetedHyperboloidWithNormals->GetOutputPort() );
  tensorGlyph->SetSourceConnection(2, twoSheetedHyperboloidWithNormals->GetOutputPort() );
  tensorGlyph->SetSourceConnection(3, cube->GetOutputPort() );
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

  vtkSmartPointer< vtkOutlineFilter > outline = vtkSmartPointer< vtkOutlineFilter >::New();
  outline->SetInputConnection( reader->GetOutputPort() );
  vtkSmartPointer< vtkPolyDataMapper > outlineMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  outlineMapper->SetInputConnection( outline->GetOutputPort() );
  vtkSmartPointer< vtkActor > outlineActor = vtkSmartPointer< vtkActor >::New();
  outlineActor->SetMapper( outlineMapper );
  //outlineActor->RotateZ( -90.0 );

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground( 0.2, 0.2, 0.2 );
  //renderer->SetActiveCamera( camera );
  renderer->AddActor( actor );
  renderer->AddActor( outlineActor );

  vtkSmartPointer< vtkRenderWindow > renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
  renderWindow->AddRenderer( renderer );
  renderWindow->SetSize( 600, 600 );

  vtkSmartPointer< vtkInteractorStyleImage > interactorStyle = vtkSmartPointer< vtkInteractorStyleImage >::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  renderWindow->Render();
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  return 0;
}
