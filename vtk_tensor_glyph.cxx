#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkExtractVOI.h"
#include "vtkLookupTable.h"
#include "vtkInteractorStyleImage.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkTensorGlyph.h"
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
