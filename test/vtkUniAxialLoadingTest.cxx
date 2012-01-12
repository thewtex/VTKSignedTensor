#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkCubeSource.h"
#include "vtkOneSheetedHyperboloidSource.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkSignedEigenvalueTensorGlyph.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkTensorGlyph.h"
#include "vtkTwoSheetedHyperboloidSource.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include <iostream>

void createGlyphs( const char * filename, vtkSmartPointer< vtkPolyData > & signedGlyphs )
{
  vtkSmartPointer< vtkStructuredPointsReader > reader = vtkSmartPointer< vtkStructuredPointsReader >::New();
  reader->SetFileName( filename );
  reader->Update();

  vtkSmartPointer< vtkStructuredPoints > structuredPoints = reader->GetOutput();
  const int * dimensions = structuredPoints->GetDimensions();

  const int resolution = 14;
  vtkSmartPointer< vtkSphereSource > sphere = vtkSmartPointer< vtkSphereSource >::New();
  sphere->SetThetaResolution( resolution );
  sphere->SetPhiResolution( resolution );

  vtkSmartPointer< vtkOneSheetedHyperboloidSource > oneSheetedHyperboloid = vtkSmartPointer< vtkOneSheetedHyperboloidSource >::New();
  oneSheetedHyperboloid->SetThetaResolution( resolution );
  oneSheetedHyperboloid->SetZResolution( resolution );
  vtkSmartPointer< vtkPolyDataNormals > oneSheetedHyperboloidWithNormals = vtkSmartPointer< vtkPolyDataNormals >::New();
  oneSheetedHyperboloidWithNormals->SetInputConnection( oneSheetedHyperboloid->GetOutputPort() );

  vtkSmartPointer< vtkTwoSheetedHyperboloidSource > twoSheetedHyperboloid = vtkSmartPointer< vtkTwoSheetedHyperboloidSource >::New();
  twoSheetedHyperboloid->SetThetaResolution( resolution );
  twoSheetedHyperboloid->SetZResolution( resolution );
  vtkSmartPointer< vtkPolyDataNormals > twoSheetedHyperboloidWithNormals = vtkSmartPointer< vtkPolyDataNormals >::New();
  twoSheetedHyperboloidWithNormals->SetInputConnection( twoSheetedHyperboloid->GetOutputPort() );

  vtkSmartPointer< vtkCubeSource > cube = vtkSmartPointer< vtkCubeSource >::New();

  const double glyphScale = 15.0;

  vtkSmartPointer< vtkSignedEigenvalueTensorGlyph > tensorGlyph = vtkSmartPointer< vtkSignedEigenvalueTensorGlyph >::New();
  tensorGlyph->SetInputConnection(1, sphere->GetOutputPort() );
  tensorGlyph->SetInputConnection(2, oneSheetedHyperboloidWithNormals->GetOutputPort() );
  tensorGlyph->SetInputConnection(3, twoSheetedHyperboloidWithNormals->GetOutputPort() );
  tensorGlyph->SetInputConnection(4, cube->GetOutputPort() );
  tensorGlyph->SetInputConnection( reader->GetOutputPort() );
  tensorGlyph->SetColorModeToEigenvalues();
  tensorGlyph->SetExtractEigenvalues( true );
  tensorGlyph->SetScaleFactor( glyphScale );
  tensorGlyph->SetScaling( true );

  // this is needed to prevent some of the glyphs going black.
  // taken from Graphics/Testing/Python/TestTensorGlyph.py
  vtkSmartPointer< vtkPolyDataNormals > normals = vtkSmartPointer< vtkPolyDataNormals >::New();
  normals->SetInputConnection( tensorGlyph->GetOutputPort() );

  normals->Update();

  signedGlyphs = normals->GetOutput();
}

int main( int argc, char * argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " uniaxial_tension.vtk uniaxial_compression.vtk" << std::endl;
    return 1;
    }

  vtkSmartPointer< vtkPolyData > signedGlyphs;
  createGlyphs( argv[1], signedGlyphs );

  vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  mapper->SetInput( signedGlyphs );

  vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
  actor->SetMapper( mapper );

  vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();
  outline->SetInput( signedGlyphs );
  vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  outlineMapper->SetInputConnection( outline->GetOutputPort() );
  vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
  outlineActor->SetMapper( outlineMapper );

  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground( 0.2, 0.2, 0.2 );
  renderer->AddActor( actor );
  renderer->AddActor( outlineActor );
  renderer->AddActor( axes );

  vtkSmartPointer< vtkRenderWindow > renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
  renderWindow->AddRenderer( renderer );
  renderWindow->SetSize( 800, 800 );

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  renderWindow->Render();
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  return 0;
}

