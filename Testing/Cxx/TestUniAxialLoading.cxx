#include "vtkActor.h"
#include "vtkArrowSource.h"
#include "vtkAxesActor.h"
#include "vtkCamera.h"
#include "vtkCubeSource.h"
#include "vtkOneSheetedHyperboloidSource.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkSignedEigenvalueTensorGlyph.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkTensorGlyph.h"
#include "vtkTwoSheetedHyperboloidSource.h"
#include "vtkTransform.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include <iostream>

void createGlyphs( const char * filename,
  vtkSmartPointer< vtkPolyData > & signedGlyphs,
  vtkSmartPointer< vtkPolyData > & glyphs )
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

  vtkSmartPointer< vtkSignedEigenvalueTensorGlyph > signedTensorGlyph = vtkSmartPointer< vtkSignedEigenvalueTensorGlyph >::New();
  signedTensorGlyph->SetInputConnection(1, sphere->GetOutputPort() );
  signedTensorGlyph->SetInputConnection(2, oneSheetedHyperboloidWithNormals->GetOutputPort() );
  signedTensorGlyph->SetInputConnection(3, twoSheetedHyperboloidWithNormals->GetOutputPort() );
  signedTensorGlyph->SetInputConnection(4, cube->GetOutputPort() );
  signedTensorGlyph->SetInputConnection( reader->GetOutputPort() );
  signedTensorGlyph->SetColorModeToEigenvalues();
  signedTensorGlyph->SetExtractEigenvalues( true );
  signedTensorGlyph->SetScaleFactor( glyphScale );
  signedTensorGlyph->SetScaling( true );

  // this is needed to prevent some of the glyphs going black.
  // taken from Graphics/Testing/Python/TestTensorGlyph.py
  vtkSmartPointer< vtkPolyDataNormals > signedNormals = vtkSmartPointer< vtkPolyDataNormals >::New();
  signedNormals->SetInputConnection( signedTensorGlyph->GetOutputPort() );
  signedNormals->Update();

  signedGlyphs = signedNormals->GetOutput();

  vtkSmartPointer< vtkTensorGlyph > tensorGlyph = vtkSmartPointer< vtkTensorGlyph >::New();
  tensorGlyph->SetSourceConnection( sphere->GetOutputPort() );
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

  glyphs = normals->GetOutput();
}

int main( int argc, char * argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " <uniaxial strain image>.vtk <tension (0) or compression (1)>" << std::endl;
    return 1;
    }

  vtkSmartPointer< vtkStructuredPointsReader > reader = vtkSmartPointer< vtkStructuredPointsReader >::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  vtkSmartPointer< vtkStructuredPoints > structuredPoints = reader->GetOutput();
  structuredPoints->ComputeBounds();
  double blockBounds[6];
  structuredPoints->GetBounds( blockBounds );

  bool tension = false;
  if( argv[2][0] == '0' )
    {
    tension = true;
    }

  std::cout << "block bounds: " << blockBounds[0] << " " << blockBounds[1] << std::endl;

  vtkSmartPointer< vtkPolyData > signedGlyphs;
  vtkSmartPointer< vtkPolyData > glyphs;
  createGlyphs( argv[1], signedGlyphs, glyphs );

  vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  mapper->SetInput( glyphs );

  vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
  actor->SetMapper( mapper );

  vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();
  outline->SetInput( glyphs );
  vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  outlineMapper->SetInputConnection( outline->GetOutputPort() );
  vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
  outlineActor->SetMapper( outlineMapper );

  vtkSmartPointer< vtkPolyDataMapper > signedMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  signedMapper->SetInput( signedGlyphs );

  vtkSmartPointer< vtkActor > signedActor = vtkSmartPointer< vtkActor >::New();
  signedActor->SetMapper( signedMapper );

  vtkSmartPointer<vtkOutlineFilter> signedOutline = vtkSmartPointer<vtkOutlineFilter>::New();
  signedOutline->SetInput( signedGlyphs );
  vtkSmartPointer<vtkPolyDataMapper> signedOutlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  signedOutlineMapper->SetInputConnection( signedOutline->GetOutputPort() );
  vtkSmartPointer<vtkActor> signedOutlineActor = vtkSmartPointer<vtkActor>::New();
  signedOutlineActor->SetMapper( signedOutlineMapper );

  // the axes
  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();

  // the plates
  vtkSmartPointer<vtkCubeSource> bottomPlate = vtkSmartPointer<vtkCubeSource>::New();
  double xLength = (blockBounds[1] - blockBounds[0]) * 1.2;
  const double yLength = (blockBounds[3] - blockBounds[2]) / 10.0;
  double zLength = (blockBounds[5] - blockBounds[4]) * 1.2;
  bottomPlate->SetXLength( xLength );
  bottomPlate->SetYLength( yLength );
  bottomPlate->SetZLength( zLength );
  bottomPlate->SetCenter( (blockBounds[1] - blockBounds[0]) / 2.0,
    -yLength,
    (blockBounds[5] - blockBounds[4])/ 2.0 );
  vtkSmartPointer<vtkPolyDataMapper> bottomPlateMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  bottomPlateMapper->SetInputConnection( bottomPlate->GetOutputPort() );
  vtkSmartPointer<vtkActor> bottomPlateActor = vtkSmartPointer<vtkActor>::New();
  bottomPlateActor->SetMapper( bottomPlateMapper );

  vtkSmartPointer<vtkCubeSource> topPlate = vtkSmartPointer<vtkCubeSource>::New();
  topPlate->SetXLength( xLength );
  topPlate->SetYLength( yLength );
  topPlate->SetZLength( zLength );
  topPlate->SetCenter( (blockBounds[1] - blockBounds[0]) / 2.0,
    blockBounds[3] + yLength,
    (blockBounds[5] - blockBounds[4])/ 2.0 );
  vtkSmartPointer<vtkPolyDataMapper> topPlateMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  topPlateMapper->SetInputConnection( topPlate->GetOutputPort() );
  vtkSmartPointer<vtkActor> topPlateActor = vtkSmartPointer<vtkActor>::New();
  topPlateActor->SetMapper( topPlateMapper );

  // the arrows
  double arrowColor[3];
  arrowColor[0] = 0.1;
  arrowColor[1] = 0.7;
  arrowColor[2] = 0.9;
  double arrowCenter[3];

  vtkSmartPointer<vtkArrowSource> topArrow = vtkSmartPointer<vtkArrowSource>::New();
  if( tension )
    {
    topArrow->InvertOff();
    }
  else
    {
    topArrow->InvertOn();
    }
  vtkSmartPointer<vtkPolyDataMapper> topArrowMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  topArrowMapper->SetInputConnection( topArrow->GetOutputPort() );
  vtkSmartPointer<vtkActor> topArrowActor = vtkSmartPointer<vtkActor>::New();
  topArrowActor->GetProperty()->SetColor( arrowColor );
  topArrowActor->SetMapper( topArrowMapper );
  vtkSmartPointer<vtkTransform> topArrowTransform = vtkSmartPointer<vtkTransform>::New();
  arrowCenter[0] = (blockBounds[1] - blockBounds[0]) / 2.0;
  arrowCenter[1] = blockBounds[3] + yLength * 3.0/2.0;
  arrowCenter[2] = (blockBounds[5] - blockBounds[4]) / 2.0;
  topArrowTransform->Translate( arrowCenter );
  topArrowTransform->RotateZ( 90.0 );
  topArrowActor->SetUserMatrix( topArrowTransform->GetMatrix() );

  vtkSmartPointer<vtkArrowSource> bottomArrow = vtkSmartPointer<vtkArrowSource>::New();
  if( tension )
    {
    bottomArrow->InvertOn();
    }
  else
    {
    bottomArrow->InvertOff();
    }
  vtkSmartPointer<vtkPolyDataMapper> bottomArrowMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  bottomArrowMapper->SetInputConnection( bottomArrow->GetOutputPort() );
  vtkSmartPointer<vtkActor> bottomArrowActor = vtkSmartPointer<vtkActor>::New();
  bottomArrowActor->GetProperty()->SetColor( arrowColor );
  bottomArrowActor->SetMapper( bottomArrowMapper );
  vtkSmartPointer<vtkTransform> bottomArrowTransform = vtkSmartPointer<vtkTransform>::New();
  arrowCenter[0] = (blockBounds[1] - blockBounds[0]) / 2.0;
  arrowCenter[1] = blockBounds[2] - yLength * 3.0/2.0 - 1.0;
  arrowCenter[2] = (blockBounds[5] - blockBounds[4]) / 2.0;
  bottomArrowTransform->Translate( arrowCenter );
  bottomArrowTransform->RotateZ( 90.0 );
  bottomArrowActor->SetUserMatrix( bottomArrowTransform->GetMatrix() );

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground( 0.2, 0.2, 0.2 );
  renderer->AddActor( actor );
  renderer->AddActor( outlineActor );
  renderer->AddActor( axes );
  renderer->AddActor( bottomPlateActor );
  renderer->AddActor( topPlateActor );
  renderer->AddActor( topArrowActor );
  renderer->AddActor( bottomArrowActor );
  double leftViewport[4] = { 0.0, 0.0, 0.5, 1.0 };
  renderer->SetViewport( leftViewport );

  vtkSmartPointer<vtkRenderer> signedRenderer = vtkSmartPointer<vtkRenderer>::New();
  signedRenderer->SetBackground( 0.2, 0.2, 0.2 );
  signedRenderer->AddActor( signedActor );
  signedRenderer->AddActor( signedOutlineActor );
  signedRenderer->AddActor( axes );
  signedRenderer->AddActor( bottomPlateActor );
  signedRenderer->AddActor( topPlateActor );
  signedRenderer->AddActor( topArrowActor );
  signedRenderer->AddActor( bottomArrowActor );
  double rightViewport[4] = { 0.5, 0.0, 1.0, 1.0 };
  signedRenderer->SetViewport( rightViewport );

  vtkSmartPointer< vtkRenderWindow > renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
  renderWindow->AddRenderer( renderer );
  renderWindow->AddRenderer( signedRenderer );
  renderWindow->SetSize( 1600, 800 );

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  renderWindow->Render();
  // couple the cameras
  vtkCamera * camera = renderer->GetActiveCamera();
  signedRenderer->SetActiveCamera( camera );
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  return 0;
}
