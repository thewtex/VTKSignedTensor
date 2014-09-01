#include "vtkActor.h"
#include "vtkMassProperties.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"
#include "vtkOneSheetedHyperboloidSource.h"
#include "vtkOutlineFilter.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include <iostream>

int main( int argc, char * argv[] )
{
  vtkSmartPointer< vtkOneSheetedHyperboloidSource > hyperboloid = vtkSmartPointer< vtkOneSheetedHyperboloidSource >::New();
  hyperboloid->SetThetaResolution( 90 );
  hyperboloid->SetZResolution( 50 );

  vtkSmartPointer< vtkMassProperties > massProp = vtkSmartPointer< vtkMassProperties >::New();
  massProp->SetInputConnection( hyperboloid->GetOutputPort() );
  massProp->Update();
  std::cout << "The hyperboloid's surface area is: " << massProp->GetSurfaceArea() << std::endl;

  vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  mapper->SetInputConnection( hyperboloid->GetOutputPort() );

  vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
  actor->SetMapper( mapper );
  vtkProperty * property = actor->GetProperty();
  property->SetColor( 0.1, 0.7, 0.9 );

  vtkSmartPointer< vtkOutlineFilter > outline = vtkSmartPointer< vtkOutlineFilter >::New();
  outline->SetInputConnection( hyperboloid->GetOutputPort() );
  vtkSmartPointer< vtkPolyDataMapper > outlineMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  outlineMapper->SetInputConnection( outline->GetOutputPort() );
  vtkSmartPointer< vtkActor > outlineActor = vtkSmartPointer< vtkActor >::New();
  outlineActor->SetMapper( outlineMapper );

  vtkSmartPointer< vtkRenderer > renderer = vtkSmartPointer< vtkRenderer >::New();
  renderer->AddActor( actor );
  renderer->AddActor( outlineActor );
  renderer->SetBackground( 0.1, 0.1, 0.1 );

  vtkSmartPointer< vtkRenderWindow > renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
  renderWindow->AddRenderer( renderer );
  renderWindow->SetSize( 600, 600 );

  vtkSmartPointer< vtkRenderWindowInteractor > renderWindowInteractor = vtkSmartPointer< vtkRenderWindowInteractor >::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );

  renderWindowInteractor->Start();

  return 0;
}
