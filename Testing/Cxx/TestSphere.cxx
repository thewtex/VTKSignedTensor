#include "vtkActor.h"
#include "vtkMassProperties.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkRenderWindow.h"
#include "vtkOutlineFilter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTesting.h"

int TestSphere( int argc, char * argv[] )
{
  vtkSmartPointer< vtkSphereSource > sphere = vtkSmartPointer< vtkSphereSource >::New();
  sphere->SetThetaResolution( 50 );
  sphere->SetPhiResolution( 50 );

  vtkSmartPointer< vtkMassProperties > massProp = vtkSmartPointer< vtkMassProperties >::New();
  massProp->SetInputConnection( sphere->GetOutputPort() );
  massProp->Update();
  std::cout << "The sphere's surface area is: " << massProp->GetSurfaceArea() << std::endl;

  vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  mapper->SetInputConnection( sphere->GetOutputPort() );

  vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
  actor->SetMapper( mapper );
  vtkProperty * property = actor->GetProperty();
  property->SetColor( 0.1, 0.7, 0.9 );

  vtkSmartPointer< vtkOutlineFilter > outline = vtkSmartPointer< vtkOutlineFilter >::New();
  outline->SetInputConnection( sphere->GetOutputPort() );
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
