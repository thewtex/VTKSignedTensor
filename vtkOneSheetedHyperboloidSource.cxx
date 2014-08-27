/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOneSheetedHyperboloidSource.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOneSheetedHyperboloidSource.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <math.h>

vtkStandardNewMacro(vtkOneSheetedHyperboloidSource);

vtkOneSheetedHyperboloidSource::vtkOneSheetedHyperboloidSource(int res)
{
  res = res < 4 ? 4 : res;
  this->ShapeParameters[0] = 0.32126;
  this->ShapeParameters[1] = 0.32126;
  this->ShapeParameters[2] = 0.32126;
  this->ZMax = 0.5;
  this->Center[0] = 0.0;
  this->Center[1] = 0.0;
  this->Center[2] = 0.0;

  this->ThetaResolution = res;
  this->ZResolution = res;
  this->QuadrilateralTessellation = 0;

  this->SetNumberOfInputPorts(0);
}

int vtkOneSheetedHyperboloidSource::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the output
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // streaming pieces
  const int piece =
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numPieces =
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  // We will stream over the theta angle.
  if (numPieces > this->ThetaResolution)
    {
    numPieces = this->ThetaResolution;
    }
  if (piece >= numPieces)
    {
    // Although the super class should take care of this,
    // it cannot hurt to check here.
    return 1;
    }

  // I want to modify the ivars resoultion start theta and end theta,
  // so I will make local copies of them.  These might be able to be merged
  // with the other copies of them, ...
  int localThetaResolution = this->ThetaResolution;
  double localStartTheta = 0.0;
  double localEndTheta = 360.0;
  double deltaTheta = (localEndTheta - localStartTheta) / localThetaResolution;

  // Change the ivars based on pieces.
  const int start = piece * localThetaResolution / numPieces;
  const int end = (piece+1) * localThetaResolution / numPieces;
  localEndTheta = localStartTheta + (double)(end) * deltaTheta;
  localStartTheta = localStartTheta + (double)(start) * deltaTheta;
  localThetaResolution = end - start;

  // Set things up; allocate memory
  //
  vtkDebugMacro("OneSheetedHyperboloidSource Executing piece index " << piece
                << " of " << numPieces << " pieces.");

  const int numPts = this->ZResolution * localThetaResolution;
  // creating triangles
  const int numPolys = this->ZResolution * 2 * localThetaResolution;

  vtkPoints *newPoints = vtkPoints::New();
  newPoints->Allocate(numPts);
  vtkCellArray *newPolys = vtkCellArray::New();
  newPolys->Allocate(newPolys->EstimateSize(numPolys, 3));

  // Create hyperboloid
  //

  // Check data, determine increments, and convert to radians

  if (fabs(localStartTheta - localEndTheta) < 360.0)
    {
    ++localThetaResolution;
    }
  localStartTheta *= vtkMath::Pi() / 180.0;
  localEndTheta *= vtkMath::Pi() / 180.0;
  deltaTheta = (localEndTheta - localStartTheta) / localThetaResolution;

  // We calculate the values with the parametric representation:
  // x(u, v) = a cosh(v) cos(u)
  // y(u, v) = b cosh(v) sin(u)
  // z(u, v) = c sinh(v)
  // the parameter u here is theta
  const double vStart = asinh(-this->ZMax / this->ShapeParameters[2]);
  const double vEnd = asinh(this->ZMax / this->ShapeParameters[2]);
  const double deltaV = (vEnd - vStart) / (this->ZResolution - 1);

  this->UpdateProgress(0.1);

  // Create intermediate points
  double theta = localStartTheta;
  double location[3];
  for (int ii = 0; ii < localThetaResolution; ++ii)
    {
    double vParam = vStart;
    for (int jj = 0; jj < this->ZResolution; ++jj)
      {
      location[0] = this->ShapeParameters[0] * cosh(vParam) * cos(theta) + this->Center[0];
      location[1] = this->ShapeParameters[1] * cosh(vParam) * sin(theta) + this->Center[1];
      location[2] = this->ShapeParameters[2] * sinh(vParam) + this->Center[2];
      newPoints->InsertNextPoint(location);
      vParam += deltaV;
      }
    theta += deltaTheta;
    this->UpdateProgress (0.10 + 0.50*ii / static_cast< double >(localThetaResolution));
    }

  // Generate mesh connectivity
  const int base = this->ZResolution * localThetaResolution;

  if (fabs(localStartTheta - localEndTheta) < 2.0 * vtkMath::Pi())
    {
    --localThetaResolution;
    }

  vtkIdType pts[4];
  // bands in-between poles
  for (int ii = 0; ii < localThetaResolution; ++ii)
    {
    for (int jj = 0; jj < this->ZResolution - 1; ++jj)
      {
      pts[0] = this->ZResolution*ii + jj;
      pts[1] = pts[0] + 1;
      pts[2] = ((this->ZResolution*(ii + 1) + jj) % base) + 1;
      if ( !this->QuadrilateralTessellation )
        {
        newPolys->InsertNextCell(3, pts);
        pts[1] = pts[2];
        pts[2] = pts[1] - 1;
        newPolys->InsertNextCell(3, pts);
        }
      else
        {
        pts[3] = pts[2] - 1;
        newPolys->InsertNextCell(4, pts);
        }
      }
    this->UpdateProgress (0.70 + 0.30*ii/static_cast<double>(localThetaResolution));
    }

  // Update ourselves and release memeory
  //
  newPoints->Squeeze();
  output->SetPoints(newPoints);
  newPoints->Delete();

  newPolys->Squeeze();
  output->SetPolys(newPolys);
  newPolys->Delete();

  return 1;
}

//----------------------------------------------------------------------------
void vtkOneSheetedHyperboloidSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Theta Resolution: " << this->ThetaResolution << "\n";
  os << indent << "Z Resolution: " << this->ZResolution << "\n";
  os << indent << "ShapeParameters: (" << this->ShapeParameters[0] << ", "
     << this->ShapeParameters[1] << ", " << this->ShapeParameters[2] << ")\n";
  os << indent << "ZMax: " << this->ZMax << "\n";
  os << indent << "Center: (" << this->Center[0] << ", "
     << this->Center[1] << ", " << this->Center[2] << ")\n";
  os << indent
     << "Quadrilateral Tessellation: " << this->QuadrilateralTessellation << "\n";
}

//----------------------------------------------------------------------------
int vtkOneSheetedHyperboloidSource::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
               -1);

  const double zTerm = 1.0 + this->ZMax * this->ZMax / (this->ShapeParameters[2] * this->ShapeParameters[2]);
  const double xMax = sqrt(this->ShapeParameters[0] * this->ShapeParameters[0] * zTerm);
  const double yMax = sqrt(this->ShapeParameters[1] * this->ShapeParameters[1] * zTerm);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_BOUNDING_BOX(),
               this->Center[0] - xMax,
               this->Center[0] + xMax,
               this->Center[1] - yMax,
               this->Center[1] + yMax,
               this->Center[2] - this->ZMax,
               this->Center[2] + this->ZMax);

  return 1;
}
