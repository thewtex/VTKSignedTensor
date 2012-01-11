/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSignedEigenvalueTensorGlyph.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSignedEigenvalueTensorGlyph.h"


#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"

vtkStandardNewMacro(vtkSignedEigenvalueTensorGlyph);

// Construct object with scaling on and scale factor 1.0. Eigenvalues are
// extracted, glyphs are colored with input scalar data, and logarithmic
// scaling is turned off.
vtkSignedEigenvalueTensorGlyph::vtkSignedEigenvalueTensorGlyph()
{
  this->Scaling = 1;
  this->ScaleFactor = 1.0;
  this->ExtractEigenvalues = 1;
  this->ColorGlyphs = 1;
  this->ColorMode = COLOR_BY_SCALARS;
  this->ClampScaling = 0;
  this->Length = 1.0;

  this->SetNumberOfInputPorts(2);

  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);

  // by default, process active point scalars
  this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
vtkSignedEigenvalueTensorGlyph::~vtkSignedEigenvalueTensorGlyph()
{
}

//----------------------------------------------------------------------------
int vtkSignedEigenvalueTensorGlyph::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
      vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *sourceInfo;
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  const int numberOfGeometrySources = inputVector[1]->GetNumberOfInformationObjects();
  for( int ii = 0; ii < numberOfGeometrySources; ++ii )
    {
    sourceInfo = inputVector[1]->GetInformationObject(ii);
    if (sourceInfo)
      {
      sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), 0);
      sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);
      sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 0);
      }
    }

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}

//----------------------------------------------------------------------------
int vtkSignedEigenvalueTensorGlyph::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  const int numberOfGeometrySources = inputVector[1]->GetNumberOfInformationObjects();
  const int numberOfExpectedGeometrySources = 4;
  if( numberOfGeometrySources != numberOfExpectedGeometrySources )
    {
    vtkErrorMacro(<< "The number of geometry sources, " << numberOfGeometrySources <<
      " did not match the number of expected number of geometry sources, " << numberOfExpectedGeometrySources
      << ". Please add more sources with SetSourceConnection().");
    }
  // An array containing the glyphs.  The elements are:
  //         0         | geometry for zero negative eigenvalues
  //         1         | geometry for one negative eigenvalue
  //         2         | geometry for two negative eigenvalue
  //         3         | geometry for three negative eigenvalue
  vtkPolyData * sources[numberOfGeometrySources];
  for( int ii = 0; ii < numberOfGeometrySources; ++ii )
    {
    vtkInformation *sourceInfo = inputVector[1]->GetInformationObject( ii );
    if( !sourceInfo )
      {
      vtkErrorMacro(<< "Geometry source " << ii << " is NULL.");
      }
    sources[ii] = vtkPolyData::SafeDownCast(
      sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
    }
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  double x[3], s;
  vtkIdList *cellPts;
  int npts;
  int dir, eigen_dir, symmetric_dir;
  double xv[3], yv[3], zv[3];
  double maxScale;


  vtkDebugMacro(<<"Generating tensor glyphs");

  vtkPointData *inputPointData = input->GetPointData();
  vtkPointData *outputPointData = output->GetPointData();
  vtkDataArray *inputTensors = this->GetInputArrayToProcess(0, inputVector);
  vtkDataArray *inputScalars = this->GetInputArrayToProcess(1, inputVector);
  const vtkIdType numberOfInputPts = input->GetNumberOfPoints();

  if ( !inputTensors || numberOfInputPts < 1 )
    {
    vtkErrorMacro(<<"No data to glyph!");
    return 1;
    }

  //
  // Allocate storage for output PolyData
  //
  vtkIdType maxNumberOfSourcePoints = 0;
  vtkIdType maxNumberOfSourceCells = 0;
  int sourcesMaxCellSize = 0;
  bool sourcesPointDataHaveNormals = false;
  for( int ii = 0; ii < numberOfGeometrySources; ++ii )
    {
    vtkPoints *sourcePts = sources[ii]->GetPoints();
    vtkIdType numberOfSourcePts = sourcePts->GetNumberOfPoints();
    if( numberOfSourcePts > maxNumberOfSourcePoints )
      {
      maxNumberOfSourcePoints = numberOfSourcePts;
      }

    const vtkIdType numberOfSourceCells = sources[ii]->GetNumberOfCells();
    if( numberOfSourceCells > maxNumberOfSourceCells )
      {
      maxNumberOfSourceCells = numberOfSourceCells;
      }

    const int maxCellSize = sources[ii]->GetMaxCellSize();
    if( maxCellSize > sourcesMaxCellSize )
      {
      sourcesMaxCellSize = maxCellSize;
      }

    vtkCellArray *sourceCells;
    vtkIdType numberOfCells;
    // Setting up for calls to PolyData::InsertNextCell()
    sourceCells = sources[ii]->GetVerts();
    numberOfCells = sourceCells->GetNumberOfCells();
    if ( numberOfCells > 0 )
      {
      vtkCellArray *cells = vtkCellArray::New();
      cells->Allocate( numberOfInputPts * sourceCells->GetSize() );
      output->SetVerts( cells );
      cells->Delete();
      }
    sourceCells = sources[ii]->GetLines();
    numberOfCells = sourceCells->GetNumberOfCells();
    if ( numberOfCells > 0 )
      {
      vtkCellArray *cells = vtkCellArray::New();
      cells->Allocate( numberOfInputPts * sourceCells->GetSize() );
      output->SetLines( cells );
      cells->Delete();
      }
    sourceCells = sources[ii]->GetPolys();
    numberOfCells = sourceCells->GetNumberOfCells();
    if ( numberOfCells > 0 )
      {
      vtkCellArray *cells = vtkCellArray::New();
      cells->Allocate( numberOfInputPts * sourceCells->GetSize() );
      output->SetPolys( cells );
      cells->Delete();
      }
    sourceCells = sources[ii]->GetStrips();
    numberOfCells = sourceCells->GetNumberOfCells();
    if ( numberOfCells > 0 )
      {
      vtkCellArray *cells = vtkCellArray::New();
      cells->Allocate( numberOfInputPts * sourceCells->GetSize() );
      output->SetStrips( cells );
      cells->Delete();
      }

    vtkPointData * pointData = sources[ii]->GetPointData();
    if( pointData->GetNormals() )
      {
      sourcesPointDataHaveNormals = true;
      }
    } // end for each glyph source
  vtkPoints *newPts = vtkPoints::New();
  newPts->Allocate( numberOfInputPts * maxNumberOfSourcePoints );
  vtkIdType *pts = new vtkIdType[sourcesMaxCellSize];

  vtkFloatArray *newNormals = NULL;
  vtkFloatArray *newScalars = NULL;
  // generate scalars if eigenvalues are chosen or if scalars exist.
  if (this->ColorGlyphs &&
      ((this->ColorMode == COLOR_BY_EIGENVALUES) ||
       (inputScalars && (this->ColorMode == COLOR_BY_SCALARS)) ) )
    {
    newScalars = vtkFloatArray::New();
    newScalars->Allocate( numberOfInputPts * maxNumberOfSourcePoints );
    if (this->ColorMode == COLOR_BY_EIGENVALUES)
      {
      newScalars->SetName("MaxEigenvalue");
      }
    else
      {
      newScalars->SetName(inputScalars->GetName());
      }
    }
  if ( sourcesPointDataHaveNormals )
    {
    newNormals = vtkFloatArray::New();
    newNormals->SetNumberOfComponents(3);
    newNormals->SetName("Normals");
    newNormals->Allocate( 3 * numberOfInputPts * maxNumberOfSourcePoints );
    }

  //
  // Traverse all Input points, transforming glyph at Source points
  //
  vtkTransform *transform = vtkTransform::New();
  transform->PreMultiply();
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  double tensor[9];
  // set up working matrices
  double * tensorMatrix[3];
  double tensorMatrixRow0[3];
  double tensorMatrixRow1[3];
  double tensorMatrixRow2[3];
  tensorMatrix[0] = tensorMatrixRow0;
  tensorMatrix[1] = tensorMatrixRow1;
  tensorMatrix[2] = tensorMatrixRow2;

  double eigenvalues[3];

  double * eigenvectors[3];
  double eigenvectorsRow0[3];
  double eigenvectorsRow1[3];
  double eigenvectorsRow2[3];
  eigenvectors[0] = eigenvectorsRow0;
  eigenvectors[1] = eigenvectorsRow1;
  eigenvectors[2] = eigenvectorsRow2;
  double xEigenvector[3];
  double yEigenvector[3];
  double zEigenvector[3];

  for (vtkIdType inputPointId = 0; inputPointId < numberOfInputPts; ++inputPointId)
    {
    // Translation is postponed

    inputTensors->GetTuple(inputPointId, tensor);

    // compute orientation vectors and scale factors from tensor
    if ( this->ExtractEigenvalues ) // extract appropriate eigenfunctions
      {
      for (int jj = 0; jj < 3; ++jj)
        {
        for (int ii = 0; ii < 3; ++ii)
          {
          tensorMatrix[ii][jj] = tensor[ii+3*jj];
          }
        }
      vtkMath::Jacobi(tensorMatrix, eigenvalues, eigenvectors);

      //copy eigenvectors
      xEigenvector[0] = eigenvectors[0][0];
      xEigenvector[1] = eigenvectors[1][0];
      xEigenvector[2] = eigenvectors[2][0];
      yEigenvector[0] = eigenvectors[0][1];
      yEigenvector[1] = eigenvectors[1][1];
      yEigenvector[2] = eigenvectors[2][1];
      zEigenvector[0] = eigenvectors[0][2];
      zEigenvector[1] = eigenvectors[1][2];
      zEigenvector[2] = eigenvectors[2][2];
      }
    else //use tensor columns as eigenvectors
      {
      for (int ii = 0; ii < 3; ++ii)
        {
        xEigenvector[ii] = tensor[ii];
        yEigenvector[ii] = tensor[ii+3];
        zEigenvector[ii] = tensor[ii+6];
        }
      eigenvalues[0] = vtkMath::Normalize(xEigenvector);
      eigenvalues[1] = vtkMath::Normalize(yEigenvector);
      eigenvalues[2] = vtkMath::Normalize(zEigenvector);
      }

    // compute scale factors
    eigenvalues[0] *= this->ScaleFactor;
    eigenvalues[1] *= this->ScaleFactor;
    eigenvalues[2] *= this->ScaleFactor;

    double maxScale = 0.0;
    if ( this->ClampScaling )
      {
      for (int ii = 0; ii < 3; ++ii)
        {
        if ( maxScale < fabs(eigenvalues[ii]) )
          {
          maxScale = fabs(eigenvalues[ii]);
          }
        }
      if ( maxScale > this->MaxScaleFactor )
        {
        maxScale = this->MaxScaleFactor / maxScale;
        for (int ii = 0; ii < 3; ++ii)
          {
          eigenvalues[ii] *= maxScale; //preserve overall shape of glyph
          }
        }
      }

    // normalization is postponed

    // make sure scale is okay (non-zero) and scale data
    maxScale = 0.0;
    for (int ii = 0; ii < 3; ++ii)
      {
      if ( eigenvalues[ii] > maxScale )
        {
        maxScale = eigenvalues[ii];
        }
      }
    if ( maxScale == 0.0 )
      {
      maxScale = 1.0;
      }
    for (int ii = 0; ii < 3; ++ii)
      {
      if ( eigenvalues[ii] == 0.0 )
        {
        eigenvalues[ii] = maxScale * 1.0e-06;
        }
      }

    int numberOfNegativeEigenvalues = 0;
    for (int ii = 0; ii < 3; ++ii)
      {
      if ( eigenvalues[ii] < 0.0 )
        {
        ++numberOfNegativeEigenvalues;
        }
      }

    vtkPolyData * source = sources[numberOfNegativeEigenvalues];
    vtkIdType numberOfSourcePoints = source->GetNumberOfPoints();
    vtkIdType ptIncr = inputPointId * numberOfSourcePoints;
    //
    // copy all topology (transformation independent)
    //
    for (vtkIdType cellId = 0; cellId < numSourceCells; cellId++)
      {
      cell = this->GetSource()->GetCell(cellId);
      cellPts = cell->GetPointIds();
      npts = cellPts->GetNumberOfIds();
      for (dir=0; dir < numDirs; dir++)
        {
        // This variable may be removed, but that
        // will not improve readability
        vtkIdType subIncr = ptIncr + dir*numSourcePts;
        for (i=0; i < npts; i++)
          {
          pts[i] = cellPts->GetId(i) + subIncr;
          }
        output->InsertNextCell(cell->GetCellType(),npts,pts);
        }
      }

    // Now do the real work for each "direction"

    for (dir=0; dir < numDirs; dir++)
      {
      eigen_dir = dir%(this->ThreeGlyphs?3:1);
      symmetric_dir = dir/(this->ThreeGlyphs?3:1);

      // Remove previous scales ...
      trans->Identity();

      // translate Source to Input point
      input->GetPoint(inPtId, x);
      trans->Translate(x[0], x[1], x[2]);

      // normalized eigenvectors rotate object for eigen direction 0
      matrix->Element[0][0] = xv[0];
      matrix->Element[0][1] = yv[0];
      matrix->Element[0][2] = zv[0];
      matrix->Element[1][0] = xv[1];
      matrix->Element[1][1] = yv[1];
      matrix->Element[1][2] = zv[1];
      matrix->Element[2][0] = xv[2];
      matrix->Element[2][1] = yv[2];
      matrix->Element[2][2] = zv[2];
      trans->Concatenate(matrix);

      if (eigen_dir == 1)
        {
        trans->RotateZ(90.0);
        }

      if (eigen_dir == 2)
        {
        trans->RotateY(-90.0);
        }

      if (this->ThreeGlyphs)
        {
        trans->Scale(w[eigen_dir], this->ScaleFactor, this->ScaleFactor);
        }
      else
        {
        trans->Scale(w[0], w[1], w[2]);
        }

      // Mirror second set to the symmetric position
      if (symmetric_dir == 1)
        {
        trans->Scale(-1.,1.,1.);
        }

      // if the eigenvalue is negative, shift to reverse direction.
      // The && is there to ensure that we do not change the
      // old behaviour of vtkSignedEigenvalueTensorGlyphs (which only used one dir),
      // in case there is an oriented glyph, e.g. an arrow.
      if (w[eigen_dir] < 0 && numDirs > 1)
        {
        trans->Translate(-this->Length, 0., 0.);
        }

      // multiply points (and normals if available) by resulting
      // matrix
      trans->TransformPoints(sourcePts,newPts);

      // Apply the transformation to a series of points,
      // and append the results to outPts.
      if ( newNormals )
        {
        // a negative determinant means the transform turns the
        // glyph surface inside out, and its surface normals all
        // point inward. The following scale corrects the surface
        // normals to point outward.
        if (trans->GetMatrix()->Determinant() < 0)
          {
          trans->Scale(-1.0,-1.0,-1.0);
          }
        trans->TransformNormals(sourceNormals,newNormals);
        }

        // Copy point data from source
      if ( this->ColorGlyphs && inScalars &&
           (this->ColorMode == COLOR_BY_SCALARS) )
        {
        s = inScalars->GetComponent(inPtId, 0);
        for (i=0; i < numSourcePts; i++)
          {
          newScalars->InsertTuple(ptIncr+i, &s);
          }
        }
      else if (this->ColorGlyphs &&
               (this->ColorMode == COLOR_BY_EIGENVALUES) )
        {
        // If ThreeGlyphs is false we use the first (largest)
        // eigenvalue as scalar.
        s = w[eigen_dir];
        for (i=0; i < numSourcePts; i++)
          {
          newScalars->InsertTuple(ptIncr+i, &s);
          }
        }
      else
        {
        for (i=0; i < numSourcePts; i++)
          {
          outPD->CopyData(pd,i,ptIncr+i);
          }
        }
      ptIncr += numSourcePts;
      }
    }
  vtkDebugMacro(<<"Generated " << numPts <<" tensor glyphs");
  //
  // Update output and release memory
  //
  delete [] pts;

  output->SetPoints(newPts);
  newPts->Delete();

  if ( newScalars )
    {
    int idx = outPD->AddArray(newScalars);
    outPD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    newScalars->Delete();
    }

  if ( newNormals )
    {
    outPD->SetNormals(newNormals);
    newNormals->Delete();
    }

  output->Squeeze();
  trans->Delete();
  matrix->Delete();

  return 1;
}

//----------------------------------------------------------------------------
void vtkSignedEigenvalueTensorGlyph::SetSourceConnection(int id, vtkAlgorithmOutput* algOutput)
{
  if (id < 0)
    {
    vtkErrorMacro("Bad index " << id << " for source.");
    return;
    }

  int numConnections = this->GetNumberOfInputConnections(1);
  if (id < numConnections)
    {
    this->SetNthInputConnection(1, id, algOutput);
    }
  else if (id == numConnections && algOutput)
    {
    this->AddInputConnection(1, algOutput);
    }
  else if (algOutput)
    {
    vtkWarningMacro("The source id provided is larger than the maximum "
                    "source id, using " << numConnections << " instead.");
    this->AddInputConnection(1, algOutput);
    }
}

//----------------------------------------------------------------------------
int vtkSignedEigenvalueTensorGlyph::FillInputPortInformation(int port, vtkInformation *info)
{
  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
    }
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
void vtkSignedEigenvalueTensorGlyph::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Source: " << this->GetSource() << "\n";
  os << indent << "Scaling: " << (this->Scaling ? "On\n" : "Off\n");
  os << indent << "Scale Factor: " << this->ScaleFactor << "\n";
  os << indent << "Extract Eigenvalues: " << (this->ExtractEigenvalues ? "On\n" : "Off\n");
  os << indent << "Color Glyphs: " << (this->ColorGlyphs ? "On\n" : "Off\n");
  os << indent << "Color Mode: " << this->ColorMode << endl;
  os << indent << "Clamp Scaling: " << (this->ClampScaling ? "On\n" : "Off\n");
  os << indent << "Max Scale Factor: " << this->MaxScaleFactor << "\n";
  os << indent << "Length: " << this->Length << "\n";
}
