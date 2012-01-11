/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOneSheetedHyperboloidSource.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOneSheetedHyperboloidSource - create a polygonal one-sheeted hyperboloid centered at the origin
// .SECTION Description
// vtkOneSheetedHyperboloidSource creates a one-sheeted hyperboloid (represented by polygons) orientated along
// the z-axis and centered at the origin. The hyperboloid is a representation of the equation
// \f[
//    \frac{x^2}{a^2} + \frac{y^2}{b^2} - \frac{z^2}{c^2} = 1
// \f]
// where \b a, \b b, and \b c parameterize its shape.  These parameters are
// specified with the SetShapeParameters method. The representation is
// truncated at +/- \b ZMax.
// The resolution (polygonal discretization)
// in both the theta and z directions can be specified.
// By default, the surface tessellation of
// the sphere uses triangles; however you can set QuadrilateralTessellation to
// produce a tessellation using quadrilaterals.
// .SECTION Caveats
// A hyperboloid is not a closed surface, and this polygonal representation only
// covers the extent surrounding the origin.  This vtkPolyDataAlgorithm does not
// compute the normals for the polygons.  If they are needed, the
// vtkPolyDataNormals algorithm is a possibility.

#ifndef __vtkOneSheetedHyperboloidSource_h
#define __vtkOneSheetedHyperboloidSource_h

#include "vtkPolyDataAlgorithm.h"

#define VTK_MAX_ONE_SHEETED_HYPERBOLOID_RESOLUTION 1024

class VTK_GRAPHICS_EXPORT vtkOneSheetedHyperboloidSource : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkOneSheetedHyperboloidSource,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct one-sheeted hyperboloid with z_max=0.5 and
  // shape parameters (a, b, and c) = 0.1.
  // Default resolution is 8 in both theta and z directions.  Theta ranges from
  // (0, 360) and z from (-z_max, z_max).
  static vtkOneSheetedHyperboloidSource *New();

  // Description:
  // Set z_max of the hyperboloid.  The hyperboloid extends from -z_max to
  // z_max.  Default is .5.
  vtkSetClampMacro(ZMax,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(ZMax,double);

  // Description:
  // Set the shape parameters for the hyperboloid  If
  // \$\frac{x^2}{a^2} + \frac{y^2}{b^2} - \frac{z^2}{c^2} = 1\$f is the equation of
  // the hyperboloid, then ShapeParameters[0] = a, ShapeParameters[1] = b, nad
  // ShapeParameters[2] = c.  Default is 0.32126,0.32126,0.32126 such that with
  // a z_max of 0.5, the surface area approximately equals the surface area of a
  // sphere with radius 0.5.
  vtkSetVector3Macro(ShapeParameters,double);
  vtkGetVectorMacro(ShapeParameters,double,3);

  // Description:
  // Set the center of the sphere. Default is 0,0,0.
  vtkSetVector3Macro(Center,double);
  vtkGetVectorMacro(Center,double,3);

  // Description:
  // Set the number of points in the rotational direction.
  vtkSetClampMacro(ThetaResolution,int,3,VTK_MAX_ONE_SHEETED_HYPERBOLOID_RESOLUTION);
  vtkGetMacro(ThetaResolution,int);

  // Description:
  // Set the number of points in the Z direction (ranging
  // from -ZMax to ZMax).
  vtkSetClampMacro(ZResolution,int,3,VTK_MAX_ONE_SHEETED_HYPERBOLOID_RESOLUTION);
  vtkGetMacro(ZResolution,int);

  // Description:
  // Cause the sphere to be tessellated with edges along the latitude
  // and longitude lines. If off, triangles are generated at non-polar
  // regions, which results in edges that are not parallel to latitude and
  // longitude lines. If on, quadrilaterals are generated everywhere
  // except at the poles. This can be useful for generating a wireframe
  // sphere with natural latitude and longitude lines.
  vtkSetMacro(QuadrilateralTessellation,int);
  vtkGetMacro(QuadrilateralTessellation,int);
  vtkBooleanMacro(QuadrilateralTessellation,int);

protected:
  vtkOneSheetedHyperboloidSource(int res=8);
  ~vtkOneSheetedHyperboloidSource() {}

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  double ZMax;
  double ShapeParameters[3];
  double Center[3];
  int ThetaResolution;
  int ZResolution;
  int QuadrilateralTessellation;

private:
  vtkOneSheetedHyperboloidSource(const vtkOneSheetedHyperboloidSource&);  // Not implemented.
  void operator=(const vtkOneSheetedHyperboloidSource&);  // Not implemented.
};

#endif
