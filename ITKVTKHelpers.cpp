/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "ITKVTKHelpers.h"

// Submodules
#include "ITKHelpers/ITKHelpers.h"

// VTK
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>

// ITK
#include "itkVectorMagnitudeImageFilter.h"

namespace ITKVTKHelpers
{

void CreateTransparentVTKImage(const itk::Size<2>& size, vtkImageData* const outputImage)
{
  outputImage->SetDimensions(size[0], size[1], 1);
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 4);

  for(unsigned int i = 0; i < size[0]; ++i)
    {
    for(unsigned int j = 0; j < size[1]; ++j)
      {
      unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(i, j ,0));
      pixel[0] = 0;
      pixel[1] = 0;
      pixel[2] = 0;
      pixel[3] = 0; // transparent
      }
    }
  outputImage->Modified();
}

void SetRegionCenterPixel(vtkImageData* const image, const itk::ImageRegion<2>& region, const unsigned char color[3])
{
  int dims[3];
  image->GetDimensions(dims);

  unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(region.GetIndex()[0] + region.GetSize()[0]/2,
                                                                             region.GetIndex()[1] + region.GetSize()[1]/2, 0));
  pixel[0] = color[0];
  pixel[1] = color[1];
  pixel[2] = color[2];
  pixel[3] = 255; // visible
}

void ITKImageToVTKVectorFieldImage(const FloatVector2ImageType* const image, vtkImageData* const outputImage)
{
  //std::cout << "ITKImagetoVTKVectorFieldImage()" << std::endl;

  // Setup and allocate the image data
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);
  outputImage->AllocateScalars(VTK_FLOAT, 3);// We really want this to be 2, but VTK complains, so we must add a 3rd component (0) to every pixel

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<FloatVector2ImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    float* pixel = static_cast<float*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                     imageIterator.GetIndex()[1],0));

    FloatVector2ImageType::PixelType inputPixel = imageIterator.Get();
    pixel[0] = inputPixel[0];
    pixel[1] = inputPixel[1];
    pixel[2] = 0;

    ++imageIterator;
    }

  outputImage->GetPointData()->SetActiveVectors("ImageScalars");
  outputImage->Modified();
}


void ConvertNonZeroPixelsToVectors(const FloatVector2ImageType* const vectorImage, vtkPolyData* const output)
{
  vtkSmartPointer<vtkFloatArray> vectors = vtkSmartPointer<vtkFloatArray>::New();
  vectors->SetNumberOfComponents(3);
  vectors->SetName("Vectors");

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<FloatVector2ImageType> imageIterator(vectorImage, vectorImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    FloatVector2ImageType::PixelType inputPixel = imageIterator.Get();
    if(inputPixel.GetNorm() > .05)
      {
      float v[3];
      v[0] = inputPixel[0];
      v[1] = inputPixel[1];
      v[2] = 0;
      vectors->InsertNextTupleValue(v);

      points->InsertNextPoint(imageIterator.GetIndex()[0], imageIterator.GetIndex()[1], 0);
      }

    ++imageIterator;
    }

  output->SetPoints(points);
  output->GetPointData()->SetVectors(vectors);
  output->Modified();
}


void BlankAndOutlineRegion(vtkImageData* const image, const itk::ImageRegion<2>& region, const unsigned char value[3])
{
  BlankRegion(image, region);
  OutlineRegion(image, region, value);
}


void BlankRegion(vtkImageData* image, const itk::ImageRegion<2>& region)
{
  // Blank the image
  for(unsigned int i = region.GetIndex()[0]; i < region.GetIndex()[0] + region.GetSize()[0]; ++i)
    {
    for(unsigned int j = region.GetIndex()[1]; j < region.GetIndex()[1] + region.GetSize()[1]; ++j)
      {
      unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(i, region.GetIndex()[1], 0));
      pixel[0] = 0;
      pixel[1] = 0;
      pixel[2] = 0;
      pixel[3] = 0; // transparent
      }
    }
  image->Modified();
}

void OutlineRegion(vtkImageData* const image, const itk::ImageRegion<2>& region, const unsigned char color[3])
{
//   std::cout << "Outlining region: " << region << std::endl;
//    std::cout << "Outline color: " << static_cast<int>(color[0])
//              << " " << static_cast<int>(color[1]) << " "
//              << static_cast<int>(color[2]) << std::endl;
//   std::cout << "Image components: " << image->GetNumberOfScalarComponents() << std::endl;

  unsigned int leftEdge = region.GetIndex()[0];
  //std::cout << "leftEdge : " << leftEdge << std::endl;

  unsigned int rightEdge = region.GetIndex()[0] + region.GetSize()[0] - 1;
  //std::cout << "rightEdge : " << rightEdge << std::endl;

  unsigned int topEdge = region.GetIndex()[1];
  //std::cout << "topEdge : " << topEdge << std::endl;

  unsigned int bottomEdge = region.GetIndex()[1] + region.GetSize()[1] - 1;
  //std::cout << "bottomEdge : " << bottomEdge << std::endl;

  // Move along the top and bottom of the region, setting the border pixels.
  unsigned int counter = 0;
  for(unsigned int xpos = leftEdge; xpos <= rightEdge; ++xpos)
    {
    //std::cout << "xpos : " << xpos << std::endl;
    unsigned char* topPixel = static_cast<unsigned char*>(image->GetScalarPointer(xpos, topEdge, 0));
    topPixel[0] = color[0];
    topPixel[1] = color[1];
    topPixel[2] = color[2];
    topPixel[3] = 255; // visible

    unsigned char* bottomPixel = static_cast<unsigned char*>(image->GetScalarPointer(xpos, bottomEdge, 0));
    bottomPixel[0] = color[0];
    bottomPixel[1] = color[1];
    bottomPixel[2] = color[2];
    bottomPixel[3] = 255; // visible

    counter++;
    }

  //std::cout << "Set " << counter << " pixels." << std::endl;
  // Move along the left and right of the region, setting the border pixels.
  for(unsigned int j = region.GetIndex()[1]; j < region.GetIndex()[1] + region.GetSize()[1]; ++j)
    {
    unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(region.GetIndex()[0], j, 0));
    pixel[0] = color[0];
    pixel[1] = color[1];
    pixel[2] = color[2];
    pixel[3] = 255; // visible

    pixel = static_cast<unsigned char*>(image->GetScalarPointer(region.GetIndex()[0] + region.GetSize()[0] - 1, j, 0));
    pixel[0] = color[0];
    pixel[1] = color[1];
    pixel[2] = color[2];
    pixel[3] = 255; // visible
    }

  image->Modified();
}

void CreateVectorImageFromStructuredGridArray(vtkStructuredGrid* const structuredGrid,
                                              const std::string& arrayName, FloatVectorImageType* const outputImage)
{
  itk::ImageRegionIteratorWithIndex<FloatVectorImageType> imageIterator(outputImage, outputImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  vtkDataArray* dataArray = structuredGrid->GetPointData()->GetArray(arrayName.c_str());
  if(!dataArray)
    {
    std::stringstream ss;
    ss << "Array \"" << arrayName << "\" does not exist!";
    throw std::runtime_error(ss.str());
    }

  int dimensions[3];
  structuredGrid->GetDimensions(dimensions);

  outputImage->SetNumberOfComponentsPerPixel(dataArray->GetNumberOfComponents());
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> imageSize = {{static_cast<itk::SizeValueType>(dimensions[0]),
                             static_cast<itk::SizeValueType>(dimensions[1])}};
  itk::ImageRegion<2> region(corner, imageSize);
  outputImage->SetRegions(region);
  outputImage->Allocate();

  while(!imageIterator.IsAtEnd())
  {
    int queryPoint[3] = {static_cast<int>(imageIterator.GetIndex()[0]),
                         static_cast<int>(imageIterator.GetIndex()[1]), 0};
    vtkIdType pointId = vtkStructuredData::ComputePointId(dimensions, queryPoint);

    FloatVectorImageType::PixelType p;
    p.SetSize(dataArray->GetNumberOfComponents());

    if(structuredGrid->IsPointVisible(pointId))
    {
      std::vector<double> value(dataArray->GetNumberOfComponents());
      dataArray->GetTuple(pointId, value.data());

      for(vtkIdType component = 0; component < dataArray->GetNumberOfComponents(); ++component)
      {
        p[component] = value[component];
      }

    }
    else
    {
      for(vtkIdType component = 0; component < dataArray->GetNumberOfComponents(); ++component)
      {
        p[component] = 0;
      }
    }

    imageIterator.Set(p);
    ++imageIterator;
  }
}


void CreateScalarImageFromStructuredGridArray(vtkStructuredGrid* const structuredGrid, const std::string& arrayName, FloatScalarImageType* const outputImage)
{
  itk::ImageRegionIteratorWithIndex<FloatScalarImageType> imageIterator(outputImage, outputImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  vtkDataArray* dataArray = structuredGrid->GetPointData()->GetArray(arrayName.c_str());
  if(!dataArray)
    {
    std::stringstream ss;
    ss << "Array \"" << arrayName << "\" does not exist!";
    throw std::runtime_error(ss.str());
    }

  int dimensions[3];
  structuredGrid->GetDimensions(dimensions);

  outputImage->SetNumberOfComponentsPerPixel(1);
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> imageSize = {{static_cast<itk::SizeValueType>(dimensions[0]),
                             static_cast<itk::SizeValueType>(dimensions[1])}};
  itk::ImageRegion<2> region(corner, imageSize);
  outputImage->SetRegions(region);
  outputImage->Allocate();

  while(!imageIterator.IsAtEnd())
    {
    int queryPoint[3] = {static_cast<int>(imageIterator.GetIndex()[0]),
                         static_cast<int>(imageIterator.GetIndex()[1]), 0};
    vtkIdType pointId = vtkStructuredData::ComputePointId(dimensions, queryPoint);

    FloatScalarImageType::PixelType p;

    if(structuredGrid->IsPointVisible(pointId))
      {
      double value[1];
      dataArray->GetTuple(pointId, value);

      p = value[0];

      imageIterator.Set(p);
      }
    else
      {
      p = 0;
      }

    imageIterator.Set(p);
    ++imageIterator;
    }
}


void SetPixels(vtkImageData* const VTKImage, const std::vector<itk::Index<2> >& pixels, const unsigned char color[3])
{
  int* dims = VTKImage->GetDimensions();

  for(unsigned int i = 0; i < pixels.size(); ++i)
    {
    if(pixels[i][0] >= dims[0] || pixels[i][1] >= dims[1]) // The pixel is out of bounds
      {
      continue;
      }
    unsigned char* pixel = static_cast<unsigned char*>(VTKImage->GetScalarPointer(pixels[i][0],pixels[i][1],0));
    pixel[0] = color[0];
    pixel[1] = color[1];
    pixel[2] = color[2];
    // Make sure the pixel is not transparent
    if(VTKImage->GetNumberOfScalarComponents() == 4)
      {
      pixel[3] = 255;
      }
    }

}

std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* const polydata)
{
  return PointsToPixelList(polydata->GetPoints());
}

std::vector<itk::Index<2> > PointsToPixelList(vtkPoints* const points)
{
  // The points of the polydata are floating point values, we must convert them to pixel indices.

  //std::cout << "Enter PolyDataToPixelList()" << std::endl;
  std::cout << "There are " << points->GetNumberOfPoints() << " points." << std::endl;

  // Convert vtkPoints to indices
  //std::cout << "Converting vtkPoints to indices..." << std::endl;
  std::vector<itk::Index<2> > linePoints;
  for(vtkIdType pointId = 0; pointId < points->GetNumberOfPoints(); ++pointId)
    {
    itk::Index<2> index;
    double p[3];
    points->GetPoint(pointId, p);
    // std::cout << "point " << pointId << " : " << p[0] << " " << p[1] << " " << p[2] << std::endl;

    // Use itk::Math::Round instead of round() for cross-platform compatibility (specifically,
    // VS2010 does not have round() in cmath)
    index[0] = static_cast<itk::Index<2>::IndexValueType>(itk::Math::Round<double,double>(p[0]));
    index[1] = static_cast<itk::Index<2>::IndexValueType>(itk::Math::Round<double,double>(p[1]));
    if(linePoints.size() == 0)
      {
      linePoints.push_back(index);
      continue;
      }

    // Don't duplicate indices of points acquired in a row that round to the same pixel.
    if(index != linePoints[linePoints.size() - 1])
      {
      linePoints.push_back(index);
      }
    }

  if(linePoints.size() < 2)
    {
    std::cerr << "Cannot draw a lines between " << linePoints.size() << " points." << std::endl;
    return linePoints;
    }

  // Compute the indices between every pair of points
  //std::cout << "Computing the indices between every pair of points..." << std::endl;
  std::vector<itk::Index<2> > allIndices;
  for(unsigned int linePointId = 1; linePointId < linePoints.size(); linePointId++)
    {
    //std::cout << "Getting the indices..." << std::endl;
    itk::Index<2> index0 = linePoints[linePointId-1];
    itk::Index<2> index1 = linePoints[linePointId];

    if(index0 == index1)
      {
      std::cout << "Can't draw a line between the same pixels (" << index0 << " and " << index1 << "!" << std::endl;
      continue;
      }

    //std::cout << "Constructing the line..." << std::endl;
    itk::BresenhamLine<2> line;
    std::vector<itk::Index<2> > indices = line.BuildLine(index0, index1);
    //std::cout << "Saving indices..." << std::endl;
    for(unsigned int i = 0; i < indices.size(); i++)
      {
      allIndices.push_back(indices[i]);
      }

    } // end for loop over line segments

  //std::cout << "Exit PolyDataToPixelList()" << std::endl;
  return allIndices;
}

void PixelListToPoints(const std::vector<itk::Index<2> >& pixels, vtkPoints* const points)
{
  typedef std::vector<itk::Index<2> > PixelContainer;

  for(PixelContainer::const_iterator iter = pixels.begin(); iter != pixels.end(); ++iter)
  {
    points->InsertNextPoint((*iter)[0], (*iter)[1], 0);
  }
}


void ITKRGBImageToVTKImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>* const image,
                           vtkImageData* const outputImage)
{
  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  // Setup and allocate the VTK image
  //outputImage->SetNumberOfScalarComponents(3);
  //outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  //outputImage->AllocateScalars();
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<RGBImageType> imageIterator(image, image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get().GetRed();
    pixel[1] = imageIterator.Get().GetGreen();
    pixel[2] = imageIterator.Get().GetBlue();

    ++imageIterator;
    }

  outputImage->Modified();
}

void InitializeVTKImage(const itk::ImageRegion<2>& region, const unsigned int channels, vtkImageData* outputImage)
{
  // Setup and allocate the VTK image
  outputImage->SetDimensions(region.GetSize()[0],
                             region.GetSize()[1],
                             1);
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, channels);
}

void SetPixelTransparency(vtkImageData* const image, const std::vector<itk::Index<2> >& pixels, const unsigned char value)
{
  if(image->GetNumberOfScalarComponents() != 4)
  {
    std::stringstream ss;
    ss << "Cannot set pixel transparency of an image with " << image->GetNumberOfScalarComponents() << " components.";
    throw std::runtime_error(ss.str());
  }

  for(unsigned int i = 0; i < pixels.size(); ++i)
  {
    unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(pixels[i][0],
                                                                               pixels[i][1],0));
    pixel[3] = value;
  }
}

} // end namespace
