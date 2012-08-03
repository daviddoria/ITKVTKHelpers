/*=========================================================================
 *
 *  Copyright David Doria 2012 daviddoria@gmail.com
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

#ifndef ITKVTKHelpers_H
#define ITKVTKHelpers_H

// VTK
class vtkImageData;
class vtkPoints;
class vtkPolyData;
class vtkPoints;
class vtkStructuredGrid;

// ITK
#include "itkCovariantVector.h"
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkRGBPixel.h"
#include "itkSize.h"
#include "itkVectorImage.h"

// Submodules
#include "ITKHelpers/ITKHelpersTypes.h"

namespace ITKVTKHelpers
{
  using namespace ITKHelpersTypes;

/** Convert the points in a polydata to a list of indices. */
std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* const polydata);

/** Convert the points in a vtkPoints object to a list of indices. */
std::vector<itk::Index<2> > PointsToPixelList(vtkPoints* const points);

/** Create a transparent VTK image with the specified size. */
void CreateTransparentVTKImage(const itk::Size<2>& size, vtkImageData* const outputImage);

/** Create an image from the values in an array of the corresponding structured grid points. */
void CreateVectorImageFromStructuredGridArray(vtkStructuredGrid* const structuredGrid, const std::string& arrayName,
                                              FloatVectorImageType* const outputImage);

/** Convert the list of pixels to a vtkPoints object. */
void PixelListToPoints(const std::vector<itk::Index<2> >& pixels, vtkPoints* const points);

/** Create an image from the values in an array of the corresponding structured grid points. */
void CreateScalarImageFromStructuredGridArray(vtkStructuredGrid* const structuredGrid, const std::string& arrayName,
                                              FloatScalarImageType* const outputImage);

/** Set the pixels in 'pixels' to 'color'. */
void SetPixels(vtkImageData* const VTKImage, const std::vector<itk::Index<2> >& pixels,
               const unsigned char color[3]);

/** Set the center pixel of a 'region' in an 'image' to the specified 'color'.
    The region is assumed to have odd dimensions.
*/
void SetRegionCenterPixel(vtkImageData* const image, const itk::ImageRegion<2>& region,
                          const unsigned char color[3]);

/** Initialize a VTK image from an ITK image's size. */
void InitializeVTKImage(const itk::ImageRegion<2>& region, const unsigned int channels,
                        vtkImageData* outputImage);

/** Convert an ITK image with vector pixels to a VTK image. */
template <typename TPixel>
void ITKImageToVTKRGBImage(const itk::Image<itk::CovariantVector<TPixel, 3>, 2>* const image,
                           vtkImageData* const outputImage, const bool alreadyInitialized = false);


/** This function simply drives ITKImagetoVTKRGBImage or ITKImagetoVTKMagnitudeImage based on
  * the number of components of the input. */
template <typename TPixel>
void ITKVectorImageToVTKImageFromDimension(const itk::VectorImage<TPixel, 2>* const image,
                                           vtkImageData* const outputImage);

/** These functions create a VTK image from a multidimensional ITK image. */
template <typename TPixel>
void ITKImageToVTKRGBImage(const itk::VectorImage<TPixel, 2>* const image, vtkImageData* const outputImage,
                           const bool alreadyInitialized = false);

/** Convert an ITK image to a VTK image where each channel is the magnitude of the ITK image. */
template <typename TPixel>
void ITKImageToVTKMagnitudeImage(const itk::VectorImage<TPixel, 2>* const image, vtkImageData* const outputImage);

/** Convert a specified 'channel' of an ITK image to a VTK image. */
template <typename TPixel>
void ITKImageChannelToVTKImage(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                               vtkImageData* const outputImage);

/** Create a VTK image filled with values representing vectors. (There is no concept of a "vector image" in VTK). */
void ITKImageToVTKVectorFieldImage(const FloatVector2ImageType* image, vtkImageData* outputImage);

/** It is often too intensive to glyph every vector in a vector image. In many cases, the vector field may have
  * very large regions of zero vectors. This function creates the vectors for only the non-zero pixels in
  * the vector image. */
void ConvertNonZeroPixelsToVectors(const FloatVector2ImageType* const vectorImage, vtkPolyData* const output);


/** Simply calls BlankRegion followed by OutlineRegion */
void BlankAndOutlineRegion(vtkImageData* const image, const itk::ImageRegion<2>& region,
                           const unsigned char value[3]);

/** Set pixels on the boundary of 'region' in 'image' to 'value'. */
void OutlineRegion(vtkImageData* const image, const itk::ImageRegion<2>& region, const unsigned char value[3]);

/** Convert an ITK image to a VTK image */
void ITKRGBImageToVTKImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>* const image,
                            vtkImageData* const outputImage);

/** Set all pixels in 'region' in 'image' to black. */
void BlankRegion(vtkImageData* const image, const itk::ImageRegion<2>& region);

/** Set the alpha channel of all pixels in 'pixels' to 'value'. */
void SetPixelTransparency(vtkImageData* const image, const std::vector<itk::Index<2> >& pixels, const unsigned char value);

/** Convert a scalar ITK image into a VTK image after scaling the magnitude to a grayscale range (0 - 255). */
template <typename TImage>
void ITKScalarImageToScaledVTKImage(const TImage* const image, vtkImageData* const outputImage);

/** Create a VTK image of a patch of an image. */
template <typename TImage>
void CreatePatchVTKImage(const TImage* image, const itk::ImageRegion<2>& region, vtkImageData* outputImage);

} // end namespace

#include "ITKVTKHelpers.hpp"

#endif
