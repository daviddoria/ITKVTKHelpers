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

// ITK
#include "itkRegionOfInterestImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"

// VTK
#include <vtkImageData.h>
#include <vtkPolyData.h>

// Submodules
#include "ITKHelpers/ITKHelpers.h"

namespace ITKVTKHelpers
{

typedef itk::Image<unsigned char, 2> UnsignedCharScalarImageType;

template <typename TImage>
void ITKScalarImageToScaledVTKImage(const TImage* const image, vtkImageData* const outputImage)
{
  //std::cout << "ITKScalarImagetoVTKImage()" << std::endl;

  // Rescale and cast for display
  typedef itk::RescaleIntensityImageFilter<TImage, UnsignedCharScalarImageType > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput(image);
  rescaleFilter->Update();

  // Setup and allocate the VTK image
  //outputImage->SetNumberOfScalarComponents(1);
  //outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);
  //outputImage->AllocateScalars();
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<UnsignedCharScalarImageType> imageIterator(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get();

    ++imageIterator;
    }

  outputImage->Modified();
}

template <typename TImage>
void CreatePatchVTKImage(const TImage* image, const itk::ImageRegion<2>& region, vtkImageData* outputImage)
{
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> ExtractFilterType;
  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetRegionOfInterest(region);
  extractFilter->SetInput(image);
  extractFilter->Update();

  ITKVectorImageToVTKImageFromDimension(extractFilter->GetOutput(), outputImage);
}


// Convert a vector ITK image to a VTK image for display
template <typename TPixel>
void ITKVectorImageToVTKImageFromDimension(const itk::VectorImage<TPixel, 2>* const image, vtkImageData* const outputImage)
{
  // If the image has 3 channels, assume it is RGB.
  //if(image->GetNumberOfComponentsPerPixel() == 3)
  if(image->GetNumberOfComponentsPerPixel() >= 3)
    {
    ITKImageToVTKRGBImage(image, outputImage);
    }
  else
    {
    ITKImageToVTKMagnitudeImage(image, outputImage);
    }

  outputImage->Modified();
}

// Convert a vector ITK image to a VTK image for display
template <typename TPixel>
void ITKImageToVTKRGBImage(const itk::VectorImage<TPixel, 2>* const image, vtkImageData* const outputImage)
{
  // This function assumes an ND (with N>3) image has the first 3 channels as RGB and extra
  // information in the remaining channels.

  typedef itk::VectorImage<TPixel, 2> VectorImageType;
  //std::cout << "ITKImagetoVTKRGBImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() < 3)
    {
    std::cerr << "The input image has " << image->GetNumberOfComponentsPerPixel()
              << " components, but at least 3 are required." << std::endl;
    return;
    }

  // Setup and allocate the image data
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<VectorImageType> imageIterator(image, image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    for(unsigned int component = 0; component < 3; component++)
      {
      pixel[component] = static_cast<unsigned char>(imageIterator.Get()[component]);
      }

    ++imageIterator;
    }

  outputImage->Modified();
}


// Convert a vector ITK image to a VTK image for display
template <typename TPixel>
void ITKImageToVTKMagnitudeImage(const itk::VectorImage<TPixel, 2>* const image, vtkImageData* const outputImage)
{
  //std::cout << "ITKImagetoVTKMagnitudeImage()" << std::endl;

  typedef itk::VectorImage<TPixel> VectorImageType;
  
  // Compute the magnitude of the ITK image
  typedef itk::VectorMagnitudeImageFilter<
                  VectorImageType, FloatScalarImageType >  VectorMagnitudeFilterType;

  // Create and setup a magnitude filter
  typename VectorMagnitudeFilterType::Pointer magnitudeFilter = VectorMagnitudeFilterType::New();
  magnitudeFilter->SetInput( image );
  magnitudeFilter->Update();

  // Rescale and cast for display
  typedef itk::RescaleIntensityImageFilter<
                  FloatScalarImageType, UnsignedCharScalarImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput( magnitudeFilter->GetOutput() );
  rescaleFilter->Update();

  // Setup and allocate the VTK image
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<UnsignedCharScalarImageType> imageIterator(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get();
    pixel[1] = imageIterator.Get();
    pixel[2] = imageIterator.Get();

    ++imageIterator;
    }

  outputImage->Modified();
}

template <typename TPixel>
void ITKImageChannelToVTKImage(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel, vtkImageData* const outputImage)
{
  FloatScalarImageType::Pointer channelImage = FloatScalarImageType::New();
  ITKHelpers::ExtractChannel<float>(image, channel, channelImage);
  ITKScalarImageToScaledVTKImage<FloatScalarImageType>(channelImage, outputImage);
}

} // end namespace
