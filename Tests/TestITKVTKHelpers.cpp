#include "ITKVTKHelpers.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>

// ITK
#include "itkImageRegion.h"

void TestOutlineRegion();

int main( int argc, char ** argv )
{
  TestOutlineRegion();
  return 0;
}

void TestOutlineRegion()
{
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  image->SetDimensions(10,10,1);
  image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
  
  itk::Index<2> corner = {{3,3}};
  itk::Size<2> size = {{3,3}};
  itk::ImageRegion<2> region(corner, size);

  unsigned char red[3] = {255,0,0};
  ITKVTKHelpers::OutlineRegion(image, region, red);
  
  vtkSmartPointer<vtkPNGWriter> writer =
    vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName("test.png");
  writer->SetInputData(image);
  writer->Write();
}
