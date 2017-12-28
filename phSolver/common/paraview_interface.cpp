#ifdef USE_CATALYST
#include "vtkCPAdaptorAPI.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"

#include <FCMangle.h>

#define addsurfids FortranCInterface_GLOBAL_(addsurfids,ADDSURFIDS)

extern "C" void addsurfids(int* surfid)
{
  vtkCPInputDataDescription* idd =
    vtkCPAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input");
  vtkUnstructuredGrid* UnstructuredGrid = vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
  if (!UnstructuredGrid)
  {
    vtkGenericWarningMacro("No unstructured grid to attach field data to.");
    return;
  }
  vtkIdType NumberOfNodes = UnstructuredGrid->GetNumberOfPoints();
  if(idd->IsFieldNeeded("surfid"))
  {
  	vtkIntArray* vtksurfid = vtkIntArray::New();
  	vtksurfid->SetName("surfid");
  	vtksurfid->SetArray(surfid,NumberOfNodes, 1);
  	UnstructuredGrid->GetPointData()->AddArray(vtksurfid);
  	vtksurfid->Delete();
  }
}
#endif