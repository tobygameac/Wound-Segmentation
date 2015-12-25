#include "wound_segmentation_form.h"

using namespace System;
using namespace System::Windows::Forms;

[STAThread]
int main(array<String^> ^args) {
  Application::EnableVisualStyles();
  Application::SetCompatibleTextRenderingDefault(false);
  WoundSegmentation::WoundSegmentationForm wound_segmentation_form;
  Application::Run(%wound_segmentation_form);
  return 0;
}