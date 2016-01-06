#pragma once

#include "image_processing.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace WoundSegmentation {

  const size_t PICTURE_BOX_LOCATION_GAP = 50;

  ImageProcessing::ImageProcesser image_processer;

  using namespace System;
  using namespace System::ComponentModel;
  using namespace System::Collections;
  using namespace System::Windows::Forms;
  using namespace System::Data;
  using namespace System::Drawing;

  /// <summary>
  /// Summary for MyForm
  /// </summary>
  public ref class WoundSegmentationForm : public System::Windows::Forms::Form {
  public:
    WoundSegmentationForm(void) {
      InitializeComponent();
      //
      //TODO: Add the constructor code here
      //

      open_file_tool_strip_menu_item_->Click += gcnew System::EventHandler(this, &WoundSegmentation::WoundSegmentationForm::OnButtonsClick);
      save_file_tool_strip_menu_item_->Click += gcnew System::EventHandler(this, &WoundSegmentation::WoundSegmentationForm::OnButtonsClick);

      button_segmentation_->Click += gcnew System::EventHandler(this, &WoundSegmentation::WoundSegmentationForm::OnButtonsClick);
      button_confidence_->Click += gcnew System::EventHandler(this, &WoundSegmentation::WoundSegmentationForm::OnButtonsClick);
      button_significance_->Click += gcnew System::EventHandler(this, &WoundSegmentation::WoundSegmentationForm::OnButtonsClick);
      button_wound_segmentation_->Click += gcnew System::EventHandler(this, &WoundSegmentation::WoundSegmentationForm::OnButtonsClick);
      button_evaluation_->Click += gcnew System::EventHandler(this, &WoundSegmentation::WoundSegmentationForm::OnButtonsClick);

      RunAllData();
    }

  protected:
    /// <summary>
    /// Clean up any resources being used.
    /// </summary>
    ~WoundSegmentationForm() {
      if (components) {
        delete components;
      }
    }

  private:

    void RunAllData() {

      std::ofstream dc_output_stream;

      dc_output_stream.open("../data/dc.txt");

      double total_dc = 0;
      size_t ground_truth_count = 0;

      for (size_t picture_id = 1; picture_id <= 7; ++picture_id) {
        std::ostringstream input_picture_name_oss;
        input_picture_name_oss << "../data/S" << std::setw(2) << std::setfill('0') << picture_id << ".jpg";

        std::cout << "Start : " << input_picture_name_oss.str() << "\n";

        System::String ^input_file_path = gcnew System::String(input_picture_name_oss.str().c_str());
        Bitmap ^source_bitmap = gcnew Bitmap(input_file_path);
        image_processer.SetPixelValuesFromBitmap(source_bitmap);

        image_processer.WoundSegmentation();

        Bitmap ^result_bitmap = image_processer.GetResultBitmapFromPixelValues();
        picture_box_result_->Image = result_bitmap;

        std::ostringstream output_picture_name_oss;
        output_picture_name_oss << "../data/result_S" << std::setw(2) << std::setfill('0') << picture_id << ".jpg";

        System::String ^output_file_path = gcnew System::String(output_picture_name_oss.str().c_str());

        picture_box_result_->Image->Save(output_file_path, System::Drawing::Imaging::ImageFormat::Bmp);

        std::cout << "File saved : " << output_picture_name_oss.str() << "\n";

        std::ostringstream ground_truth_picture_name_oss;
        ground_truth_picture_name_oss << "../data/S" << std::setw(2) << std::setfill('0') << picture_id << "_GT.jpg";

        if (std::ifstream(ground_truth_picture_name_oss.str()).good()) {
          System::String ^ground_truth_file_path = gcnew System::String(ground_truth_picture_name_oss.str().c_str());
          image_processer.SetGroundTruthPixelValuesFromBitmap(gcnew Bitmap(ground_truth_file_path));
          double dc = image_processer.DiceCoefficient();

          std::cout << "dc : " << dc << " " << input_picture_name_oss.str() << "\n";
          dc_output_stream << "dc : " << dc << " " << input_picture_name_oss.str() << "\n";

          total_dc += dc;
          ++ground_truth_count;

          //std::ostringstream evaluation_picture_name_oss;
          //evaluation_picture_name_oss << "../data/evaluation_S" << std::setw(2) << std::setfill('0') << picture_id << ".jpg";

          //output_file_path = gcnew System::String(evaluation_picture_name_oss.str().c_str());
          //image_processer.WoundSegmentationWithOutlineOverlapping();
          //result_bitmap = image_processer.GetResultBitmapFromPixelValues();
          //picture_box_result_->Image = result_bitmap;
          //picture_box_result_->Image->Save(output_file_path, System::Drawing::Imaging::ImageFormat::Bmp);
          //std::cout << "File saved : " << evaluation_picture_name_oss.str() << "\n";
        }
      }

      if (ground_truth_count) {
        dc_output_stream << "average dc : " << total_dc / (double)ground_truth_count << "\n";
      }

      dc_output_stream.close();

      exit(0);
    }

    void OnButtonsClick(System::Object ^sender, System::EventArgs ^e) {
      if (sender == open_file_tool_strip_menu_item_) {
        OpenFileDialog ^open_image_file_dialog = gcnew OpenFileDialog();
        open_image_file_dialog->Filter = "Image files | *.*";
        open_image_file_dialog->Title = "Open an image file.";

        if (open_image_file_dialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
          Bitmap ^source_bitmap = gcnew Bitmap(open_image_file_dialog->FileName);

          image_processer.SetPixelValuesFromBitmap(source_bitmap);

          picture_box_source_->Image = image_processer.GetSourceBitmapFromPixelValues();

          picture_box_result_->Image = image_processer.GetResultBitmapFromPixelValues();

          picture_box_result_->Location = System::Drawing::Point(picture_box_source_->Location.X + source_bitmap->Width + PICTURE_BOX_LOCATION_GAP, picture_box_source_->Location.Y);

          label_original_image_->Location = System::Drawing::Point(picture_box_source_->Location.X + source_bitmap->Width * 0.5, label_original_image_->Location.Y);

          label_result_image_->Location = System::Drawing::Point(picture_box_result_->Location.X + source_bitmap->Width * 0.5, label_original_image_->Location.Y);
        }

      } else if (sender == save_file_tool_strip_menu_item_) {
        SaveFileDialog ^save_bmp_file_dialog = gcnew SaveFileDialog();
        save_bmp_file_dialog->Filter = "BMP image files | *.bmp";
        save_bmp_file_dialog->Title = "Save a BMP image file.";

        if (save_bmp_file_dialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
          System::IO::Stream^ bmp_file_stream;
          if ((bmp_file_stream = save_bmp_file_dialog->OpenFile()) != nullptr && save_bmp_file_dialog->FileName != "") {
            bmp_file_stream->Close();
            picture_box_result_->Image->Save(save_bmp_file_dialog->FileName, System::Drawing::Imaging::ImageFormat::Bmp);
          }
        }
      } else if (sender == button_segmentation_) {
        image_processer.Segmentation();
        picture_box_result_->Image = image_processer.GetResultBitmapFromPixelValues();
      } else if (sender == button_confidence_) {
        image_processer.ConfidenceMap();
        picture_box_result_->Image = image_processer.GetResultBitmapFromPixelValues();
      } else if (sender == button_significance_) {
        image_processer.SignificanceMap();
        picture_box_result_->Image = image_processer.GetResultBitmapFromPixelValues();
      } else if (sender == button_wound_segmentation_) {
        image_processer.WoundSegmentation();
        picture_box_result_->Image = image_processer.GetResultBitmapFromPixelValues();
      } else if (sender == button_evaluation_) {
        OpenFileDialog ^open_image_file_dialog = gcnew OpenFileDialog();
        open_image_file_dialog->Filter = "Image files | *.*";
        open_image_file_dialog->Title = "Open an image file.";

        if (open_image_file_dialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
          image_processer.SetGroundTruthPixelValuesFromBitmap(gcnew Bitmap(open_image_file_dialog->FileName));
          image_processer.WoundSegmentationWithOutlineOverlapping();
          picture_box_result_->Image = image_processer.GetResultBitmapFromPixelValues();

          std::cout << image_processer.DiceCoefficient() << "\n";
        }
      }
    }

    /// <summary>
    /// Required designer variable.
    /// </summary>
    System::ComponentModel::Container ^components;

    System::Windows::Forms::PictureBox ^picture_box_source_;
    System::Windows::Forms::PictureBox ^picture_box_result_;

    System::Windows::Forms::MenuStrip ^menu_strip;
    System::Windows::Forms::ToolStripMenuItem ^file_tool_strip_menu_item_;
    System::Windows::Forms::ToolStripMenuItem ^open_file_tool_strip_menu_item_;
    System::Windows::Forms::ToolStripMenuItem ^save_file_tool_strip_menu_item_;

    System::Windows::Forms::Label ^label_original_image_;
    System::Windows::Forms::Label ^label_result_image_;

    System::Windows::Forms::Button ^button_segmentation_;
    System::Windows::Forms::Button ^button_confidence_;
    System::Windows::Forms::Button ^button_significance_;
    System::Windows::Forms::Button ^button_wound_segmentation_;
    System::Windows::Forms::Button ^button_evaluation_;

    System::Windows::Forms::Panel ^panel_picture_box_;

#pragma region Windows Form Designer generated code
    /// <summary>
    /// Required method for Designer support - do not modify
    /// the contents of this method with the code editor.
    /// </summary>
    void InitializeComponent(void) {
      this->picture_box_source_ = (gcnew System::Windows::Forms::PictureBox());
      this->picture_box_result_ = (gcnew System::Windows::Forms::PictureBox());
      this->menu_strip = (gcnew System::Windows::Forms::MenuStrip());
      this->file_tool_strip_menu_item_ = (gcnew System::Windows::Forms::ToolStripMenuItem());
      this->open_file_tool_strip_menu_item_ = (gcnew System::Windows::Forms::ToolStripMenuItem());
      this->save_file_tool_strip_menu_item_ = (gcnew System::Windows::Forms::ToolStripMenuItem());
      this->label_original_image_ = (gcnew System::Windows::Forms::Label());
      this->button_segmentation_ = (gcnew System::Windows::Forms::Button());
      this->button_confidence_ = (gcnew System::Windows::Forms::Button());
      this->button_significance_ = (gcnew System::Windows::Forms::Button());
      this->panel_picture_box_ = (gcnew System::Windows::Forms::Panel());
      this->label_result_image_ = (gcnew System::Windows::Forms::Label());
      this->button_wound_segmentation_ = (gcnew System::Windows::Forms::Button());
      this->button_evaluation_ = (gcnew System::Windows::Forms::Button());
      (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->picture_box_source_))->BeginInit();
      (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->picture_box_result_))->BeginInit();
      this->menu_strip->SuspendLayout();
      this->panel_picture_box_->SuspendLayout();
      this->SuspendLayout();
      // 
      // picture_box_source_
      // 
      this->picture_box_source_->Location = System::Drawing::Point(0, 0);
      this->picture_box_source_->Name = L"picture_box_source_";
      this->picture_box_source_->Size = System::Drawing::Size(60, 60);
      this->picture_box_source_->SizeMode = System::Windows::Forms::PictureBoxSizeMode::AutoSize;
      this->picture_box_source_->TabIndex = 0;
      this->picture_box_source_->TabStop = false;
      // 
      // picture_box_result_
      // 
      this->picture_box_result_->Location = System::Drawing::Point(60, 0);
      this->picture_box_result_->Name = L"picture_box_result_";
      this->picture_box_result_->Size = System::Drawing::Size(60, 60);
      this->picture_box_result_->SizeMode = System::Windows::Forms::PictureBoxSizeMode::AutoSize;
      this->picture_box_result_->TabIndex = 1;
      this->picture_box_result_->TabStop = false;
      // 
      // menu_strip
      // 
      this->menu_strip->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) {
        this->file_tool_strip_menu_item_
      });
      this->menu_strip->Location = System::Drawing::Point(0, 0);
      this->menu_strip->Name = L"menu_strip";
      this->menu_strip->Size = System::Drawing::Size(1584, 24);
      this->menu_strip->TabIndex = 2;
      this->menu_strip->Text = L"Menu strip";
      // 
      // file_tool_strip_menu_item_
      // 
      this->file_tool_strip_menu_item_->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
        this->open_file_tool_strip_menu_item_,
          this->save_file_tool_strip_menu_item_
      });
      this->file_tool_strip_menu_item_->Name = L"file_tool_strip_menu_item_";
      this->file_tool_strip_menu_item_->Size = System::Drawing::Size(37, 20);
      this->file_tool_strip_menu_item_->Text = L"File";
      // 
      // open_file_tool_strip_menu_item_
      // 
      this->open_file_tool_strip_menu_item_->Name = L"open_file_tool_strip_menu_item_";
      this->open_file_tool_strip_menu_item_->Size = System::Drawing::Size(124, 22);
      this->open_file_tool_strip_menu_item_->Text = L"Open File";
      // 
      // save_file_tool_strip_menu_item_
      // 
      this->save_file_tool_strip_menu_item_->Name = L"save_file_tool_strip_menu_item_";
      this->save_file_tool_strip_menu_item_->Size = System::Drawing::Size(124, 22);
      this->save_file_tool_strip_menu_item_->Text = L"Save File";
      // 
      // label_original_image_
      // 
      this->label_original_image_->AutoSize = true;
      this->label_original_image_->Location = System::Drawing::Point(120, 30);
      this->label_original_image_->Name = L"label_original_image_";
      this->label_original_image_->Size = System::Drawing::Size(74, 13);
      this->label_original_image_->TabIndex = 3;
      this->label_original_image_->Text = L"Original Image";
      // 
      // button_segmentation_
      // 
      this->button_segmentation_->Location = System::Drawing::Point(10, 100);
      this->button_segmentation_->Name = L"button_segmentation_";
      this->button_segmentation_->Size = System::Drawing::Size(125, 25);
      this->button_segmentation_->TabIndex = 5;
      this->button_segmentation_->Text = L"Segmentation";
      this->button_segmentation_->UseVisualStyleBackColor = true;
      // 
      // button_confidence_
      // 
      this->button_confidence_->Location = System::Drawing::Point(10, 150);
      this->button_confidence_->Name = L"button_confidence_";
      this->button_confidence_->Size = System::Drawing::Size(125, 25);
      this->button_confidence_->TabIndex = 6;
      this->button_confidence_->Text = L"Confidence";
      this->button_confidence_->UseVisualStyleBackColor = true;
      // 
      // button_significance_
      // 
      this->button_significance_->Location = System::Drawing::Point(10, 200);
      this->button_significance_->Name = L"button_significance_";
      this->button_significance_->Size = System::Drawing::Size(125, 25);
      this->button_significance_->TabIndex = 7;
      this->button_significance_->Text = L"Significance";
      this->button_significance_->UseVisualStyleBackColor = true;
      // 
      // panel_picture_box_
      // 
      this->panel_picture_box_->AutoScroll = true;
      this->panel_picture_box_->Controls->Add(this->label_original_image_);
      this->panel_picture_box_->Controls->Add(this->label_result_image_);
      this->panel_picture_box_->Controls->Add(this->picture_box_result_);
      this->panel_picture_box_->Controls->Add(this->picture_box_source_);
      this->panel_picture_box_->Location = System::Drawing::Point(150, 24);
      this->panel_picture_box_->Name = L"panel_picture_box_";
      this->panel_picture_box_->Size = System::Drawing::Size(1366, 768);
      this->panel_picture_box_->TabIndex = 8;
      // 
      // label_result_image_
      // 
      this->label_result_image_->AutoSize = true;
      this->label_result_image_->Location = System::Drawing::Point(620, 30);
      this->label_result_image_->Name = L"label_result_image_";
      this->label_result_image_->Size = System::Drawing::Size(69, 13);
      this->label_result_image_->TabIndex = 4;
      this->label_result_image_->Text = L"Result Image";
      // 
      // button_wound_segmentation_
      // 
      this->button_wound_segmentation_->Location = System::Drawing::Point(10, 250);
      this->button_wound_segmentation_->Name = L"button_wound_segmentation_";
      this->button_wound_segmentation_->Size = System::Drawing::Size(125, 25);
      this->button_wound_segmentation_->TabIndex = 9;
      this->button_wound_segmentation_->Text = L"Wound Segmentation";
      this->button_wound_segmentation_->UseVisualStyleBackColor = true;
      // 
      // button_evaluation_
      // 
      this->button_evaluation_->Location = System::Drawing::Point(10, 300);
      this->button_evaluation_->Name = L"button_evaluation_";
      this->button_evaluation_->Size = System::Drawing::Size(125, 25);
      this->button_evaluation_->TabIndex = 10;
      this->button_evaluation_->Text = L"Evaluation";
      this->button_evaluation_->UseVisualStyleBackColor = true;
      // 
      // WoundSegmentationForm
      // 
      this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
      this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
      this->AutoScroll = true;
      this->ClientSize = System::Drawing::Size(1584, 862);
      this->Controls->Add(this->button_evaluation_);
      this->Controls->Add(this->button_wound_segmentation_);
      this->Controls->Add(this->panel_picture_box_);
      this->Controls->Add(this->button_confidence_);
      this->Controls->Add(this->button_segmentation_);
      this->Controls->Add(this->button_significance_);
      this->Controls->Add(this->menu_strip);
      this->MainMenuStrip = this->menu_strip;
      this->Name = L"WoundSegmentationForm";
      this->Text = L"Wound Segmentation Form";
      (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->picture_box_source_))->EndInit();
      (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->picture_box_result_))->EndInit();
      this->menu_strip->ResumeLayout(false);
      this->menu_strip->PerformLayout();
      this->panel_picture_box_->ResumeLayout(false);
      this->panel_picture_box_->PerformLayout();
      this->ResumeLayout(false);
      this->PerformLayout();

    }
#pragma endregion
  };
}
