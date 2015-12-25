#pragma once

namespace WoundSegmentation {

  const size_t PICTURE_BOX_LOCATION_GAP = 50;

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

    void OnButtonsClick(System::Object ^sender, System::EventArgs ^e) {
      if (sender == open_file_tool_strip_menu_item_) {
        OpenFileDialog ^open_image_file_dialog = gcnew OpenFileDialog();
        open_image_file_dialog->Filter = "Image files | *.*";
        open_image_file_dialog->Title = "Open an image file.";

        if (open_image_file_dialog->ShowDialog() == System::Windows::Forms::DialogResult::OK) {
          source_image_ = gcnew Bitmap(open_image_file_dialog->FileName);

          picture_box_source_->Image = source_image_;
          picture_box_result_->Image = source_image_;

          picture_box_result_->Location = System::Drawing::Point(picture_box_source_->Location.X + source_image_->Width + PICTURE_BOX_LOCATION_GAP, picture_box_source_->Location.Y);

          label_original_image_->Location = System::Drawing::Point(picture_box_source_->Location.X + source_image_->Width * 0.5, label_original_image_->Location.Y);

          label_result_image_->Location = System::Drawing::Point(picture_box_result_->Location.X + source_image_->Width * 0.5, label_original_image_->Location.Y);
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
      }
    }

    /// <summary>
    /// Required designer variable.
    /// </summary>
    System::ComponentModel::Container ^components;

    System::Drawing::Bitmap ^source_image_;

    System::Windows::Forms::PictureBox ^picture_box_source_;
    System::Windows::Forms::PictureBox ^picture_box_result_;

    System::Windows::Forms::MenuStrip ^menu_strip;
    System::Windows::Forms::ToolStripMenuItem ^file_tool_strip_menu_item_;
    System::Windows::Forms::ToolStripMenuItem ^open_file_tool_strip_menu_item_;
    System::Windows::Forms::ToolStripMenuItem ^save_file_tool_strip_menu_item_;

    System::Windows::Forms::Label ^label_original_image_;
    System::Windows::Forms::Label ^label_result_image_;

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
      this->label_result_image_ = (gcnew System::Windows::Forms::Label());
      this->SuspendLayout();
      // 
      // picture_box_source_
      // 
      this->picture_box_source_->Location = System::Drawing::Point(150, 100);
      this->picture_box_source_->Name = L"picture_box_source_";
      this->picture_box_source_->Size = System::Drawing::Size(600, 600);
      this->picture_box_source_->SizeMode = System::Windows::Forms::PictureBoxSizeMode::AutoSize;
      this->picture_box_source_->TabIndex = 0;
      this->picture_box_source_->TabStop = false;
      // 
      // picture_box_result_
      // 
      this->picture_box_result_->Location = System::Drawing::Point(800, 100);
      this->picture_box_result_->Name = L"picture_box_result_";
      this->picture_box_result_->Size = System::Drawing::Size(600, 600);
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
      this->menu_strip->Size = System::Drawing::Size(1400, 24);
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
      this->label_original_image_->Location = System::Drawing::Point(400, 75);
      this->label_original_image_->Name = L"label_original_image_";
      this->label_original_image_->Size = System::Drawing::Size(74, 13);
      this->label_original_image_->TabIndex = 3;
      this->label_original_image_->Text = L"Original Image";
      // 
      // label_result_image_
      // 
      this->label_result_image_->AutoSize = true;
      this->label_result_image_->Location = System::Drawing::Point(1000, 75);
      this->label_result_image_->Name = L"label_result_image_";
      this->label_result_image_->Size = System::Drawing::Size(69, 13);
      this->label_result_image_->TabIndex = 4;
      this->label_result_image_->Text = L"Result Image";
      // 
      // WoundSegmentationForm
      // 
      this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
      this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
      this->AutoScroll = true;
      this->ClientSize = System::Drawing::Size(1024, 768);
      this->Controls->Add(this->label_result_image_);
      this->Controls->Add(this->label_original_image_);
      this->Controls->Add(this->picture_box_result_);
      this->Controls->Add(this->picture_box_source_);
      this->Controls->Add(this->menu_strip);
      this->MainMenuStrip = this->menu_strip;
      this->Name = L"WoundSegmentationForm";
      this->Text = L"Wound Segmentation Form";
      (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->picture_box_source_))->EndInit();
      (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->picture_box_result_))->EndInit();
      this->menu_strip->ResumeLayout(false);
      this->menu_strip->PerformLayout();
      this->ResumeLayout(false);
      this->PerformLayout();

    }
#pragma endregion
  };
}
