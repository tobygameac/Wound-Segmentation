#pragma once

#include <algorithm>
#include <map>
#include <memory>
#include <iostream>
#include <queue>

#include "graph.h"
#include "disjoint_set.h"

namespace WoundSegmentation {

  namespace ImageProcessing {

    template<class T>
    class RGBA {

    public:

      RGBA() : r_(0), g_(0), b_(0), a_(1) {
      }

      RGBA(T r, T g, T b) : r_(r), g_(g), b_(b), a_(1) {
      }

      RGBA(T r, T g, T b, T a) : RGBA(r, g, b), a_(a) {
      }

      RGBA<T> Square() {
        return RGBA<T>(r_ * r_, g_ * g_, b_ * b_, a_);
      }

      RGBA<T> Sqrt() {
        return RGBA<T>(sqrt(r_), sqrt(g_), sqrt(b_), a_);
      }

      T Magnitude() {
        return sqrt(r_ * r_ + g_ * g_ + b_ * b_);
      }

      const RGBA<T> &operator =(const RGBA<T> &other) {
        r_ = other.r_;
        g_ = other.g_;
        b_ = other.b_;
        a_ = other.a_;

        return *this;
      }

      bool operator ==(const RGBA<T> &other) {
        if (r_ != other.r_) {
          return false;
        }

        if (g_ != other.g_) {
          return false;
        }

        if (b_ != other.b_) {
          return false;
        }

        if (a_ != other.a_) {
          return false;
        }

        return true;
      }

      T r_, g_, b_, a_;
    };

    RGBA<unsigned char> SignificanceValueToSignifanceColor(double significance_value) {
      RGBA<unsigned char> signifance_color(0, 0, 0);

      if (significance_value > 1) {
        signifance_color.r_ = 1;
      }

      if (significance_value < 0) {
        signifance_color.b_ = 1;
      }

      if (significance_value < (1 / 3.0)) {
        signifance_color.g_ = (significance_value * 3.0) * 255;
        signifance_color.b_ = (1 - significance_value * 3.0) * 255;
      } else if (significance_value < (2 / 3.0)) {
        signifance_color.r_ = ((significance_value - (1 / 3.0)) * 3.0) * 255;
        signifance_color.g_ = 1.0 * 255;
      } else if (significance_value <= 1) {
        signifance_color.r_ = 1.0 * 255;
        signifance_color.g_ = (1.0 - (significance_value - (2 / 3.0)) * 3.0) * 255;
      }

      return signifance_color;
    }

    void BitmapToVectorOfPixels(System::Drawing::Bitmap ^source_image, std::vector<std::vector<RGBA<unsigned char> > > &source_image_pixel_values) {
      System::Drawing::Rectangle image_rectangle(0, 0, source_image->Width, source_image->Height);
      System::Drawing::Imaging::BitmapData ^source_image_data = source_image->LockBits(image_rectangle, System::Drawing::Imaging::ImageLockMode::ReadOnly, source_image->PixelFormat);

      source_image_pixel_values = std::vector<std::vector<RGBA<unsigned char> > >(source_image->Height, std::vector<RGBA<unsigned char> >(source_image->Width));

      unsigned char *base_pointer = (unsigned char *)(void *)source_image_data->Scan0;
      size_t bytes_per_pixel = source_image_data->Stride / source_image_data->Width;

      for (size_t r = 0; r < source_image_data->Height; ++r) {
        unsigned char *row_base_pointer = base_pointer + r * source_image_data->Stride;
        for (size_t c = 0; c < source_image_data->Width; ++c) {
          unsigned char *pixel_pointer = row_base_pointer + c * bytes_per_pixel;
          source_image_pixel_values[r][c].r_ = pixel_pointer[2];
          source_image_pixel_values[r][c].g_ = pixel_pointer[1];
          source_image_pixel_values[r][c].b_ = pixel_pointer[0];
        }
      }

      source_image->UnlockBits(source_image_data);
    }

    void VectorOfPixelsToBitmap(const std::vector<std::vector<RGBA<unsigned char> > > &source_image_pixel_values, System::Drawing::Bitmap ^destination_image) {

      if (!source_image_pixel_values.size()) {
        return;
      }

      if (source_image_pixel_values.size() != destination_image->Height || source_image_pixel_values[0].size() != destination_image->Width) {
        return;
      }

      System::Drawing::Rectangle image_rectangle(0, 0, destination_image->Width, destination_image->Height);
      System::Drawing::Imaging::BitmapData ^destination_image_data = destination_image->LockBits(image_rectangle, System::Drawing::Imaging::ImageLockMode::ReadWrite, destination_image->PixelFormat);

      unsigned char *base_pointer = (unsigned char *)(void *)destination_image_data->Scan0;
      size_t bytes_per_pixel = destination_image_data->Stride / destination_image_data->Width;

      for (size_t r = 0; r < destination_image_data->Height; ++r) {
        unsigned char *row_base_pointer = base_pointer + r * destination_image_data->Stride;
        for (size_t c = 0; c < destination_image_data->Width; ++c) {
          unsigned char *pixel_pointer = row_base_pointer + c * bytes_per_pixel;
          pixel_pointer[0] = source_image_pixel_values[r][c].b_;
          pixel_pointer[1] = source_image_pixel_values[r][c].g_;
          pixel_pointer[2] = source_image_pixel_values[r][c].r_;
        }
      }

      destination_image->UnlockBits(destination_image_data);
    }

    double DiceCoefficient(System::Drawing::Bitmap ^image_1, System::Drawing::Bitmap ^image_2) {
      if (!image_1 || !image_2) {
        return -1;
      }

      if (image_1->Width != image_2->Width) {
        return -1;
      }

      if (image_1->Height != image_2->Height) {
        return -1;
      }

      if (!image_1->Width || !image_1->Height) {
        return -1;
      }

      std::vector<std::vector<RGBA<unsigned char> > > image_1_pixel_values;
      BitmapToVectorOfPixels(image_1, image_1_pixel_values);

      std::vector<std::vector<RGBA<unsigned char> > > image_2_pixel_values;
      BitmapToVectorOfPixels(image_2, image_2_pixel_values);

      size_t intersection_pixel_count = 0;
      for (size_t r = 0; r < image_1_pixel_values.size(); ++r) {
        for (size_t c = 0; c < image_1_pixel_values[r].size(); ++c) {
          intersection_pixel_count += image_1_pixel_values[r][c] == image_2_pixel_values[r][c];
        }
      }

      return (2 * intersection_pixel_count) / (double)(2 * image_1->Width * image_1->Height);
    }

    System::Drawing::Bitmap ^RemoveIsolatedRegion(System::Drawing::Bitmap ^source_image) {
      if (!source_image) {
        return nullptr;
      }

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      std::vector<std::vector<size_t> > group_of_pixel = std::vector<std::vector<size_t> >(source_image->Height, std::vector<size_t>(source_image->Width, 0));
      std::vector<std::vector<size_t> > visited = std::vector<std::vector<size_t> >(source_image->Height, std::vector<size_t>(source_image->Width, 0));

      size_t group_count = 0;
      size_t max_group_id = 0;
      size_t max_group_size = 0;

      for (size_t row = 0; row < source_image->Height; ++row) {
        for (size_t column = 0; column < source_image->Width; ++column) {
          if (!source_image_pixel_values[row][column].r_ && !visited[row][column]) {
            visited[row][column] = 1;

            std::queue<std::pair<size_t, size_t> > q;
            q.push(std::make_pair(row, column));

            size_t group_size = 0;

            while (!q.empty()) {
              ++group_size;

              size_t current_r = q.front().first;
              size_t current_c = q.front().second;
              group_of_pixel[current_r][current_c] = group_count;

              q.pop();

              for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                  if (!dr && !dc) {
                    continue;
                  }
                  int neighbor_r = (int)current_r + dr;
                  int neighbor_c = (int)current_c + dc;
                  if (neighbor_r >= 0 && neighbor_r < source_image->Height && neighbor_c >= 0 && neighbor_c < source_image->Width) {
                    if ((source_image_pixel_values[current_r][current_c] == source_image_pixel_values[neighbor_r][neighbor_c]) && !visited[neighbor_r][neighbor_c]) {
                      visited[neighbor_r][neighbor_c] = 1;
                      q.push(std::make_pair(neighbor_r, neighbor_c));
                    }
                  }
                }
              }
            }

            if (group_size > max_group_size) {
              max_group_size = group_size;
              max_group_id = group_count;
            }

            ++group_count;
          }
        }
      }


      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;

      for (size_t row = 0; row < source_image->Height; ++row) {
        for (size_t column = 0; column < source_image->Width; ++column) {

          if (!source_image_pixel_values[row][column].r_ && group_of_pixel[row][column] != max_group_id) {
            result_image_pixel_values[row][column] = RGBA<unsigned char>(255, 255, 255);
          }

        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      return result_image;
    }

    System::Drawing::Bitmap ^MeanFilter(System::Drawing::Bitmap ^source_image, int filter_size) {
      if (!source_image) {
        return nullptr;
      }

      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;

      for (size_t row = 0; row < source_image->Height; ++row) {
        for (size_t column = 0; column < source_image->Width; ++column) {
          size_t total_pixel_count = 0;
          RGBA<size_t> rgb_sum(0, 0, 0);
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int target_r = row + delta_r;
              int target_c = column + delta_c;
              if (target_r >= 0 && target_r < source_image->Height && target_c >= 0 && target_c < source_image->Width) {
                ++total_pixel_count;
                RGBA<unsigned char> rgb_value = source_image_pixel_values[target_r][target_c];
                rgb_sum.r_ += rgb_value.r_;
                rgb_sum.g_ += rgb_value.g_;
                rgb_sum.b_ += rgb_value.b_;
              }
            }
          }
          if (total_pixel_count) {
            result_image_pixel_values[row][column] = RGBA<unsigned char>(rgb_sum.r_ / (float)total_pixel_count, rgb_sum.g_ / (float)total_pixel_count, rgb_sum.b_ / (float)total_pixel_count);
          }
        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      return result_image;
    }

    System::Drawing::Bitmap ^MeanFilter(System::Drawing::Bitmap ^source_image) {
      return MeanFilter(source_image, 3);
    }

    System::Drawing::Bitmap ^MaxFilter(System::Drawing::Bitmap ^source_image, int filter_size) {
      if (!source_image) {
        return nullptr;
      }

      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;

      for (size_t row = 0; row < source_image->Height; ++row) {
        for (size_t column = 0; column < source_image->Width; ++column) {
          RGBA<unsigned char> max_rgb(0, 0, 0);
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int target_r = row + delta_r;
              int target_c = column + delta_c;
              if (target_r >= 0 && target_r < source_image->Height && target_c >= 0 && target_c < source_image->Width) {
                if (source_image_pixel_values[target_r][target_c].Magnitude() > max_rgb.Magnitude()) {
                  max_rgb = source_image_pixel_values[target_r][target_c];
                }
              }
            }
          }

          result_image_pixel_values[row][column] = max_rgb;
        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      return result_image;
    }

    System::Drawing::Bitmap ^MaxFilter(System::Drawing::Bitmap ^source_image) {
      return MaxFilter(source_image, 3);
    }

    System::Drawing::Bitmap ^MinFilter(System::Drawing::Bitmap ^source_image, int filter_size) {
      if (!source_image) {
        return nullptr;
      }

      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;

      for (size_t row = 0; row < source_image->Height; ++row) {
        for (size_t column = 0; column < source_image->Width; ++column) {
          RGBA<unsigned char> min_rgb(255, 255, 255);
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int target_r = row + delta_r;
              int target_c = column + delta_c;
              if (target_r >= 0 && target_r < source_image->Height && target_c >= 0 && target_c < source_image->Width) {
                if (source_image_pixel_values[target_r][target_c].Magnitude() < min_rgb.Magnitude()) {
                  min_rgb = source_image_pixel_values[target_r][target_c];
                }
              }
            }
          }

          result_image_pixel_values[row][column] = min_rgb;
        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      return result_image;
    }

    System::Drawing::Bitmap ^RemoveVerticalLines(System::Drawing::Bitmap ^source_image, size_t max_line_width) {
      if (!source_image) {
        return nullptr;
      }

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;

      for (size_t row = 0; row < source_image->Height; ++row) {
        size_t current_line_width = 0;
        for (size_t column = 0; column < source_image->Width; ++column) {
          if (source_image_pixel_values[row][column].r_) { // White
            if (current_line_width < max_line_width) {
              for (size_t back_index = 1; back_index <= current_line_width; ++back_index) {
                result_image_pixel_values[row][column - back_index] = RGBA<unsigned char>(255, 255, 255);
              }
            }
            current_line_width = 0;
          } else {
            ++current_line_width;
          }
        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      return result_image;
    }

    System::Drawing::Bitmap ^MinFilter(System::Drawing::Bitmap ^source_image) {
      return MinFilter(source_image, 3);
    }

    void BuildGraphFromImage(const std::vector<std::vector<RGBA<unsigned char> > > &pixel_values, Graph<std::pair<size_t, size_t> > &target_graph) {
      target_graph = Graph<std::pair<size_t, size_t> >();

      if (!pixel_values.size()) {
        return;
      }

      const size_t IMAGE_WIDTH = pixel_values[0].size();
      const size_t IMAGE_HEIGHT = pixel_values.size();
      const size_t IMAGE_SIZE = IMAGE_WIDTH * IMAGE_HEIGHT;

      target_graph.vertices_.reserve(IMAGE_SIZE);
      target_graph.edges_.reserve(IMAGE_SIZE * 4);

      const int DELTA_R[4] = {0, 1, 1, -1};
      const int DELTA_C[4] = {1, 0, 1, 1};

      for (int r = 0; r < IMAGE_HEIGHT; ++r) {
        for (int c = 0; c < IMAGE_WIDTH; ++c) {
          target_graph.vertices_.push_back(std::make_pair(c, r));

          size_t index = r * IMAGE_WIDTH + c;

          for (size_t direction = 0; direction < 4; ++direction) {
            int neighbor_r = r + DELTA_R[direction];
            int neighbor_c = c + DELTA_C[direction];
            if (neighbor_r >= 0 && neighbor_r < IMAGE_HEIGHT && neighbor_c >= 0 && neighbor_c < IMAGE_WIDTH) {
              size_t neighbor_index = neighbor_r * IMAGE_WIDTH + neighbor_c;
              std::pair<size_t, size_t> e(index, neighbor_index);

              RGBA<double> edge_difference;
              edge_difference.r_ = (double)pixel_values[r][c].r_ - (double)pixel_values[neighbor_r][neighbor_c].r_;
              edge_difference.g_ = (double)pixel_values[r][c].g_ - (double)pixel_values[neighbor_r][neighbor_c].g_;
              edge_difference.b_ = (double)pixel_values[r][c].b_ - (double)pixel_values[neighbor_r][neighbor_c].b_;

              double edge_weight = edge_difference.Magnitude();

              target_graph.edges_.push_back(Edge(e, edge_weight));
            }
          }
        }
      }
    }

    System::Drawing::Bitmap ^Segmentation(System::Drawing::Bitmap ^source_image, Graph<std::pair<size_t, size_t> > &target_graph, std::vector<std::vector<size_t> > &group_of_pixel, const double k, const int min_patch_size, const double similar_color_patch_merge_threshold) {

      if (!source_image) {
        return nullptr;
      }

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      BuildGraphFromImage(source_image_pixel_values, target_graph);

      sort(target_graph.edges_.begin(), target_graph.edges_.end());

      DisjointSet vertex_disjoint_set(target_graph.vertices_.size());

      std::vector<double> thresholds(target_graph.vertices_.size());
      for (auto &threshold : thresholds) {
        threshold = 1 / k;
      }

      // Segmentation
      for (const auto &edge : target_graph.edges_) {
        size_t group_of_x = vertex_disjoint_set.FindGroup(edge.edge_indices_pair_.first);
        size_t group_of_y = vertex_disjoint_set.FindGroup(edge.edge_indices_pair_.second);
        if (group_of_x == group_of_y) {
          continue;
        }
        if (edge.weight_ <= std::min(thresholds[group_of_x], thresholds[group_of_y])) {
          vertex_disjoint_set.UnionGroup(group_of_x, group_of_y);
          thresholds[group_of_x] = edge.weight_ + (k / vertex_disjoint_set.SizeOfGroup(edge.edge_indices_pair_.first));
        }
      }

      // Deal with the smaller set
      for (const auto &edge : target_graph.edges_) {
        size_t group_of_x = vertex_disjoint_set.FindGroup(edge.edge_indices_pair_.first);
        size_t group_of_y = vertex_disjoint_set.FindGroup(edge.edge_indices_pair_.second);
        if (group_of_x == group_of_y) {
          continue;
        }
        if (min_patch_size > std::min(vertex_disjoint_set.SizeOfGroup(group_of_x), vertex_disjoint_set.SizeOfGroup(group_of_y))) {
          vertex_disjoint_set.UnionGroup(group_of_x, group_of_y);
        }
      }

      // Calculate the color of each group
      std::vector<RGBA<double> > group_color(target_graph.vertices_.size());
      for (size_t r = 0; r < source_image->Height; ++r) {
        for (size_t c = 0; c < source_image->Width; ++c) {
          size_t index = r * source_image->Width + c;
          size_t group = vertex_disjoint_set.FindGroup(index);
          size_t group_size = vertex_disjoint_set.SizeOfGroup(group);
          if (group_size) {
            group_color[group].r_ += source_image_pixel_values[r][c].r_ / (double)group_size;
            group_color[group].g_ += source_image_pixel_values[r][c].g_ / (double)group_size;
            group_color[group].b_ += source_image_pixel_values[r][c].b_ / (double)group_size;
          }
        }
      }

      // Deal with the similar color set
      for (const auto &edge : target_graph.edges_) {
        size_t group_of_x = vertex_disjoint_set.FindGroup(edge.edge_indices_pair_.first);
        size_t group_of_y = vertex_disjoint_set.FindGroup(edge.edge_indices_pair_.second);
        if (group_of_x == group_of_y) {
          continue;
        }
        RGBA<double> difference;

        difference.r_ = group_color[group_of_x].r_ - group_color[group_of_y].r_;
        difference.g_ = group_color[group_of_x].g_ - group_color[group_of_y].g_;
        difference.b_ = group_color[group_of_x].b_ - group_color[group_of_y].b_;

        double color_difference = difference.Magnitude();

        if (color_difference < similar_color_patch_merge_threshold) {
          vertex_disjoint_set.UnionGroup(group_of_x, group_of_y);
        }
      }

      // Calculate the color of each group again
      for (auto &color : group_color) {
        color = RGBA<double>();
      }

      for (size_t r = 0; r < source_image->Height; ++r) {
        for (size_t c = 0; c < source_image->Width; ++c) {
          size_t index = r * source_image->Width + c;
          size_t group = vertex_disjoint_set.FindGroup(index);
          size_t group_size = vertex_disjoint_set.SizeOfGroup(group);
          if (group_size) {
            group_color[group].r_ += source_image_pixel_values[r][c].r_ / (double)group_size;
            group_color[group].g_ += source_image_pixel_values[r][c].g_ / (double)group_size;
            group_color[group].b_ += source_image_pixel_values[r][c].b_ / (double)group_size;
          }
        }
      }

      // Write the pixel value
      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;
      group_of_pixel = std::vector<std::vector<size_t> >(source_image->Height, std::vector<size_t>(source_image->Width));

      for (size_t r = 0; r < source_image->Height; ++r) {
        for (size_t c = 0; c < source_image->Width; ++c) {
          size_t index = r * source_image->Width + c;
          size_t group = vertex_disjoint_set.FindGroup(index);
          group_of_pixel[r][c] = group;
          result_image_pixel_values[r][c] = RGBA<unsigned char>((unsigned char)group_color[group].r_, (unsigned char)group_color[group].g_, (unsigned char)group_color[group].b_);
        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      return result_image;
    }

    System::Drawing::Bitmap ^Segmentation(System::Drawing::Bitmap ^source_image, const double k, const int min_patch_size, const double similar_color_patch_merge_threshold) {
      Graph<std::pair<size_t, size_t> > target_graph;
      std::vector<std::vector<size_t> > group_of_pixel;
      return Segmentation(source_image, target_graph, group_of_pixel, k, min_patch_size, similar_color_patch_merge_threshold);
    }

    System::Drawing::Bitmap ^Segmentation(System::Drawing::Bitmap ^source_image, std::vector<std::vector<size_t> > &group_of_pixel) {
      Graph<std::pair<size_t, size_t> > target_graph;
      const double SEGMENTATION_K = pow(source_image->Width * source_image->Height, 0.6);
      const double SEGMENTATION_MIN_PATCH_SIZE = (source_image->Width * source_image->Height) * 0.001;
      const double SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD = 20;
      return Segmentation(source_image, target_graph, group_of_pixel, SEGMENTATION_K, SEGMENTATION_MIN_PATCH_SIZE, SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD);
    }

    System::Drawing::Bitmap ^Segmentation(System::Drawing::Bitmap ^source_image) {
      Graph<std::pair<size_t, size_t> > target_graph;
      std::vector<std::vector<size_t> > group_of_pixel;
      const double SEGMENTATION_K = pow(source_image->Width * source_image->Height, 0.6);
      const double SEGMENTATION_MIN_PATCH_SIZE = (source_image->Width * source_image->Height) * 0.001;
      const double SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD = 20;
      return Segmentation(source_image, target_graph, group_of_pixel, SEGMENTATION_K, SEGMENTATION_MIN_PATCH_SIZE, SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD);
    }

    System::Drawing::Bitmap ^ConfidenceMap(System::Drawing::Bitmap ^source_image, std::vector<std::vector<double> > &confidence_map) {
      if (!source_image) {
        return nullptr;
      }

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;

      confidence_map = std::vector<std::vector<double> >(source_image->Height, std::vector<double>(source_image->Width));

      const double RED_WEIGHT = 0.9;
      const double POSITION_WEIGHT = 0.1;
      const double SIMILIAR_WEIGHT = 0.0;

      for (size_t r = 0; r < result_image_pixel_values.size(); ++r) {
        for (size_t c = 0; c < result_image_pixel_values[r].size(); ++c) {
          double distance_with_red = 0;
          distance_with_red += pow(0.4 - (source_image_pixel_values[r][c].r_ / 255.0), 2.0);
          distance_with_red += pow(1.0 - (source_image_pixel_values[r][c].g_ / 255.0), 2.0);
          distance_with_red += pow(1.0 - (source_image_pixel_values[r][c].b_ / 255.0), 2.0);
          distance_with_red = sqrt(distance_with_red / 3.0);

          double distance_with_center = sqrt(pow(r - result_image_pixel_values.size() / 2.0, 2.0) + pow(c - result_image_pixel_values[r].size() / 2.0, 2.0));
          distance_with_center /= sqrt(pow(result_image_pixel_values.size() / 2.0, 2.0) + pow(result_image_pixel_values[r].size() / 2.0, 2.0));
          distance_with_center = 1 - distance_with_center;
          distance_with_center *= distance_with_center;

          double total_distance_with_neighbor = 0;

          //size_t neighbor_count = 0;
          //for (int dr = -3; dr <= 3; ++dr) {
          //  for (int dc = -3; dc <= 3; ++dc) {
          //    if (!dr && !dc) {
          //      continue;
          //    }
          //    int neighbor_r = (int)r + dr;
          //    int neighbor_c = (int)c + dr;
          //    if (neighbor_r >= 0 && neighbor_r < result_image_pixel_values.size() && neighbor_c >= 0 && neighbor_c < result_image_pixel_values[r].size()) {
          //      ++neighbor_count;

          //      double distance_with_neighbor = pow((source_image_pixel_values[neighbor_r][neighbor_c].r_ / 255.0) - (source_image_pixel_values[r][c].r_ / 255.0), 2.0);
          //      distance_with_neighbor += pow((source_image_pixel_values[neighbor_r][neighbor_c].g_ / 255.0) - (source_image_pixel_values[r][c].g_ / 255.0), 2.0);
          //      distance_with_neighbor += pow((source_image_pixel_values[neighbor_r][neighbor_c].b_ / 255.0) - (source_image_pixel_values[r][c].b_ / 255.0), 2.0);
          //      distance_with_neighbor = sqrt(distance_with_neighbor / 3.0);

          //      total_distance_with_neighbor += distance_with_neighbor;
          //    }
          //  }
          //}

          //if (neighbor_count) {
          //  total_distance_with_neighbor /= neighbor_count;
          //}

          //total_distance_with_neighbor = 1 - total_distance_with_neighbor;

          double confidence_value = RED_WEIGHT * distance_with_red + POSITION_WEIGHT * distance_with_center + SIMILIAR_WEIGHT * total_distance_with_neighbor;

          confidence_map[r][c] = confidence_value;

          unsigned confidence_gray_level = confidence_map[r][c] * 255;

          result_image_pixel_values[r][c].r_ = confidence_gray_level;
          result_image_pixel_values[r][c].g_ = confidence_gray_level;
          result_image_pixel_values[r][c].b_ = confidence_gray_level;
        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      return result_image;
    }

    System::Drawing::Bitmap ^ConfidenceMap(System::Drawing::Bitmap ^source_image) {
      std::vector<std::vector<double> > confidence_map;
      return ConfidenceMap(source_image, confidence_map);
    }

    System::Drawing::Bitmap ^SignificanceMap(System::Drawing::Bitmap ^source_image) {
      if (!source_image) {
        return nullptr;
      }

      source_image = MeanFilter(source_image, 3);

      std::vector<std::vector<size_t> > group_of_pixel;
      Segmentation(source_image, group_of_pixel);

      std::vector<std::vector<double> > confidence_map;
      ConfidenceMap(source_image, confidence_map);

      std::map<size_t, size_t> size_of_group;
      std::map<size_t, double> significance_of_group;

      std::vector<std::vector<RGBA<unsigned char> > > source_image_pixel_values;
      BitmapToVectorOfPixels(source_image, source_image_pixel_values);

      for (size_t r = 0; r < source_image_pixel_values.size(); ++r) {
        for (size_t c = 0; c < source_image_pixel_values[r].size(); ++c) {
          size_t pixel_group = group_of_pixel[r][c];
          ++size_of_group[pixel_group];
          significance_of_group[pixel_group] += confidence_map[r][c];
        }
      }

      double min_significance = 2e9;
      double max_significance = -2e9;

      for (auto group_iterator = significance_of_group.begin(); group_iterator != significance_of_group.end(); ++group_iterator) {
        size_t pixel_group = group_iterator->first;
        size_t group_size = size_of_group[pixel_group];

        if (group_size) {
          group_iterator->second /= (double)group_size;

          min_significance = std::min(min_significance, group_iterator->second);
          max_significance = std::max(max_significance, group_iterator->second);
        }
      }

      // Normalized
      for (auto group_iterator = significance_of_group.begin(); group_iterator != significance_of_group.end(); ++group_iterator) {
        group_iterator->second = (group_iterator->second - min_significance) / (max_significance - min_significance);
      }

      std::vector<std::vector<RGBA<unsigned char> > > result_image_pixel_values = source_image_pixel_values;

      for (size_t r = 0; r < result_image_pixel_values.size(); ++r) {
        for (size_t c = 0; c < result_image_pixel_values[r].size(); ++c) {
          size_t pixel_group = group_of_pixel[r][c];

          double group_significance = significance_of_group[pixel_group];

          RGBA<unsigned char> thresholding_color = (group_significance >= 0.8) ? RGBA<unsigned char>(0, 0, 0) : RGBA<unsigned char>(255, 255, 255);
          //result_image_pixel_values[r][c] = RGBA<unsigned char>((unsigned char)(group_significance * 255), (unsigned char)(group_significance * 255), (unsigned char)(group_significance * 255));
          result_image_pixel_values[r][c] = thresholding_color;
          //result_image_pixel_values[r][c] = SignificanceValueToSignifanceColor(group_significance);
        }
      }

      System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_image);

      VectorOfPixelsToBitmap(result_image_pixel_values, result_image);

      for (size_t t = 0; t < 0; ++t) {
        result_image = MaxFilter(result_image, 15);
      }

      for (size_t t = 0; t < 1; ++t) {
        result_image = RemoveVerticalLines(result_image, 20);
        result_image = RemoveIsolatedRegion(result_image);
      }

      for (size_t t = 0; t < 0; ++t) {
        result_image = MinFilter(result_image, 15);
      }

      return result_image;
    }

  }
}