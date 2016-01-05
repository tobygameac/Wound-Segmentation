#pragma once

#include <algorithm>
#include <map>
#include <memory>
#include <iostream>
#include <queue>

#include "graph.h"
#include "disjoint_set.h"
#include "rgba.h"

namespace WoundSegmentation {

  namespace ImageProcessing {

    class ImageProcesser {

    public:

      ImageProcesser() : image_loaded_(false), ground_truth_image_loaded_(false) {
      }

      ImageProcesser(System::Drawing::Bitmap ^source_bitmap) : ground_truth_image_loaded_(false) {
        SetPixelValuesFromBitmap(source_bitmap);
      }

      void SetPixelValuesFromBitmap(System::Drawing::Bitmap ^source_bitmap) {
        BitmapToVectorOfPixels(source_bitmap, source_pixel_values_);
        result_pixel_values_ = source_pixel_values_;
        image_loaded_ = source_pixel_values_.size() && source_pixel_values_[0].size();
      }

      void SetGroundTruthPixelValuesFromBitmap(System::Drawing::Bitmap ^source_bitmap) {
        BitmapToVectorOfPixels(source_bitmap, ground_truth_pixel_values_);
        ground_truth_image_loaded_ = ground_truth_pixel_values_.size() && ground_truth_pixel_values_[0].size();
      }

      System::Drawing::Bitmap ^GetSourceBitmapFromPixelValues() {
        if (!image_loaded_) {
          return nullptr;
        }

        System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(source_pixel_values_[0].size(), source_pixel_values_.size(), System::Drawing::Imaging::PixelFormat::Format32bppRgb);

        VectorOfPixelsToBitmap(source_pixel_values_, result_image);

        return result_image;
      }

      System::Drawing::Bitmap ^GetResultBitmapFromPixelValues() {
        if (!image_loaded_) {
          return nullptr;
        }

        System::Drawing::Bitmap ^result_image = gcnew System::Drawing::Bitmap(result_pixel_values_[0].size(), result_pixel_values_.size(), System::Drawing::Imaging::PixelFormat::Format32bppRgb);

        VectorOfPixelsToBitmap(result_pixel_values_, result_image);

        return result_image;
      }

      inline bool IsInImage(int r, int c) {
        return image_loaded_ && r >= 0 && r < source_pixel_values_.size() && c >= 0 && c < source_pixel_values_[r].size();
      }

      double DiceCoefficient() {
        if (!image_loaded_ || !ground_truth_image_loaded_) {
          return 0;
        }

        if (result_pixel_values_.size() != ground_truth_pixel_values_.size()) {
          return 0;
        }

        if (result_pixel_values_[0].size() != ground_truth_pixel_values_[0].size()) {
          return 0;
        }

        size_t intersection_pixel_count = 0;
        for (size_t r = 0; r < result_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < result_pixel_values_[r].size(); ++c) {
            intersection_pixel_count += result_pixel_values_[r][c] == ground_truth_pixel_values_[r][c];
          }
        }

        return (2 * intersection_pixel_count) / (double)(2 * ground_truth_pixel_values_.size() * ground_truth_pixel_values_[0].size());
      }

      void RemoveIsolatedRegion() {
        if (!image_loaded_) {
          return;
        }

        std::vector<std::vector<size_t> > group_of_pixel = std::vector<std::vector<size_t> >(source_pixel_values_.size(), std::vector<size_t>(source_pixel_values_[0].size(), 0));
        std::vector<std::vector<size_t> > visited = std::vector<std::vector<size_t> >(source_pixel_values_.size(), std::vector<size_t>(source_pixel_values_[0].size(), 0));

        size_t group_count = 0;
        size_t max_group_id = 0;
        size_t max_group_size = 0;

        for (size_t row = 0; row < source_pixel_values_.size(); ++row) {
          for (size_t column = 0; column < source_pixel_values_[row].size(); ++column) {
            if (!source_pixel_values_[row][column].r_ && !visited[row][column]) {
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
                    if (IsInImage(neighbor_r, neighbor_c)) {
                      if ((source_pixel_values_[current_r][current_c] == source_pixel_values_[neighbor_r][neighbor_c]) && !visited[neighbor_r][neighbor_c]) {
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

        for (size_t row = 0; row < source_pixel_values_.size(); ++row) {
          for (size_t column = 0; column < source_pixel_values_[row].size(); ++column) {

            if (!source_pixel_values_[row][column].r_ && group_of_pixel[row][column] != max_group_id) {
              result_pixel_values_[row][column] = RGBA<unsigned char>(255, 255, 255);
            }

          }
        }
      }

      void MeanFilter(int filter_size) {
        if (!image_loaded_) {
          return;
        }

        if (!(filter_size & 1)) {
          // Only accept odd size filter
          ++filter_size;
        }

        for (size_t row = 0; row < source_pixel_values_.size(); ++row) {
          for (size_t column = 0; column < source_pixel_values_[row].size(); ++column) {
            size_t total_pixel_count = 0;
            RGBA<size_t> rgb_sum(0, 0, 0);
            for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
              for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
                int neighbor_r = (int)row + delta_r;
                int neighbor_c = (int)column + delta_c;
                if (IsInImage(neighbor_r, neighbor_c)) {
                  ++total_pixel_count;
                  RGBA<unsigned char> rgb_value = source_pixel_values_[neighbor_r][neighbor_c];
                  rgb_sum.r_ += rgb_value.r_;
                  rgb_sum.g_ += rgb_value.g_;
                  rgb_sum.b_ += rgb_value.b_;
                }
              }
            }
            if (total_pixel_count) {
              result_pixel_values_[row][column] = RGBA<unsigned char>(rgb_sum.r_ / (float)total_pixel_count, rgb_sum.g_ / (float)total_pixel_count, rgb_sum.b_ / (float)total_pixel_count);
            }
          }
        }
      }

      void MeanFilter() {
        MeanFilter(3);
      }

      void MaxFilter(int filter_size) {
        if (!image_loaded_) {
          return;
        }

        if (!(filter_size & 1)) {
          // Only accept odd size filter
          ++filter_size;
        }

        for (size_t row = 0; row < source_pixel_values_.size(); ++row) {
          for (size_t column = 0; column < source_pixel_values_[row].size(); ++column) {
            RGBA<unsigned char> max_rgb(0, 0, 0);
            for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
              for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
                int target_r = row + delta_r;
                int target_c = column + delta_c;
                if (IsInImage(target_r, target_c)) {
                  if (source_pixel_values_[target_r][target_c].Magnitude() > max_rgb.Magnitude()) {
                    max_rgb = source_pixel_values_[target_r][target_c];
                  }
                }
              }
            }

            result_pixel_values_[row][column] = max_rgb;
          }
        }
      }

      void MaxFilter() {
        MaxFilter(3);
      }

      void MinFilter(int filter_size) {
        if (!image_loaded_) {
          return;
        }

        if (!(filter_size & 1)) {
          // Only accept odd size filter
          ++filter_size;
        }

        for (size_t row = 0; row < source_pixel_values_.size(); ++row) {
          for (size_t column = 0; column < source_pixel_values_[row].size(); ++column) {
            RGBA<unsigned char> min_rgb(255, 255, 255);
            for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
              for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
                int target_r = row + delta_r;
                int target_c = column + delta_c;
                if (IsInImage(target_r, target_c)) {
                  if (source_pixel_values_[target_r][target_c].Magnitude() < min_rgb.Magnitude()) {
                    min_rgb = source_pixel_values_[target_r][target_c];
                  }
                }
              }
            }

            result_pixel_values_[row][column] = min_rgb;
          }
        }
      }

      void MinFilter() {
        MinFilter(3);
      }

      void RemoveVerticalLines(size_t max_line_width) {
        if (!image_loaded_) {
          return;
        }

        for (size_t row = 0; row < source_pixel_values_.size(); ++row) {
          size_t current_line_width = 0;
          for (size_t column = 0; column < source_pixel_values_[row].size(); ++column) {
            if (source_pixel_values_[row][column].r_) { // White
              if (current_line_width < max_line_width) {
                for (size_t back_index = 1; back_index <= current_line_width; ++back_index) {
                  result_pixel_values_[row][column - back_index] = RGBA<unsigned char>(255, 255, 255);
                }
              }
              current_line_width = 0;
            } else {
              ++current_line_width;
            }
          }
        }
      }

      void Segmentation(Graph<std::pair<size_t, size_t> > &target_graph, std::vector<std::vector<size_t> > &group_of_pixel, const double k, const int min_patch_size, const double similar_color_patch_merge_threshold) {

        if (!image_loaded_) {
          return;
        }

        BuildGraphFromImage(source_pixel_values_, target_graph);

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
        for (size_t r = 0; r < source_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < source_pixel_values_[r].size(); ++c) {
            size_t index = r * source_pixel_values_[r].size() + c;
            size_t group = vertex_disjoint_set.FindGroup(index);
            size_t group_size = vertex_disjoint_set.SizeOfGroup(group);
            if (group_size) {
              group_color[group].r_ += source_pixel_values_[r][c].r_ / (double)group_size;
              group_color[group].g_ += source_pixel_values_[r][c].g_ / (double)group_size;
              group_color[group].b_ += source_pixel_values_[r][c].b_ / (double)group_size;
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

        for (size_t r = 0; r < source_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < source_pixel_values_[r].size(); ++c) {
            size_t index = r * source_pixel_values_[r].size() + c;
            size_t group = vertex_disjoint_set.FindGroup(index);
            size_t group_size = vertex_disjoint_set.SizeOfGroup(group);
            if (group_size) {
              group_color[group].r_ += source_pixel_values_[r][c].r_ / (double)group_size;
              group_color[group].g_ += source_pixel_values_[r][c].g_ / (double)group_size;
              group_color[group].b_ += source_pixel_values_[r][c].b_ / (double)group_size;
            }
          }
        }

        // Write the pixel value
        group_of_pixel = std::vector<std::vector<size_t> >(source_pixel_values_.size(), std::vector<size_t>(source_pixel_values_[0].size()));

        for (size_t r = 0; r < source_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < source_pixel_values_[r].size(); ++c) {
            size_t index = r * source_pixel_values_[r].size() + c;
            size_t group = vertex_disjoint_set.FindGroup(index);
            group_of_pixel[r][c] = group;
            result_pixel_values_[r][c] = RGBA<unsigned char>((unsigned char)group_color[group].r_, (unsigned char)group_color[group].g_, (unsigned char)group_color[group].b_);
          }
        }
      }

      void Segmentation(const double k, const int min_patch_size, const double similar_color_patch_merge_threshold) {
        Graph<std::pair<size_t, size_t> > target_graph;
        std::vector<std::vector<size_t> > group_of_pixel;
        Segmentation(target_graph, group_of_pixel, k, min_patch_size, similar_color_patch_merge_threshold);
      }

      void Segmentation(std::vector<std::vector<size_t> > &group_of_pixel) {
        if (!image_loaded_) {
          return;
        }
        Graph<std::pair<size_t, size_t> > target_graph;
        const double SEGMENTATION_K = pow(source_pixel_values_.size() * source_pixel_values_[0].size(), 0.6);
        const double SEGMENTATION_MIN_PATCH_SIZE = (source_pixel_values_.size() * source_pixel_values_[0].size()) * 0.001;
        const double SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD = 20;
        Segmentation(target_graph, group_of_pixel, SEGMENTATION_K, SEGMENTATION_MIN_PATCH_SIZE, SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD);
      }

      void Segmentation() {
        std::vector<std::vector<size_t> > group_of_pixel;
        Segmentation(group_of_pixel);
      }

      void ConfidenceMap(std::vector<std::vector<double> > &confidence_map) {
        if (!image_loaded_) {
          return;
        }

        confidence_map = std::vector<std::vector<double> >(source_pixel_values_.size(), std::vector<double>(source_pixel_values_[0].size()));

        const double RED_WEIGHT = 0.9;
        const double POSITION_WEIGHT = 0.1;
        const double SIMILIAR_WEIGHT = 0.0;

        for (size_t r = 0; r < source_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < source_pixel_values_[r].size(); ++c) {
            double distance_with_red = 0;
            distance_with_red += pow(0.4 - (source_pixel_values_[r][c].r_ / 255.0), 2.0);
            distance_with_red += pow(1.0 - (source_pixel_values_[r][c].g_ / 255.0), 2.0);
            distance_with_red += pow(1.0 - (source_pixel_values_[r][c].b_ / 255.0), 2.0);
            distance_with_red = sqrt(distance_with_red / 3.0);

            double distance_with_center = sqrt(pow(r - source_pixel_values_.size() / 2.0, 2.0) + pow(c - source_pixel_values_[r].size() / 2.0, 2.0));
            distance_with_center /= sqrt(pow(source_pixel_values_.size() / 2.0, 2.0) + pow(source_pixel_values_[r].size() / 2.0, 2.0));
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
            //    if (IsInImage(neighbor_r, neighbor_c)) {
            //      ++neighbor_count;

            //      double distance_with_neighbor = pow((source_pixel_values_[neighbor_r][neighbor_c].r_ / 255.0) - (source_pixel_values_[r][c].r_ / 255.0), 2.0);
            //      distance_with_neighbor += pow((source_pixel_values_[neighbor_r][neighbor_c].g_ / 255.0) - (source_pixel_values_[r][c].g_ / 255.0), 2.0);
            //      distance_with_neighbor += pow((source_pixel_values_[neighbor_r][neighbor_c].b_ / 255.0) - (source_pixel_values_[r][c].b_ / 255.0), 2.0);
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

            result_pixel_values_[r][c].r_ = confidence_gray_level;
            result_pixel_values_[r][c].g_ = confidence_gray_level;
            result_pixel_values_[r][c].b_ = confidence_gray_level;
          }
        }
      }

      void ConfidenceMap() {
        std::vector<std::vector<double> > confidence_map;
        ConfidenceMap(confidence_map);
      }

      void SignificanceMap(std::vector<std::vector<double> > &significance_map) {
        if (!image_loaded_) {
          return;
        }

        significance_map = std::vector<std::vector<double> >(source_pixel_values_.size(), std::vector<double>(source_pixel_values_[0].size()));

        std::vector<std::vector<RGBA<unsigned char> > > original_pixel_values_ = source_pixel_values_;

        MeanFilter(3);
        source_pixel_values_ = result_pixel_values_;

        std::vector<std::vector<size_t> > group_of_pixel;
        Segmentation(group_of_pixel);

        std::vector<std::vector<double> > confidence_map;
        ConfidenceMap(confidence_map);

        std::map<size_t, size_t> size_of_group;
        std::map<size_t, double> significance_of_group;

        for (size_t r = 0; r < source_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < source_pixel_values_[r].size(); ++c) {
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

        for (size_t r = 0; r < source_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < source_pixel_values_[r].size(); ++c) {
            size_t pixel_group = group_of_pixel[r][c];

            double group_significance = significance_of_group[pixel_group];
            significance_map[r][c] = group_significance;

            result_pixel_values_[r][c] = SignificanceValueToSignifanceColor(group_significance);
          }
        }

        source_pixel_values_ = original_pixel_values_;
      }

      void SignificanceMap() {
        std::vector<std::vector<double> > significance_map;
        SignificanceMap(significance_map);
      }

      void WoundSegmentation() {
        if (!image_loaded_) {
          return;
        }

        std::vector<std::vector<RGBA<unsigned char> > > original_pixel_values_ = source_pixel_values_;

        std::vector<std::vector<double> > significance_map;
        SignificanceMap(significance_map);

        for (size_t r = 0; r < source_pixel_values_.size(); ++r) {
          for (size_t c = 0; c < source_pixel_values_[r].size(); ++c) {
            double pixel_significance = significance_map[r][c];

            RGBA<unsigned char> thresholding_color = (pixel_significance >= 0.8) ? RGBA<unsigned char>(0, 0, 0) : RGBA<unsigned char>(255, 255, 255);
            result_pixel_values_[r][c] = thresholding_color;
          }
        }

        for (size_t t = 0; t < 0; ++t) {
          source_pixel_values_ = result_pixel_values_;
          MaxFilter(15);
          source_pixel_values_ = result_pixel_values_;
        }

        for (size_t t = 0; t < 1; ++t) {
          source_pixel_values_ = result_pixel_values_;
          RemoveVerticalLines(20);
          source_pixel_values_ = result_pixel_values_;
          RemoveIsolatedRegion();
          source_pixel_values_ = result_pixel_values_;
        }

        for (size_t t = 0; t < 0; ++t) {
          source_pixel_values_ = result_pixel_values_;
          MinFilter(15);
          source_pixel_values_ = result_pixel_values_;
        }

        source_pixel_values_ = original_pixel_values_;
      }

    private:

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

      bool image_loaded_;
      bool ground_truth_image_loaded_;

      std::vector<std::vector<RGBA<unsigned char> > > source_pixel_values_;
      std::vector<std::vector<RGBA<unsigned char> > > result_pixel_values_;
      std::vector<std::vector<RGBA<unsigned char> > > ground_truth_pixel_values_;
    };
  }
}