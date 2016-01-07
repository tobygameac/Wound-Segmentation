#pragma once

#include <algorithm>
#include <map>
#include <memory>
#include <iostream>
#include <queue>

#include "graph.h"
#include "disjoint_set.h"
#include "rgba.h"
#include "omp.h"

namespace WoundSegmentation {

  namespace ImageProcessing {

    inline RGBA<unsigned char> SignificanceValueToSignifanceColor(double significance_value) {
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

              double edge_weight = edge_difference.NormalizedMagnitude();

              target_graph.edges_.push_back(Edge(e, edge_weight));
            }
          }
        }
      }
    }

    inline bool IsInImage(const std::vector<std::vector<RGBA<unsigned char> > > &pixel_values, int r, int c) {
      return r >= 0 && r < pixel_values.size() && c >= 0 && c < pixel_values[r].size();
    }

    void RGBToGray(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      for (size_t row = 0; row < source_pixel_values.size(); ++row) {
        for (size_t column = 0; column < source_pixel_values[row].size(); ++column) {
          RGBA<unsigned char> rgb_value = source_pixel_values[row][column];
          float gray_level = (rgb_value.r_ + rgb_value.g_ + rgb_value.b_) / 3.0;
          result_pixel_values[row][column] = RGBA<unsigned char>(gray_level, gray_level, gray_level);
        }
      }
    }

    void Thresholding(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const float threshold) {
      for (size_t row = 0; row < source_pixel_values.size(); ++row) {
        for (size_t column = 0; column < source_pixel_values[row].size(); ++column) {
          RGBA<unsigned char> rgb_value = source_pixel_values[row][column];
          float gray_level = (rgb_value.r_ + rgb_value.g_ + rgb_value.b_) / 3.0;
          result_pixel_values[row][column] = (gray_level >= threshold) ? RGBA<unsigned char>(255, 255, 255) : RGBA<unsigned char>(0, 0, 0);
        }
      }
    }

    double DiceCoefficient(const std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const std::vector<std::vector<RGBA<unsigned char> > > &ground_truth_pixel_values) {
      if (result_pixel_values.size() != ground_truth_pixel_values.size()) {
        return 0;
      }

      if (result_pixel_values[0].size() != ground_truth_pixel_values[0].size()) {
        return 0;
      }

      std::vector<std::vector<RGBA<unsigned char> > > thresholded_result_pixel_values = result_pixel_values;
      std::vector<std::vector<RGBA<unsigned char> > > thresholded_ground_truth_pixel_values = ground_truth_pixel_values;

      Thresholding(result_pixel_values, thresholded_result_pixel_values, 127);
      Thresholding(ground_truth_pixel_values, thresholded_ground_truth_pixel_values, 127);

      size_t intersection_pixel_count = 0;
      size_t total_pixel_count = 0;

      for (size_t r = 0; r < result_pixel_values.size(); ++r) {
        for (size_t c = 0; c < result_pixel_values[r].size(); ++c) {
          // Only check black pixel
          if (!thresholded_result_pixel_values[r][c].r_ || !thresholded_ground_truth_pixel_values[r][c].r_) {
            total_pixel_count += !thresholded_result_pixel_values[r][c].r_;
            total_pixel_count += !thresholded_ground_truth_pixel_values[r][c].r_;
            intersection_pixel_count += thresholded_result_pixel_values[r][c] == thresholded_ground_truth_pixel_values[r][c];
          }
        }
      }

      if (!total_pixel_count) {
        return 0;
      }

      return (2 * intersection_pixel_count) / (double)(total_pixel_count);
    }

    void RemoveIsolatedRegion(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {

      std::vector<std::vector<size_t> > group_of_pixel = std::vector<std::vector<size_t> >(source_pixel_values.size(), std::vector<size_t>(source_pixel_values[0].size(), 0));
      std::vector<std::vector<size_t> > visited = std::vector<std::vector<size_t> >(source_pixel_values.size(), std::vector<size_t>(source_pixel_values[0].size(), 0));

      size_t group_count = 0;
      size_t max_group_id = 0;
      size_t max_group_size = 0;

      for (size_t row = 0; row < source_pixel_values.size(); ++row) {
        for (size_t column = 0; column < source_pixel_values[row].size(); ++column) {
          if (!source_pixel_values[row][column].r_ && !visited[row][column]) {
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
                  if (IsInImage(source_pixel_values, neighbor_r, neighbor_c)) {
                    if ((source_pixel_values[current_r][current_c] == source_pixel_values[neighbor_r][neighbor_c]) && !visited[neighbor_r][neighbor_c]) {
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

#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
#pragma omp parallel for
        for (int column = 0; column < source_pixel_values[row].size(); ++column) {
          if (!source_pixel_values[row][column].r_ && group_of_pixel[row][column] != max_group_id) {
            result_pixel_values[row][column] = RGBA<unsigned char>(255, 255, 255);
          }

        }
      }
    }

    void GaussianFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, int filter_size, double sigma) {

      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

      std::vector<std::vector<double> > filter(filter_size, std::vector<double>(filter_size, 1 / (2 * acos(-1) * sigma * sigma)));
      for (int r = 0; r < filter_size; ++r) {
        for (int c = 0; c < filter_size; ++c) {
          int x = (filter_size / 2) - c;
          int y = (filter_size / 2) - r;
          filter[r][c] *= exp(-(x * x + y * y) / (2.0 * sigma * sigma));
        }
      }

#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
#pragma omp parallel for
        for (int column = 0; column < source_pixel_values[row].size(); ++column) {
          RGBA<double> rgb_sum(0, 0, 0);
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int neighbor_r = (int)row + delta_r;
              int neighbor_c = (int)column + delta_c;
              double filter_value = filter[delta_r + filter_size / 2][delta_c + filter_size / 2];
              if (IsInImage(source_pixel_values, neighbor_r, neighbor_c)) {
                RGBA<unsigned char> rgb_value = source_pixel_values[neighbor_r][neighbor_c];
                rgb_sum.r_ += rgb_value.r_ * filter_value;
                rgb_sum.g_ += rgb_value.g_ * filter_value;
                rgb_sum.b_ += rgb_value.b_ * filter_value;
              }
            }
          }
          result_pixel_values[row][column] = RGBA<unsigned char>(rgb_sum.r_, rgb_sum.g_, rgb_sum.b_);
        }
      }
    }

    void GaussianFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      GaussianFilter(source_pixel_values, result_pixel_values, 3, 0.84089642);
    }

    void MeanFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, int filter_size) {

      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
#pragma omp parallel for
        for (int column = 0; column < source_pixel_values[row].size(); ++column) {
          size_t total_pixel_count = 0;
          RGBA<size_t> rgb_sum(0, 0, 0);
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int neighbor_r = (int)row + delta_r;
              int neighbor_c = (int)column + delta_c;
              if (IsInImage(source_pixel_values, neighbor_r, neighbor_c)) {
                ++total_pixel_count;
                RGBA<unsigned char> rgb_value = source_pixel_values[neighbor_r][neighbor_c];
                rgb_sum.r_ += rgb_value.r_;
                rgb_sum.g_ += rgb_value.g_;
                rgb_sum.b_ += rgb_value.b_;
              }
            }
          }
          if (total_pixel_count) {
            result_pixel_values[row][column] = RGBA<unsigned char>(rgb_sum.r_ / (float)total_pixel_count, rgb_sum.g_ / (float)total_pixel_count, rgb_sum.b_ / (float)total_pixel_count);
          }
        }
      }
    }

    void MeanFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      MeanFilter(source_pixel_values, result_pixel_values, 3);
    }

    void MedianFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, int filter_size) {

      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

      // Gray image only
#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
#pragma omp parallel for
        for (int column = 0; column < source_pixel_values[row].size(); ++column) {
          std::vector<unsigned char> gray_values;
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int target_r = row + delta_r;
              int target_c = column + delta_c;
              if (IsInImage(source_pixel_values, target_r, target_c)) {
                RGBA<unsigned char> rgb_value = source_pixel_values[target_r][target_c];
                gray_values.push_back(((double)rgb_value.r_ + (double)rgb_value.g_ + (double)rgb_value.b_) / 3.0);
              }
            }
          }

          if (gray_values.size()) {
            std::sort(gray_values.begin(), gray_values.end());
            unsigned char gray_value = gray_values[gray_values.size() / 2];
            result_pixel_values[row][column] = RGBA<unsigned char>(gray_value, gray_value, gray_value);
          }
        }
      }
    }

    void MedianFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      MedianFilter(source_pixel_values, result_pixel_values, 3);
    }

    void MaxFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, int filter_size) {
      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
#pragma omp parallel for
        for (int column = 0; column < source_pixel_values[row].size(); ++column) {
          RGBA<unsigned char> max_rgb(0, 0, 0);
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int target_r = row + delta_r;
              int target_c = column + delta_c;
              if (IsInImage(source_pixel_values, target_r, target_c)) {
                if (source_pixel_values[target_r][target_c].Magnitude() > max_rgb.Magnitude()) {
                  max_rgb = source_pixel_values[target_r][target_c];
                }
              }
            }
          }

          result_pixel_values[row][column] = max_rgb;
        }
      }
    }

    void MaxFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      MaxFilter(source_pixel_values, result_pixel_values, 3);
    }

    void MinFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, int filter_size) {
      if (!(filter_size & 1)) {
        // Only accept odd size filter
        ++filter_size;
      }

#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
#pragma omp parallel for
        for (int column = 0; column < source_pixel_values[row].size(); ++column) {
          RGBA<unsigned char> min_rgb(255, 255, 255);
          for (int delta_r = -filter_size / 2; delta_r <= filter_size / 2; ++delta_r) {
            for (int delta_c = -filter_size / 2; delta_c <= filter_size / 2; ++delta_c) {
              int target_r = row + delta_r;
              int target_c = column + delta_c;
              if (IsInImage(source_pixel_values, target_r, target_c)) {
                if (source_pixel_values[target_r][target_c].Magnitude() < min_rgb.Magnitude()) {
                  min_rgb = source_pixel_values[target_r][target_c];
                }
              }
            }
          }

          result_pixel_values[row][column] = min_rgb;
        }
      }
    }

    void MinFilter(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      MinFilter(source_pixel_values, result_pixel_values, 3);
    }

    void RemoveVerticalLines(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const size_t max_line_width) {
#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
        size_t current_line_width = 0;
        for (size_t column = 0; column < source_pixel_values[row].size(); ++column) {
          if (source_pixel_values[row][column].r_) { // White
            if (current_line_width < max_line_width) {
              for (size_t back_index = 1; back_index <= current_line_width; ++back_index) {
                result_pixel_values[row][column - back_index] = RGBA<unsigned char>(255, 255, 255);
              }
            }
            current_line_width = 0;
          } else {
            ++current_line_width;
          }
        }
      }
    }

    void RemoveHorizontalLines(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const size_t max_line_height) {
#pragma omp parallel for
      for (int column = 0; column < source_pixel_values[0].size(); ++column) {
        for (size_t row = 0; row < source_pixel_values.size(); ++row) {
          size_t current_line_height = 0;
          if (source_pixel_values[row][column].r_) { // White
            if (current_line_height < max_line_height) {
              for (size_t back_index = 1; back_index <= current_line_height; ++back_index) {
                result_pixel_values[row - back_index][column] = RGBA<unsigned char>(255, 255, 255);
              }
            }
            current_line_height = 0;
          } else {
            ++current_line_height;
          }
        }
      }
    }

    void FillVerticalGaps(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const size_t max_gap_width) {
#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
        size_t current_gap_width = 0;
        for (size_t column = 0; column < source_pixel_values[row].size(); ++column) {
          if (!source_pixel_values[row][column].r_) { // Black
            if (current_gap_width < max_gap_width) {
              for (size_t back_index = 1; back_index <= current_gap_width; ++back_index) {
                result_pixel_values[row][column - back_index] = RGBA<unsigned char>(0, 0, 0);
              }
            }
            current_gap_width = 0;
          } else {
            ++current_gap_width;
          }
        }
      }
    }

    void FillHorizontalGaps(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const size_t max_gap_height) {
#pragma omp parallel for
      for (int column = 0; column < source_pixel_values[0].size(); ++column) {
        for (size_t row = 0; row < source_pixel_values.size(); ++row) {
          size_t current_gap_height = 0;
          if (!source_pixel_values[row][column].r_) { // Black
            if (current_gap_height < max_gap_height) {
              for (size_t back_index = 1; back_index <= current_gap_height; ++back_index) {
                result_pixel_values[row - back_index][column] = RGBA<unsigned char>(0, 0, 0);
              }
            }
            current_gap_height = 0;
          } else {
            ++current_gap_height;
          }
        }
      }
    }

    void SobelEdgeDetection(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      const static int KERNAL_SIZE = 3;
      const static float VERTICAL_KERNAL[KERNAL_SIZE][KERNAL_SIZE] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
      const static float HORIZONTAL_KERNAL[KERNAL_SIZE][KERNAL_SIZE] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

#pragma omp parallel for
      for (int row = 0; row < source_pixel_values.size(); ++row) {
#pragma omp parallel for
        for (int column = 0; column < source_pixel_values[row].size(); ++column) {
          RGBA<float> rgb_sum_vertical(0, 0, 0);
          RGBA<float> rgb_sum_horizontal(0, 0, 0);
          for (int delta_r = -KERNAL_SIZE / 2; delta_r <= KERNAL_SIZE / 2; ++delta_r) {
            for (int delta_c = -KERNAL_SIZE / 2; delta_c <= KERNAL_SIZE / 2; ++delta_c) {
              int target_r = row + delta_r;
              int target_c = column + delta_c;
              if (IsInImage(source_pixel_values, target_r, target_c)) {
                RGBA<unsigned char> rgb_value = source_pixel_values[target_r][target_c];
                float vertical_weight = VERTICAL_KERNAL[delta_r + KERNAL_SIZE / 2][delta_c + KERNAL_SIZE / 2];
                float horizontal_weight = HORIZONTAL_KERNAL[delta_r + KERNAL_SIZE / 2][delta_c + KERNAL_SIZE / 2];

                rgb_sum_vertical.r_ += rgb_value.r_ * vertical_weight;
                rgb_sum_vertical.g_ += rgb_value.g_ * vertical_weight;
                rgb_sum_vertical.b_ += rgb_value.b_ * vertical_weight;

                rgb_sum_horizontal.r_ += rgb_value.r_ * horizontal_weight;
                rgb_sum_horizontal.g_ += rgb_value.g_ * horizontal_weight;
                rgb_sum_horizontal.b_ += rgb_value.b_ * horizontal_weight;
              }
            }
          }

          rgb_sum_vertical = RGBA<float>(std::abs(rgb_sum_vertical.r_), std::abs(rgb_sum_vertical.g_), std::abs(rgb_sum_vertical.b_));
          rgb_sum_horizontal = RGBA<float>(std::abs(rgb_sum_horizontal.r_), std::abs(rgb_sum_horizontal.g_), std::abs(rgb_sum_horizontal.b_));

          RGBA<float> rgb_sum = RGBA<float>(sqrt(rgb_sum_vertical.r_ * rgb_sum_vertical.r_ + rgb_sum_horizontal.r_ * rgb_sum_horizontal.r_),
            sqrt(rgb_sum_vertical.g_ * rgb_sum_vertical.g_ + rgb_sum_horizontal.g_ * rgb_sum_horizontal.g_),
            sqrt(rgb_sum_vertical.b_ * rgb_sum_vertical.b_ + rgb_sum_horizontal.b_ * rgb_sum_horizontal.b_));

          RGBA<float> real_rgb_value(rgb_sum.r_, rgb_sum.g_, rgb_sum.b_);
          real_rgb_value.r_ = std::max(0.0f, real_rgb_value.r_);
          real_rgb_value.g_ = std::max(0.0f, real_rgb_value.g_);
          real_rgb_value.b_ = std::max(0.0f, real_rgb_value.b_);

          real_rgb_value.r_ = std::min(255.0f, real_rgb_value.r_);
          real_rgb_value.g_ = std::min(255.0f, real_rgb_value.g_);
          real_rgb_value.b_ = std::min(255.0f, real_rgb_value.b_);

          result_pixel_values[row][column] = RGBA<unsigned char>(real_rgb_value.r_, real_rgb_value.g_, real_rgb_value.b_);

        }
      }
    }

    void Segmentation(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, Graph<std::pair<size_t, size_t> > &target_graph, std::vector<std::vector<size_t> > &group_of_pixel, const double k, const int min_patch_size, const double similar_color_patch_merge_threshold) {

      std::vector<std::vector<RGBA<unsigned char> > > intermediate_pixel_values = source_pixel_values;

      //GaussianFilter(source_pixel_values, intermediate_pixel_values);
      //MedianFilter(source_pixel_values, intermediate_pixel_values);
      MeanFilter(source_pixel_values, intermediate_pixel_values);

      BuildGraphFromImage(intermediate_pixel_values, target_graph);

      std::sort(target_graph.edges_.begin(), target_graph.edges_.end());

      DisjointSet vertex_disjoint_set(target_graph.vertices_.size());

      std::vector<double> thresholds(target_graph.vertices_.size(), 1 / k);

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
#pragma omp parallel for
      for (int r = 0; r < intermediate_pixel_values.size(); ++r) {
#pragma omp parallel for
        for (int c = 0; c < intermediate_pixel_values[r].size(); ++c) {
          size_t index = r * intermediate_pixel_values[r].size() + c;
          size_t group = vertex_disjoint_set.FindGroup(index);
          size_t group_size = vertex_disjoint_set.SizeOfGroup(group);
          if (group_size) {
            group_color[group].r_ += intermediate_pixel_values[r][c].r_ / (double)group_size;
            group_color[group].g_ += intermediate_pixel_values[r][c].g_ / (double)group_size;
            group_color[group].b_ += intermediate_pixel_values[r][c].b_ / (double)group_size;
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

        double color_difference = difference.NormalizedMagnitude();

        if (color_difference < similar_color_patch_merge_threshold) {
          vertex_disjoint_set.UnionGroup(group_of_x, group_of_y);
        }
      }

      // Calculate the color of each group again
      for (auto &color : group_color) {
        color = RGBA<double>();
      }

#pragma omp parallel for
      for (int r = 0; r < intermediate_pixel_values.size(); ++r) {
#pragma omp parallel for
        for (int c = 0; c < intermediate_pixel_values[r].size(); ++c) {
          size_t index = r * intermediate_pixel_values[r].size() + c;
          size_t group = vertex_disjoint_set.FindGroup(index);
          size_t group_size = vertex_disjoint_set.SizeOfGroup(group);
          if (group_size) {
            group_color[group].r_ += intermediate_pixel_values[r][c].r_ / (double)group_size;
            group_color[group].g_ += intermediate_pixel_values[r][c].g_ / (double)group_size;
            group_color[group].b_ += intermediate_pixel_values[r][c].b_ / (double)group_size;
          }
        }
      }

      // Write the pixel value
      group_of_pixel = std::vector<std::vector<size_t> >(source_pixel_values.size(), std::vector<size_t>(source_pixel_values[0].size()));

#pragma omp parallel for
      for (int r = 0; r < source_pixel_values.size(); ++r) {
#pragma omp parallel for
        for (int c = 0; c < source_pixel_values[r].size(); ++c) {
          size_t index = r * source_pixel_values[r].size() + c;
          size_t group = vertex_disjoint_set.FindGroup(index);
          group_of_pixel[r][c] = group;
          result_pixel_values[r][c] = RGBA<unsigned char>((unsigned char)group_color[group].r_, (unsigned char)group_color[group].g_, (unsigned char)group_color[group].b_);
        }
      }
    }

    void Segmentation(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const double k, const int min_patch_size, const double similar_color_patch_merge_threshold) {
      Graph<std::pair<size_t, size_t> > target_graph;
      std::vector<std::vector<size_t> > group_of_pixel;
      Segmentation(source_pixel_values, result_pixel_values, target_graph, group_of_pixel, k, min_patch_size, similar_color_patch_merge_threshold);
    }

    void Segmentation(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, std::vector<std::vector<size_t> > &group_of_pixel) {
      Graph<std::pair<size_t, size_t> > target_graph;
      const double SEGMENTATION_K = pow(source_pixel_values.size() * source_pixel_values[0].size(), 0.6);
      const double SEGMENTATION_MIN_PATCH_SIZE = (source_pixel_values.size() * source_pixel_values[0].size()) * 0.001;
      const double SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD = 20;
      Segmentation(source_pixel_values, result_pixel_values, target_graph, group_of_pixel, SEGMENTATION_K, SEGMENTATION_MIN_PATCH_SIZE, SEGMENTATION_SIMILAR_COLOR_MERGE_THRESHOLD);
    }

    void Segmentation(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      std::vector<std::vector<size_t> > group_of_pixel;
      Segmentation(source_pixel_values, result_pixel_values, group_of_pixel);
    }

    void ConfidenceMap(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, std::vector<std::vector<double> > &confidence_map, const RGBA<unsigned char> &target_color) {
      confidence_map = std::vector<std::vector<double> >(source_pixel_values.size(), std::vector<double>(source_pixel_values[0].size()));

      const double TARGET_COLOR_WEIGHT = 1.0;
      const double POSITION_WEIGHT = 0.0;
      const double SIMILIAR_WEIGHT = 0.0;

#pragma omp parallel for
      for (int r = 0; r < source_pixel_values.size(); ++r) {
#pragma omp parallel for
        for (int c = 0; c < source_pixel_values[r].size(); ++c) {
          double distance_with_target_color = 0;

          RGBA<double> difference_with_target_color;

          difference_with_target_color.r_ = target_color.r_ - (double)source_pixel_values[r][c].r_;
          difference_with_target_color.g_ = target_color.g_ - (double)source_pixel_values[r][c].g_;
          difference_with_target_color.b_ = target_color.b_ - (double)source_pixel_values[r][c].b_;

          distance_with_target_color = 1.0 - (difference_with_target_color.NormalizedMagnitude() / 255.0);

          //distance_with_target_color += pow(0.5 - (source_pixel_values[r][c].r_ / 255.0), 2.0);
          //distance_with_target_color += pow(1.0 - (source_pixel_values[r][c].g_ / 255.0), 2.0);
          //distance_with_target_color += pow(1.0 - (source_pixel_values[r][c].b_ / 255.0), 2.0);
          //distance_with_target_color = sqrt(distance_with_target_color / 3.0);

          double distance_with_center = sqrt(pow(r - source_pixel_values.size() / 2.0, 2.0) + pow(c - source_pixel_values[r].size() / 2.0, 2.0));
          distance_with_center /= sqrt(pow(source_pixel_values.size() / 2.0, 2.0) + pow(source_pixel_values[r].size() / 2.0, 2.0));
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
          //    if (IsInImage(source_pixel_values, neighbor_r, neighbor_c)) {
          //      ++neighbor_count;

          //      double distance_with_neighbor = pow((source_pixel_values[neighbor_r][neighbor_c].r_ / 255.0) - (source_pixel_values[r][c].r_ / 255.0), 2.0);
          //      distance_with_neighbor += pow((source_pixel_values[neighbor_r][neighbor_c].g_ / 255.0) - (source_pixel_values[r][c].g_ / 255.0), 2.0);
          //      distance_with_neighbor += pow((source_pixel_values[neighbor_r][neighbor_c].b_ / 255.0) - (source_pixel_values[r][c].b_ / 255.0), 2.0);
          //      distance_with_neighbor = sqrt(distance_with_neighbor / 3.0);

          //      total_distance_with_neighbor += distance_with_neighbor;
          //    }
          //  }
          //}

          //if (neighbor_count) {
          //  total_distance_with_neighbor /= neighbor_count;
          //}

          //total_distance_with_neighbor = 1 - total_distance_with_neighbor;

          double confidence_value = TARGET_COLOR_WEIGHT * distance_with_target_color + POSITION_WEIGHT * distance_with_center + SIMILIAR_WEIGHT * total_distance_with_neighbor;

          //confidence_value = pow(confidence_value, 4.5);

          confidence_map[r][c] = confidence_value;

          unsigned char confidence_gray_level = confidence_map[r][c] * 255;

          result_pixel_values[r][c].r_ = confidence_gray_level;
          result_pixel_values[r][c].g_ = confidence_gray_level;
          result_pixel_values[r][c].b_ = confidence_gray_level;
        }
      }
    }

    void ConfidenceMap(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const RGBA<unsigned char> &target_color) {
      std::vector<std::vector<double> > confidence_map;
      ConfidenceMap(source_pixel_values, result_pixel_values, confidence_map, target_color);
    }

    void SignificanceMap(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, std::vector<std::vector<double> > &significance_map, const RGBA<unsigned char> &target_color) {

      significance_map = std::vector<std::vector<double> >(source_pixel_values.size(), std::vector<double>(source_pixel_values[0].size()));

      std::vector<std::vector<size_t> > group_of_pixel;
      Segmentation(source_pixel_values, result_pixel_values, group_of_pixel);

      std::vector<std::vector<double> > confidence_map;
      ConfidenceMap(source_pixel_values, result_pixel_values, confidence_map, target_color);

      std::map<size_t, size_t> size_of_group;
      std::map<size_t, double> significance_of_group;

      for (size_t r = 0; r < source_pixel_values.size(); ++r) {
        for (size_t c = 0; c < source_pixel_values[r].size(); ++c) {
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

#pragma omp parallel for
      for (int r = 0; r < source_pixel_values.size(); ++r) {
#pragma omp parallel for
        for (int c = 0; c < source_pixel_values[r].size(); ++c) {
          size_t pixel_group = group_of_pixel[r][c];

          double group_significance = significance_of_group[pixel_group];
          significance_map[r][c] = group_significance;

          result_pixel_values[r][c] = SignificanceValueToSignifanceColor(group_significance);
        }
      }
    }

    void SignificanceMap(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values, const RGBA<unsigned char> &target_color) {
      std::vector<std::vector<double> > significance_map;
      SignificanceMap(source_pixel_values, result_pixel_values, significance_map, target_color);
    }

    void WoundSegmentation(const std::vector<std::vector<RGBA<unsigned char> > > &source_pixel_values, std::vector<std::vector<RGBA<unsigned char> > > &result_pixel_values) {
      std::vector<std::vector<double> > significance_map;
      SignificanceMap(source_pixel_values, result_pixel_values, significance_map, RGBA<unsigned char>(80, 0, 0));

#pragma omp parallel for
      for (int r = 0; r < source_pixel_values.size(); ++r) {
#pragma omp parallel for
        for (int c = 0; c < source_pixel_values[r].size(); ++c) {
          double pixel_significance = significance_map[r][c];
          RGBA<unsigned char> thresholding_color = (pixel_significance >= 0.8) ? RGBA<unsigned char>(0, 0, 0) : RGBA<unsigned char>(255, 255, 255);
          result_pixel_values[r][c] = thresholding_color;
        }
      }

      std::vector<std::vector<RGBA<unsigned char> > > intermediate_pixel_values = result_pixel_values;

      // Post processing

      for (size_t t = 0; t < 0; ++t) {
        MaxFilter(intermediate_pixel_values, result_pixel_values, 3);
        intermediate_pixel_values = result_pixel_values;
      }

      for (size_t t = 0; t < 0; ++t) {
        MinFilter(intermediate_pixel_values, result_pixel_values, 3);
        intermediate_pixel_values = result_pixel_values;
      }

      for (size_t t = 0; t < 1; ++t) {
        RemoveVerticalLines(intermediate_pixel_values, result_pixel_values, 40);
        intermediate_pixel_values = result_pixel_values;
        RemoveHorizontalLines(intermediate_pixel_values, result_pixel_values, 40);
        intermediate_pixel_values = result_pixel_values;
      }

      RemoveIsolatedRegion(intermediate_pixel_values, result_pixel_values);
      intermediate_pixel_values = result_pixel_values;

      for (size_t t = 0; t < 1; ++t) {
        FillHorizontalGaps(intermediate_pixel_values, result_pixel_values, 40);
        intermediate_pixel_values = result_pixel_values;
        FillVerticalGaps(intermediate_pixel_values, result_pixel_values, 40);
        intermediate_pixel_values = result_pixel_values;
      }

      for (size_t t = 0; t < 0; ++t) {
        RemoveVerticalLines(intermediate_pixel_values, result_pixel_values, 40);
        intermediate_pixel_values = result_pixel_values;
        RemoveHorizontalLines(intermediate_pixel_values, result_pixel_values, 40);
        intermediate_pixel_values = result_pixel_values;
        RemoveIsolatedRegion(intermediate_pixel_values, result_pixel_values);
        intermediate_pixel_values = result_pixel_values;
      }

      RGBA<double> average_wound_color(0, 0, 0);
      size_t wound_pixel_count = 0;
      for (int r = 0; r < result_pixel_values.size(); ++r) {
        for (int c = 0; c < result_pixel_values[r].size(); ++c) {
          if (!result_pixel_values[r][c].r_) {
            average_wound_color.r_ += source_pixel_values[r][c].r_;
            average_wound_color.g_ += source_pixel_values[r][c].g_;
            average_wound_color.b_ += source_pixel_values[r][c].b_;
            ++wound_pixel_count;
          }
        }
      }

      if (wound_pixel_count) {
        average_wound_color.r_ /= (double)wound_pixel_count;
        average_wound_color.g_ /= (double)wound_pixel_count;
        average_wound_color.b_ /= (double)wound_pixel_count;
      }
      
      //std::cout << average_wound_color.r_ << " " << average_wound_color.g_ << " " << average_wound_color.b_ << "\n";
      // Again
//      SignificanceMap(source_pixel_values, result_pixel_values, significance_map, RGBA<unsigned char>(average_wound_color.r_, average_wound_color.g_, average_wound_color.b_));
//
//#pragma omp parallel for
//      for (int r = 0; r < source_pixel_values.size(); ++r) {
//#pragma omp parallel for
//        for (int c = 0; c < source_pixel_values[r].size(); ++c) {
//          double pixel_significance = significance_map[r][c];
//          RGBA<unsigned char> thresholding_color = (pixel_significance >= 0.8) ? RGBA<unsigned char>(0, 0, 0) : RGBA<unsigned char>(255, 255, 255);
//          result_pixel_values[r][c] = thresholding_color;
//        }
//      }
//
//      intermediate_pixel_values = result_pixel_values;
//
//      for (size_t t = 0; t < 1; ++t) {
//        RemoveVerticalLines(intermediate_pixel_values, result_pixel_values, 40);
//        intermediate_pixel_values = result_pixel_values;
//        RemoveHorizontalLines(intermediate_pixel_values, result_pixel_values, 40);
//        intermediate_pixel_values = result_pixel_values;
//      }
//
//      RemoveIsolatedRegion(intermediate_pixel_values, result_pixel_values);
//      intermediate_pixel_values = result_pixel_values;
//
//      for (size_t t = 0; t < 1; ++t) {
//        FillHorizontalGaps(intermediate_pixel_values, result_pixel_values, 40);
//        intermediate_pixel_values = result_pixel_values;
//        FillVerticalGaps(intermediate_pixel_values, result_pixel_values, 40);
//        intermediate_pixel_values = result_pixel_values;
//      }

    }

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

      void Segmentation() {
        if (!image_loaded_) {
          return;
        }
        ImageProcessing::Segmentation(source_pixel_values_, result_pixel_values_);
      }

      void ConfidenceMap() {
        if (!image_loaded_) {
          return;
        }
        ImageProcessing::ConfidenceMap(source_pixel_values_, result_pixel_values_, RGBA<unsigned char>(127, 0, 0));
      }

      void SignificanceMap() {
        if (!image_loaded_) {
          return;
        }
        ImageProcessing::SignificanceMap(source_pixel_values_, result_pixel_values_, RGBA<unsigned char>(127, 0, 0));
      }

      void WoundSegmentation() {
        if (!image_loaded_) {
          return;
        }
        ImageProcessing::WoundSegmentation(source_pixel_values_, result_pixel_values_);
      }

      void WoundSegmentationWithOutlineOverlapping() {
        WoundSegmentation();

        std::vector<std::vector<RGBA<unsigned char> > > intermediate_pixel_values = source_pixel_values_;

        std::vector<std::vector<RGBA<unsigned char> > > result_outline_pixel_values = source_pixel_values_;

        ImageProcessing::SobelEdgeDetection(result_pixel_values_, intermediate_pixel_values);
        ImageProcessing::Thresholding(intermediate_pixel_values, result_outline_pixel_values, 127);

        result_pixel_values_ = source_pixel_values_;

        for (size_t r = 0; r < result_outline_pixel_values.size(); ++r) {
          for (size_t c = 0; c < result_outline_pixel_values[r].size(); ++c) {
            if (result_outline_pixel_values[r][c].r_) {
              result_pixel_values_[r][c] = RGBA<unsigned char>(255, 255, 0);
            }
          }
        }

        if (!ground_truth_image_loaded_) {
          return;
        }

        std::vector<std::vector<RGBA<unsigned char> > > ground_truth_outline_pixel_values = source_pixel_values_;

        ImageProcessing::SobelEdgeDetection(ground_truth_pixel_values_, intermediate_pixel_values);
        ImageProcessing::Thresholding(intermediate_pixel_values, ground_truth_outline_pixel_values, 127);

        for (size_t r = 0; r < result_outline_pixel_values.size(); ++r) {
          for (size_t c = 0; c < result_outline_pixel_values[r].size(); ++c) {
            if (ground_truth_outline_pixel_values[r][c].r_) {
              result_pixel_values_[r][c] = RGBA<unsigned char>(0, 255, 0);
            }
          }
        }
      }

      double DiceCoefficient() {
        if (!image_loaded_ || !ground_truth_image_loaded_) {
          return 0;
        }
        return ImageProcessing::DiceCoefficient(result_pixel_values_, ground_truth_pixel_values_);
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
            if (bytes_per_pixel == 1) {
              source_image_pixel_values[r][c].r_ = (255 - pixel_pointer[0]);
              source_image_pixel_values[r][c].g_ = (255 - pixel_pointer[0]);
              source_image_pixel_values[r][c].b_ = (255 - pixel_pointer[0]);
            } else {
              source_image_pixel_values[r][c].r_ = pixel_pointer[2];
              source_image_pixel_values[r][c].g_ = pixel_pointer[1];
              source_image_pixel_values[r][c].b_ = pixel_pointer[0];
            }
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

      bool image_loaded_;
      bool ground_truth_image_loaded_;

      std::vector<std::vector<RGBA<unsigned char> > > source_pixel_values_;
      std::vector<std::vector<RGBA<unsigned char> > > result_pixel_values_;
      std::vector<std::vector<RGBA<unsigned char> > > ground_truth_pixel_values_;
    };
  }
}