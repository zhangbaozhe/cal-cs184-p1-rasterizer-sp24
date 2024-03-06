#include "rasterizer.h"

using namespace std;

namespace CGL {

RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                             size_t width, size_t height,
                             unsigned int sample_rate) {
  this->psm = psm;
  this->lsm = lsm;
  this->width = width;
  this->height = height;
  this->sample_rate = sample_rate;

  sample_buffer.resize(width * height * sample_rate, Color::White);
}

// Used by rasterize_point and rasterize_line
void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
  // TODO: Task 2: You might need to this function to fix points and lines (such
  // as the black rectangle border in test4.svg) NOTE: You are not required to
  // implement proper supersampling for points and lines It is sufficient to use
  // the same color for all supersamples of a pixel for points and lines (not
  // triangles)

  // sample_buffer[y * width + x] = c;
  auto start_iter = sample_buffer.begin() + (sample_rate * (y * width + x) + 0);
  auto end_iter = start_iter + sample_rate;
  std::fill(start_iter, end_iter, c);
}

// Rasterize a point: simple example to help you start familiarizing
// yourself with the starter code.
//
void RasterizerImp::rasterize_point(float x, float y, Color color) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width)
    return;
  if (sy < 0 || sy >= height)
    return;

  fill_pixel(sx, sy, color);
  return;
}

// Rasterize a line.
void RasterizerImp::rasterize_line(float x0, float y0, float x1, float y1,
                                   Color color) {
  if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  float pt[] = {x0, y0};
  float m = (y1 - y0) / (x1 - x0);
  float dpt[] = {1, m};
  int steep = abs(m) > 1;
  if (steep) {
    dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
    dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
  }

  while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
    rasterize_point(pt[0], pt[1], color);
    pt[0] += dpt[0];
    pt[1] += dpt[1];
  }
}

// Rasterize a triangle.
void RasterizerImp::rasterize_triangle(float x0, float y0, float x1, float y1,
                                       float x2, float y2, Color color) {
  // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
  // TODO: Task 2: Update to implement super-sampled rasterization

  Vector2D a(x0, y0);
  Vector2D b(x1, y1);
  Vector2D c(x2, y2);
  size_t supersample_width = sqrt(sample_rate);
  double offset = 1.0 / (supersample_width * 2.0);

  // winding a -> b -> c
  auto z = cross(b - a, c - a);
  if (z < 0) {
    swap(b, c);
  }

  const double denominator = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
  const auto barycentric_coord = [denominator, a, b, c] (const Vector2D &p) -> Vector3D {
    Vector3D result;
    result.x = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / denominator;
    result.y = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / denominator;
    result.z = 1 - result.x - result.y;
    return result;
  };

  auto [x_min, x_max] = minmax({x0, x1, x2});
  auto [y_min, y_max] = minmax({y0, y1, y2});

  size_t px_min = x_min < 0.0F ? 0 : (size_t)x_min;
  size_t px_max = x_max > (float)width ? width - 1 : (size_t)x_max;
  size_t py_min = y_min < 0.0F ? 0 : (size_t)y_min;
  size_t py_max = y_max > (float)height ? height - 1 : (size_t)y_max;

  // TODO(Baozhe Zhang): make this faster 
  // now is slow
  for (size_t px = px_min; px <= px_max; px++) {
    for (size_t py = py_min; py <= py_max; py++) {

      // trivial super sampling iterations
      for (size_t x_idx_super = 0; x_idx_super < supersample_width; x_idx_super++) {
        for (size_t y_idx_super = 0; y_idx_super < supersample_width; y_idx_super++) {

          // sample point
          Vector2D p(px + (double)x_idx_super / supersample_width + offset,
                     py + (double)y_idx_super / supersample_width + offset);

          Vector3D lambda = barycentric_coord(p);

          // (super) sample idx
          size_t idx = sample_rate * (py * width + px) + y_idx_super * supersample_width + x_idx_super;

          // auto z1 = cross(b - a, p - a);
          // auto z2 = cross(c - b, p - b);
          // auto z3 = cross(a - c, p - c);
          if (lambda.x >= 0.0 
              && lambda.x <= 1.0
              && lambda.y >= 0.0
              && lambda.y <= 1.0
              && lambda.z >= 0.0
              && lambda.z <= 1.0) {
            // rasterize_point((float)p.x, (float)p.y, color);
            sample_buffer.at(idx) = color;
          }
        }
      }
    }
  }
}

void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0,
                                                          Color c0, float x1,
                                                          float y1, Color c1,
                                                          float x2, float y2,
                                                          Color c2) {
  // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates
  // and using them to interpolate vertex colors across the triangle Hint: You
  // can reuse code from rasterize_triangle]

  Vector2D a(x0, y0);
  Vector2D b(x1, y1);
  Vector2D c(x2, y2);
  size_t supersample_width = sqrt(sample_rate);
  double offset = 1.0 / (supersample_width * 2.0);

  // winding a -> b -> c
  auto z = cross(b - a, c - a);
  if (z < 0) {
    swap(b, c);
  }

  const double denominator = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
  const auto barycentric_coord = [denominator, a, b, c] (const Vector2D &p) -> Vector3D {
    Vector3D result;
    result.x = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / denominator;
    result.y = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / denominator;
    result.z = 1 - result.x - result.y;
    return result;
  };

  auto [x_min, x_max] = minmax({x0, x1, x2});
  auto [y_min, y_max] = minmax({y0, y1, y2});

  size_t px_min = x_min < 0.0F ? 0 : (size_t)x_min;
  size_t px_max = x_max > (float)width ? width - 1 : (size_t)x_max;
  size_t py_min = y_min < 0.0F ? 0 : (size_t)y_min;
  size_t py_max = y_max > (float)height ? height - 1 : (size_t)y_max;

  // TODO(Baozhe Zhang): make this faster 
  // now is slow
  for (size_t px = px_min; px <= px_max; px++) {
    for (size_t py = py_min; py <= py_max; py++) {

      // trivial super sampling iterations
      for (size_t x_idx_super = 0; x_idx_super < supersample_width; x_idx_super++) {
        for (size_t y_idx_super = 0; y_idx_super < supersample_width; y_idx_super++) {

          // sample point
          Vector2D p(px + (double)x_idx_super / supersample_width + offset,
                     py + (double)y_idx_super / supersample_width + offset);

          Vector3D lambda = barycentric_coord(p);
          Color color = lambda.x * c0 + lambda.y * c1 + lambda.z * c2;

          // (super) sample idx
          size_t idx = sample_rate * (py * width + px) + y_idx_super * supersample_width + x_idx_super;

          // auto z1 = cross(b - a, p - a);
          // auto z2 = cross(c - b, p - b);
          // auto z3 = cross(a - c, p - c);
          if (lambda.x >= 0.0 
              && lambda.x <= 1.0
              && lambda.y >= 0.0
              && lambda.y <= 1.0
              && lambda.z >= 0.0
              && lambda.z <= 1.0) {
            // rasterize_point((float)p.x, (float)p.y, color);
            sample_buffer.at(idx) = color;
          }
        }
      }
    }
  }
}

void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0,
                                                float v0, float x1, float y1,
                                                float u1, float v1, float x2,
                                                float y2, float u2, float v2,
                                                Texture &tex) {
  // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample
  // function.
  // TODO: Task 6: Set the correct barycentric differentials in the SampleParams
  // struct. Hint: You can reuse code from
  // rasterize_triangle/rasterize_interpolated_color_triangle


  Vector2D a(x0, y0);
  Vector2D b(x1, y1);
  Vector2D c(x2, y2);
  Vector2D uv0(u0, v0);
  Vector2D uv1(u1, v1);
  Vector2D uv2(u2, v2);
  size_t supersample_width = sqrt(sample_rate);
  double offset = 1.0 / (supersample_width * 2.0);

  // winding a -> b -> c
  auto z = cross(b - a, c - a);
  if (z < 0) {
    swap(b, c);
  }

  const double denominator = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
  const auto barycentric_coord = [denominator, a, b, c] (const Vector2D &p) -> Vector3D {
    Vector3D result;
    result.x = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / denominator;
    result.y = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / denominator;
    result.z = 1 - result.x - result.y;
    return result;
  };

  auto [x_min, x_max] = minmax({x0, x1, x2});
  auto [y_min, y_max] = minmax({y0, y1, y2});

  size_t px_min = x_min < 0.0F ? 0 : (size_t)x_min;
  size_t px_max = x_max > (float)width ? width - 1 : (size_t)x_max;
  size_t py_min = y_min < 0.0F ? 0 : (size_t)y_min;
  size_t py_max = y_max > (float)height ? height - 1 : (size_t)y_max;

  // TODO(Baozhe Zhang): make this faster 
  // now is slow
  for (size_t px = px_min; px <= px_max; px++) {
    for (size_t py = py_min; py <= py_max; py++) {

      // trivial super sampling iterations
      for (size_t x_idx_super = 0; x_idx_super < supersample_width; x_idx_super++) {
        for (size_t y_idx_super = 0; y_idx_super < supersample_width; y_idx_super++) {

          // sample point
          Vector2D p(px + (double)x_idx_super / supersample_width + offset,
                     py + (double)y_idx_super / supersample_width + offset);

          Vector3D lambda = barycentric_coord(p);
          Vector2D uv = lambda.x * uv0 + lambda.y * uv1 + lambda.z * uv2;

          // (super) sample idx
          size_t idx = sample_rate * (py * width + px) + y_idx_super * supersample_width + x_idx_super;

          // auto z1 = cross(b - a, p - a);
          // auto z2 = cross(c - b, p - b);
          // auto z3 = cross(a - c, p - c);
          if (lambda.x >= 0.0 
              && lambda.x <= 1.0
              && lambda.y >= 0.0
              && lambda.y <= 1.0
              && lambda.z >= 0.0
              && lambda.z <= 1.0) {
            // rasterize_point((float)p.x, (float)p.y, color);
            Color color = tex.sample({
                uv, Vector2D(), Vector2D(), this->psm, this->lsm
            });
            sample_buffer[idx] = color;
          }
        }
      }
    }
  }
}

void RasterizerImp::set_sample_rate(unsigned int rate) {
  // TODO: Task 2: You may want to update this function for supersampling
  // support

  this->sample_rate = rate;

  this->sample_buffer.resize(width * height * sample_rate, Color::White);
}

void RasterizerImp::set_framebuffer_target(unsigned char *rgb_framebuffer,
                                           size_t width, size_t height) {
  // TODO: Task 2: You may want to update this function for supersampling
  // support

  this->width = width;
  this->height = height;
  this->rgb_framebuffer_target = rgb_framebuffer;

  this->sample_buffer.resize(width * height * sample_rate, Color::White);
}

void RasterizerImp::clear_buffers() {
  std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height,
            255);
  std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
}

// This function is called at the end of rasterizing all elements of the
// SVG file.  If you use a supersample buffer to rasterize SVG elements
// for antialising, you could use this call to fill the target framebuffer
// pixels from the supersample buffer data.
//
void RasterizerImp::resolve_to_framebuffer() {
  // TODO: Task 2: You will likely want to update this function for
  // supersampling support

  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      // Avg down sampling
      Color col(0, 0, 0);
      for (size_t i = 0; i < sample_rate; i++) {
        col += sample_buffer.at(sample_rate * (y * width + x) + i);
      }
      col = std::move(
          Color(col.r / sample_rate, col.g / sample_rate, col.b / sample_rate));

      for (int k = 0; k < 3; ++k) {
        // FIXME: floating point error?
        // When I test interpolation (test7), there are 4 pixels that are not identical to 
        // the reference one. The R, G, or B channel's values are larger by one.
        this->rgb_framebuffer_target[3 * (y * width + x) + k] =
            static_cast<unsigned char>((&col.r)[k] * 255.0F); 
      }
    }
  }
}

Rasterizer::~Rasterizer() {}

} // namespace CGL
