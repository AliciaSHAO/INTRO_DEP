#ifndef VEC2_H_
#define VEC2_H_


// ------------------------------------------------

typedef struct Vec2 {
  float x;
  float y;
} Vec2;

// ------------------------------------------------
double dot(Vec2 u, Vec2 v);
Vec2 add(Vec2 u, Vec2 v);
Vec2 opposite(Vec2 u);
Vec2 substract(Vec2 u, Vec2 v);
Vec2 multiply(Vec2 u, float x);
float norm(Vec2 u);
Vec2 normalize(Vec2 u);
void print_vec(Vec2 u);
Vec2 ortho_projection(Vec2 u, Vec2 v);
Vec2 orthogonal(Vec2 u);
Vec2 mirrorY(Vec2 u);
// ------------------------------------------------

#endif