#ifndef CONSTRAINT_H_
#define CONSTRAINT_H_

// ------------------------------------------------

typedef struct SphereCollider {
  Vec2 center;
  float radius;
} SphereCollider;

typedef struct PlaneCollider {
  // n is normal to the plane.
  // p is a point of the plane.
  Vec2 n;
  Vec2 p;
} PlaneCollider;

typedef struct AABBRectangleCollider {
  // AABB defined by upper left corner, and bottom right corner.
  Vec2 up_left;
  Vec2 bottom_right;
} AABBRectangleCollider;

typedef struct OBBRectangleCollider {
  // OBB defined by the rectangle's vertices.
  Vec2 u1;
  Vec2 u2;
  Vec2 u3;
  Vec2 u4;
} OBBRectangleCollider;
// ------------------------------------------------

#endif