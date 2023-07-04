#pragma once

#include "Vec2.h"
#include "Particle.h"
#include "Constraint.h"



// ------------------------------------------------

typedef struct Context {
  int num_particles;
  int capacity_particles;
  Particle* particles;

  // Ground colliders 
  int num_ground_sphere;
  SphereCollider* ground_spheres;

  int num_planes;
  PlaneCollider* planes;

  // AABB represents Galton board "boxes".
  int num_AABBrectangles;
  AABBRectangleCollider* AABBrectangles;

  // OBB are used to make a funnel for the particles at the top of the board.
  int num_OBBrectangles;
  OBBRectangleCollider* OBBrectangles;

} Context;

// ------------------------------------------------

Context* initializeContext(int capacity);

// ------------------------------------------------

void addParticle(Context* context, float x, float y, float radius, float mass, int draw_id);

// ------------------------------------------------

Particle getParticle(Context* context, int id);

SphereCollider getGroundSphereCollider(Context* context, int id);

// ------------------------------------------------

void setDrawId(Context* context, int sphere_id, int draw_id);

// ------------------------------------------------
void updatePhysicalSystem(Context* context, float dt);

// ------------------------------------------------
// Methods below are called by updatePhysicalSystem
// ------------------------------------------------

void applyExternalForce(Context* context, float dt);
void dampVelocities(Context* context);
void updateExpectedPosition(Context* context, float dt);
void addDynamicContactConstraints(Context* context);
void addStaticContactConstraints(Context* context);
void updateVelocityAndPosition(Context* context, float dt);
void applyFriction(Context* context);

void checkContactWithPlane(Context* context, int particle_id);
void checkContactWithSphere(Context* context, int particle_id);
void checkContactWithParticle(Context* context, int particle_id);

AABBRectangleCollider getAABBRectangleCollider(Context* context, int id);
Vec2 getAABBRectangleCenter(AABBRectangleCollider rectangle);
float getAABBRectangleWidth(AABBRectangleCollider rectangle);
float getAABBRectangleHeight(AABBRectangleCollider rectangle);
float clamp(float x, float lower_bound, float upper_bound);
Vec2 ClosestPointAABBRectangle(AABBRectangleCollider rectangle, Vec2 center);
void checkContactWithAABBRectangle(Context* context, int particle_id);

OBBRectangleCollider getOBBRectangleCollider(Context* context, int id);
Vec2 getOBBRectangleCenter(OBBRectangleCollider rectangle);
float getOBBRectangleWidth(OBBRectangleCollider rectangle);
float getOBBRectangleHeight(OBBRectangleCollider rectangle);
Vec2 ClosestPointOBBRectangle(OBBRectangleCollider rectangle, Vec2 center);
void checkContactWithOBBRectangle(Context* context, int particle_id);

int countPins(int rows, int steps);
void GaltonPins(Context* context, int rows, int steps);


// ------------------------------------------------

