#include "Context.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))


//int main() {
//    return 0;
//}

// ------------------------------------------------

Particle getParticle(Context* context, int id)
{
  return context->particles[id];
}

// ------------------------------------------------

void addParticle(Context* context, float x, float y, float radius, float mass, int draw_id)
{
    assert(context->num_particles<context->capacity_particles); // currently no resize in context
    context->particles[context->num_particles].position.x = x;
    context->particles[context->num_particles].position.y = y;
    context->particles[context->num_particles].velocity.x = 0.0F;
    context->particles[context->num_particles].velocity.y = 0.0F;
    context->particles[context->num_particles].inv_mass = 1.0F/mass;
    context->particles[context->num_particles].radius = radius;
    context->particles[context->num_particles].draw_id = draw_id;
    context->num_particles += 1;
}

// ------------------------------------------------

void setDrawId(Context* context, int sphere_id, int draw_id)
{
  context->particles[sphere_id].draw_id = draw_id;
}

// ------------------------------------------------

SphereCollider getGroundSphereCollider(Context* context, int id)
{
  return context->ground_spheres[id];
}

AABBRectangleCollider getAABBRectangleCollider(Context* context, int id){
  return context->AABBrectangles[id];
}

OBBRectangleCollider getOBBRectangleCollider(Context* context, int id){
  return context->OBBrectangles[id];
}


// ------------------------------------------------

int countPins(int nb_rows, int nb_steps){
  //Returns the number of pins in a checkerboard pattern of pins with nb_rows rows and nb_steps steps,
  //where uneven rows have a pin at x = 0, and even rows do not.
  //There is a step of 1.0 along the X axis between each pin, nb_steps is the number of steps in the x < 0  (or x > 0) space.
  //There is a space of 0.5 along the Y axis between each row.
  int n;
  //There are 2*nb_steps + 1 pins for an uneven row.
  //2* nb_steps for an even row.
  if (nb_rows%2){
    //If there is a uneven number of rows there are (nb_rows + 1)/2 uneven rows and (nb_rows - 1)/2 even rows.
    n = (nb_rows + 1)/2 * ((nb_steps*2) + 1) + (nb_rows - 1) * nb_steps;
  } 
  else{
    //If there is a uneven number of rows there are nb_rows/2 uneven rows and even rows.
    n = (nb_rows/2) * ((nb_steps*2) + 1) + nb_rows*nb_steps;
  }
  return n;
}

void GaltonPins(Context* context, int rows, int steps){
  //Forms the board's pins by adding or substracting the unitary vector {1, 0} to;
  //the centered pin for uneven rows,
  //the symetrical pins at (-0.5, y_n), (+0.5, y_n) for the n_th even row.
  int counter = 0;

  for (int i = 0; i < rows; ++i){
    Vec2 u = {0.0f, 4.75 - 0.5*i};
    Vec2 e_x = {1.0f, 0.0f};
    if (1 - i%2){
      context->ground_spheres[counter].center = u;
      context->ground_spheres[counter].radius = 0.08;
      ++counter;
      for (int k = 0; k < steps; ++k){
        context->ground_spheres[counter].center = add(u, multiply(e_x, k + 1));
        context->ground_spheres[counter].radius = 0.08;
        ++counter;
        context->ground_spheres[counter].center = add(u, multiply(e_x, -(k + 1)));
        context->ground_spheres[counter].radius = 0.08;
        ++counter;
      }
    }
    else {
      Vec2 half_e_x = {0.5f, 0.0f};
      Vec2 u_p = add(u, half_e_x);
      Vec2 u_m = substract(u, half_e_x);
      context->ground_spheres[counter].center = u_p;
      context->ground_spheres[counter].radius = 0.08;
      ++counter;
      context->ground_spheres[counter].center = u_m;
      context->ground_spheres[counter].radius = 0.08;
      ++counter;
      for (int k = 0; k < steps - 1; ++k){
        context->ground_spheres[counter].center = add(u_p, multiply(e_x, (k + 1)));
        context->ground_spheres[counter].radius = 0.08;
        ++counter;
        context->ground_spheres[counter].center = add(u_m, multiply(e_x, -(k + 1)));
        context->ground_spheres[counter].radius = 0.08;
        ++counter;
      }
    }
  }
}


Context* initializeContext(int capacity)
{
  Context* context = malloc(sizeof(Context));

  //Particles
  context->num_particles = 0;
  context->capacity_particles = capacity;
  context->particles = malloc(capacity*sizeof(Particle));
  memset(context->particles,0,capacity*sizeof(Particle));

  //Galton Pins = Sphere Colliders
  int rows = 19;
  int steps = 10;
 
  context->num_ground_sphere = countPins(rows, steps);
  context->ground_spheres = malloc(context->num_ground_sphere*sizeof(SphereCollider));
  
  GaltonPins(context, rows, steps);

  //Planes
  context -> num_planes = 4;
  context -> planes = malloc(context->num_planes* sizeof(PlaneCollider));

  Vec2 n1 = {0.0f, 1.0f};
  Vec2 p1 = {0.0f, -10.0f};
  context -> planes[0].n = n1;
  context -> planes[0].p = p1;

  Vec2 n2 = {-1.0f, 0.0f};
  Vec2 p2 = {10.0f, 0.0f};
  context -> planes[1].n = n2;
  context -> planes[1].p = p2;

  Vec2 n3 = {0.0f, -1.0f};
  Vec2 p3 = {0.0f, 10.0f};
  context -> planes[2].n = n3;
  context -> planes[2].p = p3;

  Vec2 n4 = {1.0f, 0.0f};
  Vec2 p4 = {-10.0f, 0.0f};
  context -> planes[3].n = n4;
  context -> planes[3].p = p4;

  //AABB: Galton boxes
  context -> num_AABBrectangles = 19;
  context -> AABBrectangles = malloc(context->num_AABBrectangles*sizeof(AABBRectangleCollider));
  for (int i = 0; i < context->num_AABBrectangles; ++i){
    //Walls are 5 tall, 0.1 thick and every 1.0.
    Vec2 u = {(float)((i+1) * 1) - 10.0f, -5.0f};
    Vec2 v = {(float)((i+1) * 1) - 10.0f + 0.1, -10.0f};
    context -> AABBrectangles[i].up_left = u;
    context -> AABBrectangles[i].bottom_right = v;
  }
 
  //OBB: funnel
  context -> num_OBBrectangles = 2;
  context -> OBBrectangles = malloc(context->num_OBBrectangles*sizeof(OBBRectangleCollider));
  Vec2 top_center = {0.0f, 5.0f};
  Vec2 top_left = {8.0f, 10.0f};

  //The Right OBB follows u's direction.
  Vec2 u = substract(top_left, top_center);
  context -> OBBrectangles[0].u1 = add(top_center, add(multiply(normalize(u), 0.5f), multiply(orthogonal(u), 0.1)));
  context -> OBBrectangles[0].u2 = add(top_center, add(multiply(normalize(u), 0.5f), multiply(orthogonal(u), -0.1)));
  context -> OBBrectangles[0].u3 = add(top_center, add(multiply(normalize(u), norm(u) - 1.0f), multiply(orthogonal(u), -0.1)));
  context -> OBBrectangles[0].u4 = add(top_center, add(multiply(normalize(u), norm(u) - 1.0f), multiply(orthogonal(u), 0.1)));

  //The left OBB is the right's mirror.
  context -> OBBrectangles[1].u1 = mirrorY(add(top_center, add(multiply(normalize(u), 0.5f), multiply(orthogonal(u), 0.1))));
  context -> OBBrectangles[1].u2 = mirrorY(add(top_center, add(multiply(normalize(u), 0.5f), multiply(orthogonal(u), -0.1))));
  context -> OBBrectangles[1].u3 = mirrorY(add(top_center, add(multiply(normalize(u), norm(u) - 1.0f), multiply(orthogonal(u), -0.1))));
  context -> OBBrectangles[1].u4 = mirrorY(add(top_center, add(multiply(normalize(u), norm(u) - 1.0f), multiply(orthogonal(u), 0.1))));

  return context;
}

// ------------------------------------------------

void updatePhysicalSystem(Context* context, float dt)
{
  applyExternalForce(context, dt);
  dampVelocities(context);
  updateExpectedPosition(context, dt);
  addDynamicContactConstraints(context);
  addStaticContactConstraints(context);

  updateVelocityAndPosition(context, dt);
  applyFriction(context);
}

// ------------------------------------------------
// Position Based Dynamics.
void applyExternalForce(Context* context, float dt)
{
    if (context->num_particles == 0) return;
    for (int i = 0; i < context->num_particles; ++i) {
        context->particles[i].velocity.y -= dt * 9.81;
    }
}

void dampVelocities(Context* context)
{
  
}

void updateExpectedPosition(Context* context, float dt)
{
    if (context->num_particles == 0) return;
    for (int i = 0; i < context->num_particles; ++i) {
        context->particles[i].next_pos = add(multiply(context->particles[i].velocity, dt), context->particles[i].position);
    }
}


// Dynamic constraints.
void checkContactWithParticle(Context* context, int particle_id){
  for (int i = 0; i < context -> num_particles; ++i){
    if (particle_id != i && norm(substract(context -> particles[i].next_pos, context -> particles[particle_id].next_pos)) < 1){
      Vec2 x_ji = substract(context -> particles[i].next_pos, context -> particles[particle_id].next_pos);
      float C = norm(x_ji) - (context -> particles[i].radius + context -> particles[particle_id].radius);
      if (C < 0){
        float sigma_i = ((context -> particles[i].inv_mass) / (context -> particles[i].inv_mass + context -> particles[particle_id].inv_mass)) * C;
        Vec2 delta_i = multiply(normalize(x_ji), -sigma_i);
        float sigma_j = ((context -> particles[particle_id].inv_mass) / (context -> particles[i].inv_mass + context -> particles[particle_id].inv_mass)) * C;
        Vec2 delta_j = multiply(normalize(x_ji), sigma_j);
        context -> particles[i].next_pos = add(context -> particles[i].next_pos, delta_i);
      }
    }
  }
}

void addDynamicContactConstraints(Context* context)
{
  if (context->num_particles == 0) return;
  for (int i = 0; i < context->num_particles; ++i) {
    checkContactWithParticle(context, i);
  }
}

//Static constraints
//Plane collisions (equations provided in the instructions)
void checkContactWithPlane(Context* context, int particle_id)
{
  for (int i = 0; i < context->num_planes; ++i){
    float dotp = dot(substract(context->particles[particle_id].next_pos, context->planes[i].p), context->planes[i].n);
    if (dotp - context->particles[particle_id].radius < 0){
      Vec2 qc = substract(context->particles[particle_id].next_pos, multiply(context->planes[i].n, dotp));
      float C = dot(substract(context->particles[particle_id].next_pos, qc),context->planes[i].n) - context->particles[particle_id].radius;
      Vec2 delta = opposite(multiply(context->planes[i].n, C));
      context->particles[particle_id].next_pos = add(context->particles[particle_id].next_pos, delta);
    }
  }
}

//AABB and OBB Collisions with spherical object.

//The following two functions return center point of a rectangle.
Vec2 getAABBRectangleCenter(AABBRectangleCollider rectangle){
  return multiply(add(rectangle.up_left, rectangle.bottom_right), 0.5);
}

Vec2 getOBBRectangleCenter(OBBRectangleCollider rectangle){
  return multiply(add(rectangle.u1, rectangle.u3), 0.5);
}

//The following two functions return the width of a rectangle.
float getAABBRectangleWidth(AABBRectangleCollider rectangle){
  Vec2 u = {rectangle.bottom_right.x, rectangle.up_left.y};
  return norm(substract(u, rectangle.up_left));
}

float getOBBRectangleWidth(OBBRectangleCollider rectangle){
  return norm(substract(rectangle.u1, rectangle.u2));
}

//The following two functions return the height of a rectangle.
float getAABBRectangleHeight(AABBRectangleCollider rectangle){
  Vec2 u = {rectangle.bottom_right.x, rectangle.up_left.y};
  return norm(substract(u, rectangle.bottom_right));
}

float getOBBRectangleHeight(OBBRectangleCollider rectangle){
  return norm(substract(rectangle.u2, rectangle.u3));
}


float clamp(float x, float lower_bound, float upper_bound){
  //Clamps a value x in the [lower_bound, upper_bound] interval
  return max(lower_bound, min(upper_bound, x));
}

//The following two functions return the closest point to center in the rectangle.
Vec2 ClosestPointAABBRectangle(AABBRectangleCollider rectangle, Vec2 center){
  
  Vec2 B = getAABBRectangleCenter(rectangle);
  //D gives the direction of the line formed by center and the rectangle's center.
  //B is considered to be (0, 0), we need to clamp D_x to [-width/2, width/2] and D_y to [-height/2, height/2]
  Vec2 D = substract(center, B);
  float w = getAABBRectangleWidth(rectangle);
  float h = getAABBRectangleHeight(rectangle);
  float p_x = clamp(D.x, -0.5f*w, 0.5f*w);
  float p_y = clamp(D.y, -0.5f*h, 0.5f*h);
  Vec2 p = {p_x, p_y};
  //Let's not forget to add B to get p's real coordinates in the plane.
  return add(B,p);
}

Vec2 ClosestPointOBBRectangle(OBBRectangleCollider rectangle, Vec2 center){
  Vec2 B = getOBBRectangleCenter(rectangle);
  Vec2 D = substract(center, B);
  float w = getOBBRectangleWidth(rectangle);
  float h = getOBBRectangleHeight(rectangle);
  //Similar process, except OBB requires changing vector bases.
  Vec2 e_x = normalize(substract(rectangle.u2, rectangle.u1));
  Vec2 e_y = normalize(substract(rectangle.u2, rectangle.u3));
  float p_x = clamp(dot(D, e_x), -0.5f*w, 0.5f*w);
  float p_y = clamp(dot(D, e_y), -0.5f*h, 0.5f*h);
  Vec2 p = add(multiply(e_x, p_x), multiply(e_y, p_y));
  return add(B,p);
}

void checkContactWithAABBRectangle(Context* context, int particle_id){
  //Equations found in the instructions still apply.
  //pc is found through previous functions.
  //nc is calculated with pc and next_pos.
  for (int i = 0; i < context -> num_AABBrectangles; ++i){
    if (context -> particles[particle_id].next_pos.y < -4.0){
      Vec2 p = ClosestPointAABBRectangle(context->AABBrectangles[i], context->particles[particle_id].next_pos);
      float sdf = norm(substract(context->particles[particle_id].next_pos, p)) - context-> particles[particle_id].radius;
      if (sdf < 0){
        Vec2 nc = normalize(substract(context->particles[particle_id].next_pos, p));
        float dotp = dot(substract(context->particles[particle_id].next_pos, p), nc);
        Vec2 qc = substract(context->particles[particle_id].next_pos, multiply(nc, dotp));
        float C = dot(substract(context->particles[particle_id].next_pos, qc), nc) - context->particles[particle_id].radius;
        Vec2 delta = multiply(nc,-C);
        context->particles[particle_id].next_pos = add(context->particles[particle_id].next_pos, delta);
      }
    }
  }
}

void checkContactWithOBBRectangle(Context* context, int particle_id){
  for (int i = 0; i < context -> num_OBBrectangles; ++i){
    Vec2 p = ClosestPointOBBRectangle(context->OBBrectangles[i], context->particles[particle_id].next_pos);
    float sdf = norm(substract(context->particles[particle_id].next_pos, p)) - context-> particles[particle_id].radius;
    if (sdf < 0){
      Vec2 nc = normalize(substract(context->particles[particle_id].next_pos, p));
      float dotp = dot(substract(context->particles[particle_id].next_pos, p), nc);
      Vec2 qc = substract(context->particles[particle_id].next_pos, multiply(nc, dotp));
      float C = dot(substract(context->particles[particle_id].next_pos, qc), nc) - context->particles[particle_id].radius;
      Vec2 delta = multiply(nc,-C);
      context->particles[particle_id].next_pos = add(context->particles[particle_id].next_pos, delta);
    }
  }
}

void checkContactWithSphere(Context* context, int particle_id)
//Used equations in instructions.
{
  for (int i = 0; i < context -> num_ground_sphere; ++i){
    if (norm(substract(context->particles[particle_id].next_pos, context->ground_spheres[i].center)) < 1.5){
      float sdf = norm(substract(context -> particles[particle_id].next_pos, context -> ground_spheres[i].center)) - (context -> ground_spheres[i].radius + context -> particles[particle_id].radius);
      if (sdf < 0){
        Vec2 nc = normalize(substract(context -> particles[particle_id].next_pos, context -> ground_spheres[i].center));
        Vec2 pc = substract(context -> particles[particle_id].next_pos, multiply(nc, sdf));

        float dotp = dot(substract(context->particles[particle_id].next_pos, pc), nc);
        Vec2 qc = substract(context->particles[particle_id].next_pos, multiply(nc, dotp));
        float C = dot(substract(context->particles[particle_id].next_pos, qc), nc); //- context->particles[particle_id].radius;
        Vec2 delta = multiply(nc,-C);
        context->particles[particle_id].next_pos = add(context->particles[particle_id].next_pos, delta);
      }
    }
  }
}


void addStaticContactConstraints(Context* context)
{
    if (context->num_particles == 0) return;
      for (int i = 0; i < context->num_particles; ++i) {
        checkContactWithPlane(context, i);
        checkContactWithSphere(context, i);
        checkContactWithAABBRectangle(context, i);
        checkContactWithOBBRectangle(context, i);
      }
}

void updateVelocityAndPosition(Context* context, float dt)
{
    if (context->num_particles == 0) return;
    for (int i = 0; i < context->num_particles; ++i) {
        context->particles[i].velocity = multiply(substract(context->particles[i].next_pos, context->particles[i].position), 1/dt);
        context->particles[i].position = context->particles[i].next_pos;
    }
}

void applyFriction(Context* context)
{
}


// ------------------------------------------------
