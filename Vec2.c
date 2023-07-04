#include "Vec2.h"
#include <math.h>
#include <stdio.h>

// ------------------------------------------------

double dot(Vec2 u, Vec2 v) {
	double x = u.x * v.x + u.y * v.y;
	return x;
}

Vec2 add(Vec2 u, Vec2 v) {
	Vec2 w;
	w.x = u.x + v.x;
	w.y = u.y + v.y;
	return w;
}

Vec2 opposite(Vec2 u){
	Vec2 v;
	v.x = - u.x;
	v.y = - u.y;
	return v;
}

Vec2 substract(Vec2 u, Vec2 v) {
	Vec2 w;
	w.x = u.x - v.x;
	w.y = u.y - v.y;
	return w;
}

Vec2 multiply(Vec2 u, float x) {
	Vec2 v;
	v.x = x * u.x;
	v.y = x * u.y;
	return v;
}

float norm(Vec2 u){
	float n = sqrt(dot(u, u));
	return n;
}

Vec2 normalize(Vec2 u) {
	float n = sqrt(dot(u, u));
	Vec2 v;
	v.x = u.x / n;
	v.y = u.y / n;
	return v;
}

void print_vec(Vec2 u){
	printf("%f ê_x + %f ê_y\n", u.x, u.y);
}

Vec2 ortho_projection(Vec2 u, Vec2 v){
	//Returns the orthogonal projection of u on v.
return multiply(normalize(v), dot(u, v)/norm(v));
}

Vec2 orthogonal(Vec2 u){
	//Returns an unitary orthogonal vector of u.
	Vec2 v = {-1.0f, u.x/u.y};
	return normalize(v);
}

Vec2 mirrorY(Vec2 u){
	//Returns symetrical vector of u by Y axis.
	Vec2 v = {-u.x, u.y};
	return v;
}