/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, October 17, 2017 - 10:24:56
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labraytracer/sphere.h>
#include <labraytracer/util.h>
#include <cmath>

namespace inviwo {

Sphere::Sphere(const vec3& center, const double& radius) {
    center_ = center;
    radius_ = radius;
}

bool Sphere::closestIntersection(const Ray& ray, double maxLambda,
                                 RayIntersection& intersection) const {
    //Programming TASK 1: implement this method
    //Your code should compute the intersection between a ray and a sphere;

    //If you detect an intersection, the return type should look similar to this:
    //if(rayIntersectsSphere)
    //{
    //  intersection = RayIntersection(ray,shared_from_this(),lambda,normalVec,uvw);
    //  return true;
    //}
    //
    //Hints:
    // lambda is the distance form the ray origin an the intersection point.
    // Ray origin p_r : ray.getOrigin()
    // Ray direction t_r : ray.getDirection()
    // If you need the intersection point, use ray.pointOnRay(lambda)
    // You can ignore the uvw (texture coordinates)

	vec3 p_r = ray.getOrigin();
	vec3 t_r = ray.getDirection();
	vec3 c = center_;
	double r = radius_;

	//We use the equation on slide 20 from lecture 8.
	double squared_coefficient = dot(t_r, t_r);
	double linear_coefficient = (2.0 * dot(t_r, (p_r - c))) / squared_coefficient;
	double constant_term = (dot((p_r - c), (p_r - c)) - r*r)  / squared_coefficient;

	if (pow(linear_coefficient / 2, 2) - constant_term < 0) {
		return false;
	}
	double lambda1 = -linear_coefficient / 2 - sqrt(pow(linear_coefficient / 2, 2) - constant_term);
	//double lambda2 = -linear_coefficient / 2 + sqrt(pow(linear_coefficient / 2, 2) - constant_term);
	if (lambda1 < 0) {
		return false;
	}
	if (lambda1 + Util::epsilon > maxLambda) {
		return false;
	}

	//The normal is given by the vector between the center of the sphere c
	//and the intercetion point p on the sphere.
	vec3 intersection_point = p_r + vec3(lambda1 * t_r.x, lambda1 * t_r.y, lambda1 * t_r.z);
	vec3 normalVec = intersection_point - c;
	normalVec = normalVec / sqrt(pow(normalVec.x, 2) + pow(normalVec.y, 2) + pow(normalVec.z, 2));

	const vec3 uvw(0, 0, 0);
	intersection = RayIntersection(ray, shared_from_this(), lambda1, normalVec, uvw);
	return true;
}

bool Sphere::anyIntersection(const Ray& ray, double maxLambda) const {
    RayIntersection temp;
    return closestIntersection(ray, maxLambda, temp);
}

void Sphere::drawGeometry(std::shared_ptr<BasicMesh> mesh,
                          std::vector<BasicMesh::Vertex>& vertices) const {
    auto indexBuffer = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    int lat = 8;
    int lon = 10;

    for (int i = 0; i < lat - 1; i++) {
        float theta1 = float(i * M_PI) / (lat - 1);
        float theta2 = float((i + 1) * M_PI) / (lat - 1);

        for (int j = 0; j < lon - 1; j++) {
            float phi1 = float(j * 2 * M_PI) / (lon - 1);
            float phi2 = float((j + 1) * 2 * M_PI) / (lon - 1);

            vec3 v1 = vec3(radius_ * sin(theta1) * cos(phi1), radius_ * sin(theta1) * sin(phi1),
                           radius_ * cos(theta1)) + center_;
            vec3 v2 = vec3(radius_ * sin(theta2) * cos(phi1), radius_ * sin(theta2) * sin(phi1),
                           radius_ * cos(theta2)) + center_;
            vec3 v3 = vec3(radius_ * sin(theta2) * cos(phi2), radius_ * sin(theta2) * sin(phi2),
                           radius_ * cos(theta2)) + center_;

            Util::drawLineSegment(v1, v2, vec4(0.2, 0.2, 0.2, 1), indexBuffer, vertices);
            Util::drawLineSegment(v2, v3, vec4(0.2, 0.2, 0.2, 1), indexBuffer, vertices);
        }
    }
}

} // namespace
