/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, October 17, 2017 - 10:24:30
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labraytracer/triangle.h>
#include <labraytracer/util.h>
#include <memory>

namespace inviwo {

Triangle::Triangle() {
}

Triangle::Triangle(const vec3& v0, const vec3& v1, const vec3& v2, const vec3& uvw0,
                   const vec3& uvw1, const vec3& uvw2) {
    mVertices[0] = v0;
    mVertices[1] = v1;
    mVertices[2] = v2;
    mUVW[0] = uvw0;
    mUVW[1] = uvw1;
    mUVW[2] = uvw2;
}

bool Triangle::closestIntersection(const Ray& ray, double maxLambda,
                                   RayIntersection& intersection) const {
    //Programming TASK 1: Implement this method
    //Your code should compute the intersection between a ray and a triangle.
    //
    //If you detect an intersection, the return type should look similar to this:
    //if(rayIntersectsTriangle)
    //{
    //  intersection = RayIntersection(ray,shared_from_this(),lambda,ray.normal,uvw);
    //  return true;
    //} 
    //
    // Hints :
    // Ray origin p_r : ray.getOrigin()
    // Ray direction t_r : ray.getDirection()
    // Compute the intersection point using ray.pointOnRay(lambda)

	vec3 t_r = ray.getDirection();
	vec3 p_r = ray.getOrigin();
	vec3 p_0 = mVertices[0];
	vec3 t1 = mVertices[1] - p_0;
	vec3 t2 = mVertices[2] - p_0;
	vec3 normalVec = glm::cross(t1, t2);
	normalVec = glm::normalize(normalVec);

	double numerator = dot((p_0 - p_r), normalVec);
	double denominator = dot(t_r, normalVec);
	double lambda = numerator / denominator;

	vec3 p = ray.pointOnRay(lambda);
	//From p = p_0 + lambda1 * t_1 + lambda2 * t2, we get
	//lambda1 AND lambda2 ARE NOT RELATED TO lambda.
	double lambda2 = (t1.y * (p.x - p_0.x) - t1.x * (p.y - p_0.y)) / (t1.y * t2.x - t1.x * t2.y);
	double lambda1 = (p.y - p_0.y - lambda2 * t2.y) / (t1.y);

	if ((lambda1 < 0 || lambda2 < 0) || (lambda1 > 1 || lambda2 > 1) || ((lambda1 + lambda2) > 1)) {
		return false;
	}

	if (lambda + Util::epsilon > maxLambda) {
		return false;
	}

	const vec3 uvw(0, 0, 0);
	intersection = RayIntersection(ray, shared_from_this(), lambda, normalVec, uvw);
	return true;
}

bool Triangle::anyIntersection(const Ray& ray, double maxLambda) const {
    RayIntersection temp;
    return closestIntersection(ray, maxLambda, temp);
}

void Triangle::drawGeometry(std::shared_ptr<BasicMesh> mesh,
                            std::vector<BasicMesh::Vertex>& vertices) const {
    auto indexBuffer = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    Util::drawLineSegment(mVertices[0], mVertices[1], vec4(0.2, 0.2, 0.2, 1), indexBuffer,
                          vertices);
    Util::drawLineSegment(mVertices[1], mVertices[2], vec4(0.2, 0.2, 0.2, 1), indexBuffer,
                          vertices);
    Util::drawLineSegment(mVertices[2], mVertices[0], vec4(0.2, 0.2, 0.2, 1), indexBuffer,
                          vertices);
}

} // namespace
