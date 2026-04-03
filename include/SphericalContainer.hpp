/**
 * @file SphericalContainer.hpp
 * @brief Spherical boundary container for confining chromatin beads.
 *
 * Provides containment checks and energy penalties to keep the
 * reconstructed 3D structure within a spherical volume (e.g. a nucleus).
 * The boundary is optional: when radius <= 0 it has no effect.
 */

#ifndef SPHERICALCONTAINER_H_
#define SPHERICALCONTAINER_H_

#include <mtxlib.h>

/**
 * @class SphericalContainer
 * @brief Spherical boundary that penalizes beads outside a given radius.
 *
 * Usage:
 * @code
 *   SphericalContainer sphere(vector3(0,0,0), 10.0f);
 *   if (!sphere.contains(pos))
 *     energy += sphere.boundaryPenalty(pos);
 * @endcode
 */
class SphericalContainer {
public:
  /**
   * @brief Construct a spherical container.
   * @param center Center of the sphere in 3D space.
   * @param radius Radius of the sphere (must be > 0 for meaningful use).
   */
  SphericalContainer(const vector3 &center, float radius);

  /**
   * @brief Check whether a point lies inside (or on) the sphere.
   * @param pos 3D position to test.
   * @return True if the point is inside or on the sphere surface.
   */
  bool contains(const vector3 &pos) const;

  /**
   * @brief Project a point onto the sphere surface if it lies outside.
   *
   * If the point is inside the sphere, it is returned unchanged.
   * If the point is at the center, a point on the +x surface is returned.
   *
   * @param pos 3D position to clamp.
   * @return The clamped position (on or inside the sphere).
   */
  vector3 clampToSphere(const vector3 &pos) const;

  /**
   * @brief Compute a quadratic energy penalty for being outside the sphere.
   *
   * Returns 0 if the point is inside or on the sphere.
   * Otherwise returns k * (distance_from_surface)^2, where k = 1.0.
   *
   * @param pos 3D position to evaluate.
   * @return Energy penalty (>= 0).
   */
  float boundaryPenalty(const vector3 &pos) const;

  /** @brief Get the sphere center. */
  const vector3 &getCenter() const { return center_; }

  /** @brief Get the sphere radius. */
  float getRadius() const { return radius_; }

private:
  vector3 center_;  /**< Center of the bounding sphere. */
  float radius_;    /**< Radius of the bounding sphere. */
};


// ============================================================================
// Implementation
// ============================================================================

#include <cmath>

inline SphericalContainer::SphericalContainer(const vector3 &center, float radius)
    : center_(center), radius_(radius) {}

inline bool SphericalContainer::contains(const vector3 &pos) const {
  vector3 diff = pos - center_;
  return diff.lengthSqr() <= radius_ * radius_;
}

inline vector3 SphericalContainer::clampToSphere(const vector3 &pos) const {
  vector3 diff = pos - center_;
  float dist_sq = diff.lengthSqr();

  if (dist_sq <= radius_ * radius_) {
    // Already inside the sphere
    return pos;
  }

  float dist = std::sqrt(dist_sq);
  if (dist < 1e-12f) {
    // Point is at the center; pick an arbitrary surface point
    return vector3(center_.x + radius_, center_.y, center_.z);
  }

  // Project onto sphere surface: center + radius * unit_direction
  return center_ + diff * (radius_ / dist);
}

inline float SphericalContainer::boundaryPenalty(const vector3 &pos) const {
  vector3 diff = pos - center_;
  float dist_sq = diff.lengthSqr();
  float r_sq = radius_ * radius_;

  if (dist_sq <= r_sq) {
    return 0.0f;
  }

  // Quadratic penalty: (distance_from_surface)^2
  float dist = std::sqrt(dist_sq);
  float overshoot = dist - radius_;
  return overshoot * overshoot;
}

#endif /* SPHERICALCONTAINER_H_ */
