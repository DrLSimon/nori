/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <nori/mesh.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class Accel : public Shape{
public:
    /// Create a new and empty BVH
    Accel() {m_meshOffset.push_back(0u);}

    /// Release all resources
    virtual ~Accel() { clear(); };

    /// Release all resources
    void clear();

    /**
     * \brief Register a triangle mesh for inclusion in the BVH.
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the BVH
    void activate() override;

    /**
     * \brief Intersect a ray against all triangle meshes registered
     * with the BVH
     *
     * Detailed information about the intersection, if any, will be
     * stored in the provided \ref Intersection data record. 
     *
     * The <tt>shadowRay</tt> parameter specifies whether this detailed
     * information is really needed. When set to \c true, the 
     * function just checks whether or not there is occlusion, but without
     * providing any more detail (i.e. \c its will not be filled with
     * contents). This is usually much faster.
     *
     * \return \c true If an intersection was found
     */
    bool rayIntersect(Ray3f &ray, Intersection &its, 
        bool shadowRay = false) const override;

    /// Return the total number of meshes registered with the BVH
    uint32_t getMeshCount() const { return (uint32_t) m_meshes.size(); }

    /// Return the total number of internally represented triangles 
    uint32_t getTriangleCount() const { return m_meshOffset.back(); }

    /// Return one of the registered meshes
    Mesh *getMesh(uint32_t idx) { return m_meshes[idx]; }
    
    /// Return one of the registered meshes (const version)
    const Mesh *getMesh(uint32_t idx) const { return m_meshes[idx]; }

    /// Return a human-readable summary of this instance
    std::string toString() const override;

    //// Return an axis-aligned bounding box containing the entire tree
    const BoundingBox3f &getBoundingBox() const override {
        return m_bbox;
    }

protected:
    /**
     * \brief Compute the mesh and triangle indices corresponding to 
     * a primitive index used by the underlying generic BVH implementation. 
     */
    uint32_t findMesh(uint32_t &idx) const {
        auto it = std::lower_bound(m_meshOffset.begin(), m_meshOffset.end(), idx+1) - 1;
        idx -= *it;
        return (uint32_t) (it - m_meshOffset.begin());
    }

    //// Return an axis-aligned bounding box containing the given triangle
    BoundingBox3f getBoundingBox(uint32_t index) const {
        uint32_t meshIdx = findMesh(index);
        return m_meshes[meshIdx]->getBoundingBox(index);
    }
    
    //// Return the centroid of the given triangle
    Point3f getCentroid(uint32_t index) const {
        uint32_t meshIdx = findMesh(index);
        return m_meshes[meshIdx]->getCentroid(index);
    }


private:
    std::vector<Mesh *> m_meshes;       ///< List of meshes 
    std::vector<uint32_t> m_meshOffset; ///< Index of the first triangle for each shape
    BoundingBox3f m_bbox;               ///< Bounding box of the entire BVH
};

NORI_NAMESPACE_END

