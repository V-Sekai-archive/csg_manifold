/*************************************************************************/
/*  csg_shape.h                                                          */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2022 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2022 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#ifndef CSG_MANIFOLD_SHAPE_H
#define CSG_MANIFOLD_SHAPE_H

#include "csg.h"
#include "scene/3d/mesh_instance_3d.h"
#include "scene/3d/path_3d.h"
#include "scene/3d/visual_instance_3d.h"
#include "scene/resources/concave_polygon_shape_3d.h"
#include "thirdparty/misc/mikktspace.h"

class CSGManifold3D : public MeshInstance3D {
	GDCLASS(CSGManifold3D, MeshInstance3D);

public:
	enum Operation {
		OPERATION_UNION,
		OPERATION_INTERSECTION,
		OPERATION_SUBTRACTION,

	};

private:
	Operation operation = OPERATION_UNION;
	CSGManifold3D *parent_shape = nullptr;

	CSGManifoldBrush *brush = nullptr;

	AABB node_aabb;

	bool dirty = false;
	bool last_visible = false;
	float snap = 0.001;

	bool use_collision = false;
	uint32_t collision_layer = 1;
	uint32_t collision_mask = 1;
	Ref<ConcavePolygonShape3D> root_collision_shape;
	RID root_collision_instance;

	bool calculate_tangents = true;

	Ref<ArrayMesh> root_mesh;

	struct Vector3Hasher {
		_ALWAYS_INLINE_ uint32_t hash(const Vector3 &p_vec3) const {
			uint32_t h = hash_murmur3_one_float(p_vec3.x);
			h = hash_murmur3_one_float(p_vec3.y, h);
			h = hash_murmur3_one_float(p_vec3.z, h);
			return h;
		}
	};

	struct ShapeUpdateSurface {
		Vector<Vector3> vertices;
		Vector<Vector3> normals;
		Vector<Vector2> uvs;
		Vector<real_t> tans;
		Ref<Material> material;
		int last_added = 0;

		Vector3 *verticesw = nullptr;
		Vector3 *normalsw = nullptr;
		Vector2 *uvsw = nullptr;
		real_t *tansw = nullptr;
	};

	//mikktspace callbacks
	static int mikktGetNumFaces(const SMikkTSpaceContext *pContext);
	static int mikktGetNumVerticesOfFace(const SMikkTSpaceContext *pContext, const int iFace);
	static void mikktGetPosition(const SMikkTSpaceContext *pContext, float fvPosOut[], const int iFace, const int iVert);
	static void mikktGetNormal(const SMikkTSpaceContext *pContext, float fvNormOut[], const int iFace, const int iVert);
	static void mikktGetTexCoord(const SMikkTSpaceContext *pContext, float fvTexcOut[], const int iFace, const int iVert);
	static void mikktSetTSpaceDefault(const SMikkTSpaceContext *pContext, const float fvTangent[], const float fvBiTangent[], const float fMagS, const float fMagT,
			const tbool bIsOrientationPreserving, const int iFace, const int iVert);

	void _update_shape();
	void _update_collision_faces();

protected:
	void _notification(int p_what);
	virtual CSGManifoldBrush *_build_brush();
	void _make_dirty(bool p_parent_removing = false);

	static void _bind_methods();

	CSGManifoldBrush *_get_brush();

	virtual void _validate_property(PropertyInfo &property) const override;

public:
	TypedArray<String> get_configuration_warnings() const override;

	Array get_meshes() const;

	void set_operation(Operation p_operation);
	Operation get_operation() const;

	virtual Vector<Vector3> get_brush_faces();

	virtual AABB get_aabb() const override;

	void set_use_collision(bool p_enable);
	bool is_using_collision() const;

	void set_collision_layer(uint32_t p_layer);
	uint32_t get_collision_layer() const;

	void set_collision_mask(uint32_t p_mask);
	uint32_t get_collision_mask() const;

	void set_collision_layer_value(int p_layer_number, bool p_value);
	bool get_collision_layer_value(int p_layer_number) const;

	void set_collision_mask_value(int p_layer_number, bool p_value);
	bool get_collision_mask_value(int p_layer_number) const;

	void set_snap(float p_snap);
	float get_snap() const;

	void set_calculate_tangents(bool p_calculate_tangents);
	bool is_calculating_tangents() const;

	bool is_root_shape() const;

protected:
	bool flip_faces = false;
	CSGManifoldBrush *_create_brush_from_arrays(const Vector<Vector3> &p_vertices, const Vector<Vector2> &p_uv, const Vector<bool> &p_smooth, const Vector<Ref<Material>> &p_materials);

public:
	void set_flip_faces(bool p_invert);
	bool get_flip_faces();

	void _mesh_changed();

public:
	CSGManifold3D() {
		set_notify_local_transform(true);
	}

	~CSGManifold3D() {
		if (brush) {
			memdelete(brush);
			brush = nullptr;
		}
	}
};

VARIANT_ENUM_CAST(CSGManifold3D::Operation)
#endif // CSG_MANIFOLD_SHAPE_H
