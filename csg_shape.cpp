/*************************************************************************/
/*  csg_shape.cpp                                                        */
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

#include "csg_shape.h"

#include "core/core_string_names.h"
#include "core/math/geometry_2d.h"
#include "core/math/math_defs.h"
#include "csg.h"
#include "scene/scene_string_names.h"
#ifdef TOOLS_ENABLED
#include "editor/editor_node.h"
#endif

#include "thirdparty/manifold/src/manifold/include/manifold.h"

void CSGManifold3D::set_use_collision(bool p_enable) {
	if (use_collision == p_enable) {
		return;
	}

	use_collision = p_enable;

	if (!is_inside_tree() || !is_root_shape()) {
		return;
	}

	if (use_collision) {
		root_collision_shape.instantiate();
		root_collision_instance = PhysicsServer3D::get_singleton()->body_create();
		PhysicsServer3D::get_singleton()->body_set_mode(root_collision_instance, PhysicsServer3D::BODY_MODE_STATIC);
		PhysicsServer3D::get_singleton()->body_set_state(root_collision_instance, PhysicsServer3D::BODY_STATE_TRANSFORM, get_global_transform());
		PhysicsServer3D::get_singleton()->body_add_shape(root_collision_instance, root_collision_shape->get_rid());
		PhysicsServer3D::get_singleton()->body_set_space(root_collision_instance, get_world_3d()->get_space());
		PhysicsServer3D::get_singleton()->body_attach_object_instance_id(root_collision_instance, get_instance_id());
		set_collision_layer(collision_layer);
		set_collision_mask(collision_mask);
		_make_dirty(); //force update
	} else {
		PhysicsServer3D::get_singleton()->free(root_collision_instance);
		root_collision_instance = RID();
		root_collision_shape.unref();
	}
	notify_property_list_changed();
}

bool CSGManifold3D::is_using_collision() const {
	return use_collision;
}

void CSGManifold3D::set_collision_layer(uint32_t p_layer) {
	collision_layer = p_layer;
	if (root_collision_instance.is_valid()) {
		PhysicsServer3D::get_singleton()->body_set_collision_layer(root_collision_instance, p_layer);
	}
}

uint32_t CSGManifold3D::get_collision_layer() const {
	return collision_layer;
}

void CSGManifold3D::set_collision_mask(uint32_t p_mask) {
	collision_mask = p_mask;
	if (root_collision_instance.is_valid()) {
		PhysicsServer3D::get_singleton()->body_set_collision_mask(root_collision_instance, p_mask);
	}
}

uint32_t CSGManifold3D::get_collision_mask() const {
	return collision_mask;
}

void CSGManifold3D::set_collision_layer_value(int p_layer_number, bool p_value) {
	ERR_FAIL_COND_MSG(p_layer_number < 1, "Collision layer number must be between 1 and 32 inclusive.");
	ERR_FAIL_COND_MSG(p_layer_number > 32, "Collision layer number must be between 1 and 32 inclusive.");
	uint32_t collision_layer = get_collision_layer();
	if (p_value) {
		collision_layer |= 1 << (p_layer_number - 1);
	} else {
		collision_layer &= ~(1 << (p_layer_number - 1));
	}
	set_collision_layer(collision_layer);
}

bool CSGManifold3D::get_collision_layer_value(int p_layer_number) const {
	ERR_FAIL_COND_V_MSG(p_layer_number < 1, false, "Collision layer number must be between 1 and 32 inclusive.");
	ERR_FAIL_COND_V_MSG(p_layer_number > 32, false, "Collision layer number must be between 1 and 32 inclusive.");
	return get_collision_layer() & (1 << (p_layer_number - 1));
}

void CSGManifold3D::set_collision_mask_value(int p_layer_number, bool p_value) {
	ERR_FAIL_COND_MSG(p_layer_number < 1, "Collision layer number must be between 1 and 32 inclusive.");
	ERR_FAIL_COND_MSG(p_layer_number > 32, "Collision layer number must be between 1 and 32 inclusive.");
	uint32_t mask = get_collision_mask();
	if (p_value) {
		mask |= 1 << (p_layer_number - 1);
	} else {
		mask &= ~(1 << (p_layer_number - 1));
	}
	set_collision_mask(mask);
}

bool CSGManifold3D::get_collision_mask_value(int p_layer_number) const {
	ERR_FAIL_COND_V_MSG(p_layer_number < 1, false, "Collision layer number must be between 1 and 32 inclusive.");
	ERR_FAIL_COND_V_MSG(p_layer_number > 32, false, "Collision layer number must be between 1 and 32 inclusive.");
	return get_collision_mask() & (1 << (p_layer_number - 1));
}

bool CSGManifold3D::is_root_shape() const {
	return !parent_shape;
}

void CSGManifold3D::set_snap(float p_snap) {
	snap = p_snap;
}

float CSGManifold3D::get_snap() const {
	return snap;
}

void CSGManifold3D::_make_dirty(bool p_parent_removing) {
	if ((p_parent_removing || is_root_shape()) && !dirty) {
		call_deferred(SNAME("_update_shape")); // Must be deferred; otherwise, is_root_shape() will use the previous parent
	}

	if (!is_root_shape()) {
		parent_shape->_make_dirty();
	} else if (!dirty) {
		call_deferred(SNAME("_update_shape"));
	}

	dirty = true;
}

CSGManifoldBrush *CSGManifold3D::_get_brush() {
	if (dirty) {
		if (brush) {
			memdelete(brush);
		}
		brush = nullptr;

		CSGManifoldBrush *n = _build_brush();

		for (int i = 0; i < get_child_count(); i++) {
			CSGManifold3D *child = Object::cast_to<CSGManifold3D>(get_child(i));
			if (!child) {
				continue;
			}
			if (!child->is_visible()) {
				continue;
			}

			CSGManifoldBrush *n2 = child->_get_brush();
			if (!n2) {
				continue;
			}
			if (!n) {
				n = memnew(CSGManifoldBrush);

				n->copy_from(*n2, child->get_transform());

			} else {
				CSGManifoldBrush *nn = memnew(CSGManifoldBrush);
				CSGManifoldBrush *nn2 = memnew(CSGManifoldBrush);
				nn2->copy_from(*n2, child->get_transform());

				CSGManifoldBrushOperation bop;

				switch (child->get_operation()) {
					case CSGManifold3D::OPERATION_UNION:
						bop.merge_brushes(CSGManifoldBrushOperation::OPERATION_UNION, *n, *nn2, *nn, snap);
						break;
					case CSGManifold3D::OPERATION_INTERSECTION:
						bop.merge_brushes(CSGManifoldBrushOperation::OPERATION_INTERSECTION, *n, *nn2, *nn, snap);
						break;
					case CSGManifold3D::OPERATION_SUBTRACTION:
						bop.merge_brushes(CSGManifoldBrushOperation::OPERATION_SUBTRACTION, *n, *nn2, *nn, snap);
						break;
				}
				memdelete(n);
				memdelete(nn2);
				n = nn;
			}
		}

		if (n) {
			AABB aabb;
			for (int i = 0; i < n->faces.size(); i++) {
				for (int j = 0; j < 3; j++) {
					if (i == 0 && j == 0) {
						aabb.position = n->faces[i].vertices[j];
					} else {
						aabb.expand_to(n->faces[i].vertices[j]);
					}
				}
			}
			node_aabb = aabb;
		} else {
			node_aabb = AABB();
		}

		brush = n;

		dirty = false;
	}

	return brush;
}

int CSGManifold3D::mikktGetNumFaces(const SMikkTSpaceContext *pContext) {
	ShapeUpdateSurface &surface = *((ShapeUpdateSurface *)pContext->m_pUserData);

	return surface.vertices.size() / 3;
}

int CSGManifold3D::mikktGetNumVerticesOfFace(const SMikkTSpaceContext *pContext, const int iFace) {
	// always 3
	return 3;
}

void CSGManifold3D::mikktGetPosition(const SMikkTSpaceContext *pContext, float fvPosOut[], const int iFace, const int iVert) {
	ShapeUpdateSurface &surface = *((ShapeUpdateSurface *)pContext->m_pUserData);

	Vector3 v = surface.verticesw[iFace * 3 + iVert];
	fvPosOut[0] = v.x;
	fvPosOut[1] = v.y;
	fvPosOut[2] = v.z;
}

void CSGManifold3D::mikktGetNormal(const SMikkTSpaceContext *pContext, float fvNormOut[], const int iFace, const int iVert) {
	ShapeUpdateSurface &surface = *((ShapeUpdateSurface *)pContext->m_pUserData);

	Vector3 n = surface.normalsw[iFace * 3 + iVert];
	fvNormOut[0] = n.x;
	fvNormOut[1] = n.y;
	fvNormOut[2] = n.z;
}

void CSGManifold3D::mikktGetTexCoord(const SMikkTSpaceContext *pContext, float fvTexcOut[], const int iFace, const int iVert) {
	ShapeUpdateSurface &surface = *((ShapeUpdateSurface *)pContext->m_pUserData);

	Vector2 t = surface.uvsw[iFace * 3 + iVert];
	fvTexcOut[0] = t.x;
	fvTexcOut[1] = t.y;
}

void CSGManifold3D::mikktSetTSpaceDefault(const SMikkTSpaceContext *pContext, const float fvTangent[], const float fvBiTangent[], const float fMagS, const float fMagT,
		const tbool bIsOrientationPreserving, const int iFace, const int iVert) {
	ShapeUpdateSurface &surface = *((ShapeUpdateSurface *)pContext->m_pUserData);

	int i = iFace * 3 + iVert;
	Vector3 normal = surface.normalsw[i];
	Vector3 tangent = Vector3(fvTangent[0], fvTangent[1], fvTangent[2]);
	Vector3 bitangent = Vector3(-fvBiTangent[0], -fvBiTangent[1], -fvBiTangent[2]); // for some reason these are reversed, something with the coordinate system in Godot
	float d = bitangent.dot(normal.cross(tangent));

	i *= 4;
	surface.tansw[i++] = tangent.x;
	surface.tansw[i++] = tangent.y;
	surface.tansw[i++] = tangent.z;
	surface.tansw[i++] = d < 0 ? -1 : 1;
}

void CSGManifold3D::_update_shape() {
	if (!is_root_shape()) {
		return;
	}

	set_base(RID());
	root_mesh.unref(); //byebye root mesh

	CSGManifoldBrush *n = _get_brush();
	ERR_FAIL_COND_MSG(!n, "Cannot get CSGManifoldBrush.");

	OAHashMap<Vector3, Vector3> vec_map;

	Vector<int> face_count;
	face_count.resize(n->materials.size() + 1);
	for (int i = 0; i < face_count.size(); i++) {
		face_count.write[i] = 0;
	}

	for (int i = 0; i < n->faces.size(); i++) {
		int mat = n->faces[i].material;
		ERR_CONTINUE(mat < -1 || mat >= face_count.size());
		int idx = mat == -1 ? face_count.size() - 1 : mat;

		if (n->faces[i].smooth) {
			Plane p(n->faces[i].vertices[0], n->faces[i].vertices[1], n->faces[i].vertices[2]);

			for (int j = 0; j < 3; j++) {
				Vector3 v = n->faces[i].vertices[j];
				Vector3 add;
				if (vec_map.lookup(v, add)) {
					add += p.normal;
				} else {
					add = p.normal;
				}
				vec_map.set(v, add);
			}
		}

		face_count.write[idx]++;
	}

	Vector<ShapeUpdateSurface> surfaces;

	surfaces.resize(face_count.size());

	//create arrays
	for (int i = 0; i < surfaces.size(); i++) {
		surfaces.write[i].vertices.resize(face_count[i] * 3);
		surfaces.write[i].normals.resize(face_count[i] * 3);
		surfaces.write[i].uvs.resize(face_count[i] * 3);
		if (calculate_tangents) {
			surfaces.write[i].tans.resize(face_count[i] * 3 * 4);
		}
		surfaces.write[i].last_added = 0;

		if (i != surfaces.size() - 1) {
			surfaces.write[i].material = n->materials[i];
		}

		surfaces.write[i].verticesw = surfaces.write[i].vertices.ptrw();
		surfaces.write[i].normalsw = surfaces.write[i].normals.ptrw();
		surfaces.write[i].uvsw = surfaces.write[i].uvs.ptrw();
		if (calculate_tangents) {
			surfaces.write[i].tansw = surfaces.write[i].tans.ptrw();
		}
	}

	//fill arrays
	{
		for (int i = 0; i < n->faces.size(); i++) {
			int order[3] = { 0, 1, 2 };

			if (n->faces[i].invert) {
				SWAP(order[1], order[2]);
			}

			int mat = n->faces[i].material;
			ERR_CONTINUE(mat < -1 || mat >= face_count.size());
			int idx = mat == -1 ? face_count.size() - 1 : mat;

			int last = surfaces[idx].last_added;

			Plane p(n->faces[i].vertices[0], n->faces[i].vertices[1], n->faces[i].vertices[2]);

			for (int j = 0; j < 3; j++) {
				Vector3 v = n->faces[i].vertices[j];

				Vector3 normal = p.normal;

				if (n->faces[i].smooth && vec_map.lookup(v, normal)) {
					normal.normalize();
				}

				if (n->faces[i].invert) {
					normal = -normal;
				}

				int k = last + order[j];
				surfaces[idx].verticesw[k] = v;
				surfaces[idx].uvsw[k] = n->faces[i].uvs[j];
				surfaces[idx].normalsw[k] = normal;

				if (calculate_tangents) {
					// zero out our tangents for now
					k *= 4;
					surfaces[idx].tansw[k++] = 0.0;
					surfaces[idx].tansw[k++] = 0.0;
					surfaces[idx].tansw[k++] = 0.0;
					surfaces[idx].tansw[k++] = 0.0;
				}
			}

			surfaces.write[idx].last_added += 3;
		}
	}

	root_mesh.instantiate();
	//create surfaces

	for (int i = 0; i < surfaces.size(); i++) {
		// calculate tangents for this surface
		bool have_tangents = calculate_tangents;
		if (have_tangents) {
			SMikkTSpaceInterface mkif;
			mkif.m_getNormal = mikktGetNormal;
			mkif.m_getNumFaces = mikktGetNumFaces;
			mkif.m_getNumVerticesOfFace = mikktGetNumVerticesOfFace;
			mkif.m_getPosition = mikktGetPosition;
			mkif.m_getTexCoord = mikktGetTexCoord;
			mkif.m_setTSpace = mikktSetTSpaceDefault;
			mkif.m_setTSpaceBasic = nullptr;

			SMikkTSpaceContext msc;
			msc.m_pInterface = &mkif;
			msc.m_pUserData = &surfaces.write[i];
			have_tangents = genTangSpaceDefault(&msc);
		}

		if (surfaces[i].last_added == 0) {
			continue;
		}

		// and convert to surface array
		Array array;
		array.resize(Mesh::ARRAY_MAX);

		array[Mesh::ARRAY_VERTEX] = surfaces[i].vertices;
		array[Mesh::ARRAY_NORMAL] = surfaces[i].normals;
		array[Mesh::ARRAY_TEX_UV] = surfaces[i].uvs;
		if (have_tangents) {
			array[Mesh::ARRAY_TANGENT] = surfaces[i].tans;
		}

		int idx = root_mesh->get_surface_count();
		root_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, array);
		root_mesh->surface_set_material(idx, surfaces[i].material);
	}

	set_base(root_mesh->get_rid());

	_update_collision_faces();
}

void CSGManifold3D::_update_collision_faces() {
	if (use_collision && is_root_shape() && root_collision_shape.is_valid()) {
		CSGManifoldBrush *n = _get_brush();
		ERR_FAIL_COND_MSG(!n, "Cannot get CSGManifoldBrush.");
		Vector<Vector3> physics_faces;
		physics_faces.resize(n->faces.size() * 3);
		Vector3 *physicsw = physics_faces.ptrw();

		for (int i = 0; i < n->faces.size(); i++) {
			int order[3] = { 0, 1, 2 };

			if (n->faces[i].invert) {
				SWAP(order[1], order[2]);
			}

			physicsw[i * 3 + 0] = n->faces[i].vertices[order[0]];
			physicsw[i * 3 + 1] = n->faces[i].vertices[order[1]];
			physicsw[i * 3 + 2] = n->faces[i].vertices[order[2]];
		}

		root_collision_shape->set_faces(physics_faces);
	}
}

AABB CSGManifold3D::get_aabb() const {
	return node_aabb;
}

Vector<Vector3> CSGManifold3D::get_brush_faces() {
	ERR_FAIL_COND_V(!is_inside_tree(), Vector<Vector3>());
	CSGManifoldBrush *b = _get_brush();
	if (!b) {
		return Vector<Vector3>();
	}

	Vector<Vector3> faces;
	int fc = b->faces.size();
	faces.resize(fc * 3);
	{
		Vector3 *w = faces.ptrw();
		for (int i = 0; i < fc; i++) {
			w[i * 3 + 0] = b->faces[i].vertices[0];
			w[i * 3 + 1] = b->faces[i].vertices[1];
			w[i * 3 + 2] = b->faces[i].vertices[2];
		}
	}

	return faces;
}

void CSGManifold3D::_notification(int p_what) {
	switch (p_what) {
		case NOTIFICATION_PARENTED: {
			Node *parentn = get_parent();
			if (parentn) {
				parent_shape = Object::cast_to<CSGManifold3D>(parentn);
				if (parent_shape) {
					set_base(RID());
					root_mesh.unref();
				}
			}
			if (!brush || parent_shape) {
				// Update this node if uninitialized, or both this node and its new parent if it gets added to another CSG shape
				_make_dirty();
			}
			last_visible = is_visible();
		} break;

		case NOTIFICATION_UNPARENTED: {
			if (!is_root_shape()) {
				// Update this node and its previous parent only if it's currently being removed from another CSG shape
				_make_dirty(true); // Must be forced since is_root_shape() uses the previous parent
			}
			parent_shape = nullptr;
		} break;

		case NOTIFICATION_VISIBILITY_CHANGED: {
			if (!is_root_shape() && last_visible != is_visible()) {
				// Update this node's parent only if its own visibility has changed, not the visibility of parent nodes
				parent_shape->_make_dirty();
			}
			last_visible = is_visible();
		} break;

		case NOTIFICATION_LOCAL_TRANSFORM_CHANGED: {
			if (!is_root_shape()) {
				// Update this node's parent only if its own transformation has changed, not the transformation of parent nodes
				parent_shape->_make_dirty();
			}
		} break;
		case NOTIFICATION_ENTER_TREE: {
			if (use_collision && is_root_shape()) {
				root_collision_shape.instantiate();
				root_collision_instance = PhysicsServer3D::get_singleton()->body_create();
				PhysicsServer3D::get_singleton()->body_set_mode(root_collision_instance, PhysicsServer3D::BODY_MODE_STATIC);
				PhysicsServer3D::get_singleton()->body_set_state(root_collision_instance, PhysicsServer3D::BODY_STATE_TRANSFORM, get_global_transform());
				PhysicsServer3D::get_singleton()->body_add_shape(root_collision_instance, root_collision_shape->get_rid());
				PhysicsServer3D::get_singleton()->body_set_space(root_collision_instance, get_world_3d()->get_space());
				PhysicsServer3D::get_singleton()->body_attach_object_instance_id(root_collision_instance, get_instance_id());
				set_collision_layer(collision_layer);
				set_collision_mask(collision_mask);
				_update_collision_faces();
				_mesh_changed();
			}
		} break;

		case NOTIFICATION_EXIT_TREE: {
			if (use_collision && is_root_shape() && root_collision_instance.is_valid()) {
				PhysicsServer3D::get_singleton()->free(root_collision_instance);
				root_collision_instance = RID();
				root_collision_shape.unref();
			}
		} break;

		case NOTIFICATION_TRANSFORM_CHANGED: {
			if (use_collision && is_root_shape() && root_collision_instance.is_valid()) {
				PhysicsServer3D::get_singleton()->body_set_state(root_collision_instance, PhysicsServer3D::BODY_STATE_TRANSFORM, get_global_transform());
			}
		} break;
	}
}

void CSGManifold3D::set_operation(Operation p_operation) {
	operation = p_operation;
	_make_dirty();
	update_gizmos();
}

CSGManifold3D::Operation CSGManifold3D::get_operation() const {
	return operation;
}

void CSGManifold3D::set_calculate_tangents(bool p_calculate_tangents) {
	calculate_tangents = p_calculate_tangents;
	_make_dirty();
}

bool CSGManifold3D::is_calculating_tangents() const {
	return calculate_tangents;
}

void CSGManifold3D::_validate_property(PropertyInfo &property) const {
	bool is_collision_prefixed = property.name.begins_with("collision_");
	if ((is_collision_prefixed || property.name.begins_with("use_collision")) && is_inside_tree() && !is_root_shape()) {
		//hide collision if not root
		property.usage = PROPERTY_USAGE_NO_EDITOR;
	} else if (is_collision_prefixed && !bool(get("use_collision"))) {
		property.usage = PROPERTY_USAGE_NO_EDITOR | PROPERTY_USAGE_INTERNAL;
	}
	GeometryInstance3D::_validate_property(property);
}

Array CSGManifold3D::get_meshes() const {
	if (root_mesh.is_valid()) {
		Array arr;
		arr.resize(2);
		arr[0] = Transform3D();
		arr[1] = root_mesh;
		return arr;
	}

	return Array();
}

void CSGManifold3D::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_flip_faces", "flip_faces"), &CSGManifold3D::set_flip_faces);
	ClassDB::bind_method(D_METHOD("get_flip_faces"), &CSGManifold3D::get_flip_faces);

	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "flip_faces"), "set_flip_faces", "get_flip_faces");

	ClassDB::bind_method(D_METHOD("_update_shape"), &CSGManifold3D::_update_shape);
	ClassDB::bind_method(D_METHOD("is_root_shape"), &CSGManifold3D::is_root_shape);

	ClassDB::bind_method(D_METHOD("set_operation", "operation"), &CSGManifold3D::set_operation);
	ClassDB::bind_method(D_METHOD("get_operation"), &CSGManifold3D::get_operation);

	ClassDB::bind_method(D_METHOD("set_snap", "snap"), &CSGManifold3D::set_snap);
	ClassDB::bind_method(D_METHOD("get_snap"), &CSGManifold3D::get_snap);

	ClassDB::bind_method(D_METHOD("set_use_collision", "operation"), &CSGManifold3D::set_use_collision);
	ClassDB::bind_method(D_METHOD("is_using_collision"), &CSGManifold3D::is_using_collision);

	ClassDB::bind_method(D_METHOD("set_collision_layer", "layer"), &CSGManifold3D::set_collision_layer);
	ClassDB::bind_method(D_METHOD("get_collision_layer"), &CSGManifold3D::get_collision_layer);

	ClassDB::bind_method(D_METHOD("set_collision_mask", "mask"), &CSGManifold3D::set_collision_mask);
	ClassDB::bind_method(D_METHOD("get_collision_mask"), &CSGManifold3D::get_collision_mask);

	ClassDB::bind_method(D_METHOD("set_collision_mask_value", "layer_number", "value"), &CSGManifold3D::set_collision_mask_value);
	ClassDB::bind_method(D_METHOD("get_collision_mask_value", "layer_number"), &CSGManifold3D::get_collision_mask_value);

	ClassDB::bind_method(D_METHOD("set_collision_layer_value", "layer_number", "value"), &CSGManifold3D::set_collision_layer_value);
	ClassDB::bind_method(D_METHOD("get_collision_layer_value", "layer_number"), &CSGManifold3D::get_collision_layer_value);

	ClassDB::bind_method(D_METHOD("set_calculate_tangents", "enabled"), &CSGManifold3D::set_calculate_tangents);
	ClassDB::bind_method(D_METHOD("is_calculating_tangents"), &CSGManifold3D::is_calculating_tangents);

	ClassDB::bind_method(D_METHOD("get_meshes"), &CSGManifold3D::get_meshes);

	ADD_PROPERTY(PropertyInfo(Variant::INT, "operation", PROPERTY_HINT_ENUM, "Union,Intersection,Subtraction"), "set_operation", "get_operation");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "snap", PROPERTY_HINT_RANGE, "0.0001,1,0.001,suffix:m"), "set_snap", "get_snap");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "calculate_tangents"), "set_calculate_tangents", "is_calculating_tangents");

	ADD_GROUP("Collision", "collision_");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "use_collision"), "set_use_collision", "is_using_collision");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "collision_layer", PROPERTY_HINT_LAYERS_3D_PHYSICS), "set_collision_layer", "get_collision_layer");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "collision_mask", PROPERTY_HINT_LAYERS_3D_PHYSICS), "set_collision_mask", "get_collision_mask");

	BIND_ENUM_CONSTANT(OPERATION_UNION);
	BIND_ENUM_CONSTANT(OPERATION_INTERSECTION);
	BIND_ENUM_CONSTANT(OPERATION_SUBTRACTION);
}

/////////////////////

CSGManifoldBrush *CSGManifold3D::_create_brush_from_arrays(const Vector<Vector3> &p_vertices, const Vector<Vector2> &p_uv, const Vector<bool> &p_smooth, const Vector<Ref<Material>> &p_materials) {
	CSGManifoldBrush *brush = memnew(CSGManifoldBrush);

	Vector<bool> invert;
	invert.resize(p_vertices.size() / 3);
	{
		int ic = invert.size();
		bool *w = invert.ptrw();
		for (int i = 0; i < ic; i++) {
			w[i] = flip_faces;
		}
	}
	brush->build_from_faces(p_vertices, p_uv, p_smooth, p_materials, invert);

	return brush;
}

void CSGManifold3D::set_flip_faces(bool p_invert) {
	if (flip_faces == p_invert) {
		return;
	}

	flip_faces = p_invert;

	_make_dirty();
}

bool CSGManifold3D::get_flip_faces() {
	return flip_faces;
}

/////////////////////

CSGManifoldBrush *CSGManifold3D::_build_brush() {
	if (!mesh.is_valid()) {
		return memnew(CSGManifoldBrush);
	}

	Vector<Vector3> vertices;
	Vector<bool> smooth;
	Vector<Ref<Material>> materials;
	Vector<Vector2> uvs;

	for (int i = 0; i < get_mesh()->get_surface_count(); i++) {
		if (get_mesh()->surface_get_primitive_type(i) != Mesh::PRIMITIVE_TRIANGLES) {
			continue;
		}
		Ref<Material> material = get_active_material(i);

		Array arrays = get_mesh()->surface_get_arrays(i);

		if (arrays.size() == 0) {
			_make_dirty();
			ERR_FAIL_COND_V(arrays.size() == 0, memnew(CSGManifoldBrush));
		}

		Vector<Vector3> avertices = arrays[Mesh::ARRAY_VERTEX];
		if (avertices.size() == 0) {
			continue;
		}

		const Vector3 *vr = avertices.ptr();

		Vector<Vector3> anormals = arrays[Mesh::ARRAY_NORMAL];
		const Vector3 *nr = nullptr;
		if (anormals.size()) {
			nr = anormals.ptr();
		}

		Vector<Vector2> auvs = arrays[Mesh::ARRAY_TEX_UV];
		const Vector2 *uvr = nullptr;
		if (auvs.size()) {
			uvr = auvs.ptr();
		}

		Ref<Material> mat;
		if (material.is_valid()) {
			mat = material;
		} else {
			mat = mesh->surface_get_material(i);
		}

		Vector<int> aindices = arrays[Mesh::ARRAY_INDEX];
		if (aindices.size()) {
			int as = vertices.size();
			int is = aindices.size();

			vertices.resize(as + is);
			smooth.resize((as + is) / 3);
			materials.resize((as + is) / 3);
			uvs.resize(as + is);

			Vector3 *vw = vertices.ptrw();
			bool *sw = smooth.ptrw();
			Vector2 *uvw = uvs.ptrw();
			Ref<Material> *mw = materials.ptrw();

			const int *ir = aindices.ptr();

			for (int j = 0; j < is; j += 3) {
				Vector3 vertex[3];
				Vector3 normal[3];
				Vector2 uv[3];

				for (int k = 0; k < 3; k++) {
					int idx = ir[j + k];
					vertex[k] = vr[idx];
					if (nr) {
						normal[k] = nr[idx];
					}
					if (uvr) {
						uv[k] = uvr[idx];
					}
				}

				bool flat = normal[0].is_equal_approx(normal[1]) && normal[0].is_equal_approx(normal[2]);

				vw[as + j + 0] = vertex[0];
				vw[as + j + 1] = vertex[1];
				vw[as + j + 2] = vertex[2];

				uvw[as + j + 0] = uv[0];
				uvw[as + j + 1] = uv[1];
				uvw[as + j + 2] = uv[2];

				sw[(as + j) / 3] = !flat;
				mw[(as + j) / 3] = mat;
			}
		} else {
			int as = vertices.size();
			int is = avertices.size();

			vertices.resize(as + is);
			smooth.resize((as + is) / 3);
			uvs.resize(as + is);
			materials.resize((as + is) / 3);

			Vector3 *vw = vertices.ptrw();
			bool *sw = smooth.ptrw();
			Vector2 *uvw = uvs.ptrw();
			Ref<Material> *mw = materials.ptrw();

			for (int j = 0; j < is; j += 3) {
				Vector3 vertex[3];
				Vector3 normal[3];
				Vector2 uv[3];

				for (int k = 0; k < 3; k++) {
					vertex[k] = vr[j + k];
					if (nr) {
						normal[k] = nr[j + k];
					}
					if (uvr) {
						uv[k] = uvr[j + k];
					}
				}

				bool flat = normal[0].is_equal_approx(normal[1]) && normal[0].is_equal_approx(normal[2]);

				vw[as + j + 0] = vertex[0];
				vw[as + j + 1] = vertex[1];
				vw[as + j + 2] = vertex[2];

				uvw[as + j + 0] = uv[0];
				uvw[as + j + 1] = uv[1];
				uvw[as + j + 2] = uv[2];

				sw[(as + j) / 3] = !flat;
				mw[(as + j) / 3] = mat;
			}
		}
	}

	if (vertices.size() == 0) {
		return memnew(CSGManifoldBrush);
	}

	return _create_brush_from_arrays(vertices, uvs, smooth, materials);
}

void CSGManifold3D::_mesh_changed() {
	_make_dirty();
	update_gizmos();
	if (Engine::get_singleton()->is_editor_hint()) {
		update_configuration_warnings();
	}
}

TypedArray<String> CSGManifold3D::get_configuration_warnings() const {
	TypedArray<String> warnings = Node::get_configuration_warnings();
	HashMap<int64_t, std::vector<float>> mesh_id_properties;
	HashMap<int64_t, HashMap<int32_t, Ref<Material>>> mesh_materials;
	HashMap<int64_t, int> mesh_face_count;
	manifold::Manifold manifold;
	warnings.push_back(RTR("CSGManifold3D is not a 3d manifold, cannot run merge, substraction or intersection operations."));
	if (!brush) {
		return warnings;
	}
	CSGManifoldBrush::pack_manifold(brush, manifold, mesh_id_properties, mesh_materials, mesh_face_count, CMP_EPSILON);
	if (!manifold.IsManifold()) {
		return warnings;
	}
	return Node::get_configuration_warnings();
}
